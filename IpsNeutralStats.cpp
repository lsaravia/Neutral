#include "IpsNeutral.h"
#include "mf.h"
#include "hk.h" 
#include <fstream>      

using namespace std;


void IPSNeutral::CalcDensity(simplmat <double> &den)
{
	if(den.getRows()!=NumSpecies)
		den.resize(NumSpecies);
	
	den.fill(0.0);
	
	// Calculate species' densities
	//
	for(int i=0; i<DimY; i++)
		for(int j=0;j<DimX;j++)
			{
			int a = C(j,i).Elem();
			if( a>0 )
				den( a-1 )++;
			}
}

int IPSNeutral::Convert(simplmat <double> &data)
{
	int dx,dy;
	dx = data.getRows();
	dy = data.getCols();

	if( dx!=DimX || dy!=DimY )
		data.resize(DimX, DimY, 0.0);
		
	for(dy=0;dy<DimY; dy++)
		for(dx=0;dx<DimX; dx++)
			data(dx,dy) = C(dx,dy).Specie;

	return 1;
}

int IPSNeutral::Convert(simplmat <int> &data)
{
	int dx,dy;
	dx = data.getRows();
	dy = data.getCols();

	if( dx!=DimX || dy!=DimY )
		data.resize(DimX, DimY, 0.0);
		
	for(dy=0;dy<DimY; dy++)
		for(dx=0;dx<DimX; dx++)
			data(dx,dy) = C(dx,dy).Specie;

	return 1;
}

// Converting to Biomass with using the inverse of Damuth exponent
//                M=(N/a)^(-4/3) between a range bioMin & bioMax
// 
// The Biomass spectrum is calculated dynamically according to the actual densities 
// Assuming the minimun biomass to a species with a frequency = 0.9 (very dominant)
//
// The species follow a universal size-distribution according to Giometto (2013) with
// lognormal distribution.
//
double IPSNeutral::ConvertToBio(simplmat <double> &data, simplmat <double> &den,float bioMax,float bioMin)
{
	double dar=-4.0/3.0;		 // inverse of Damuth Power exponent
	double a=0,totBio=0;  

	double sigma = 0.2223; 			// From giometto 2013
	double mu = -(sigma*sigma)/2;


	// set minimun value bioMin to a density of 90% of the total
	double maxN = DimX*DimY*0.9;
	// Calculate the constant 
	a = maxN/(pow(bioMin,(1/dar)));

	// Calculate the average biomass of each species and store it in BirthRate
	for(int i=1; i<=NumSpecies; i++ )
	{
		double bio  = pow(den(i-1)/a,dar);
	    if(bio>bioMax) 
	    	bio=bioMax;
	    else if(bio<bioMin) 
	    	bio=bioMin;
		Sp[i].BirthRate = bio;
	}

  	int dx,dy;
 	dx = data.getRows();
 	dy = data.getCols();


	if( dx!=DimX || dy!=DimY )
		data.resize(DimX, DimY, 0.0);

	for(dy=0;dy<DimY; dy++)
		for(dx=0;dx<DimX; dx++)
 		{
 
			int spc = C(dx,dy).Specie;
 			if( spc>0 )
 			{
				double bio = exp(mu+sigma*normran.dev())*Sp[spc].BirthRate;
				data(dx,dy) = bio;
				totBio += bio;
 			}
 		}
	return totBio;
	
}



// Converting to Biomass with N=aM^b between a range minBio & maxBio
// 
// The Biomass spectrum is determined by the densities in the metacommunity
// The power exponent is calculated assuming that the most frequent species has 
// the minimun biomass and the least frequent has the maximun 
//
double IPSNeutral::ConvertToBio(simplmat <double> &data, float bioMax, float bioMin)
{
	double dar=-4.0/3.0;		 //  inverse of Damuth Power exponent
	double a=0,totBio=0; 
	double sigma = 0.2223; 			// From Giometto 2013 PNAS
	double mu = -(sigma*sigma)/2;

	static bool privez=true;

	if(privez)
	{

		// set minimun value of Biomass to a density of the most abundant specie
		// assumes the most abundant is the last specie
		double maxN = Sp[NumSpecies].BirthRate;
		double minN = Sp[1].BirthRate;
		for(int i=1; i<=NumSpecies; i++ )
		{
			if( Sp[i].BirthRate > maxN)  maxN=Sp[i].BirthRate;
			if( Sp[i].BirthRate < minN)  minN=Sp[i].BirthRate;
		}

		// Calculate the exponent 
		dar = log(minN/maxN)/log(bioMax/bioMin);

		// Calculate the constant 
		a = maxN/(pow(bioMin,dar));

		// BirthRate have the density in the metacommunity I replace it with mean biomass
		for(int i=1; i<=NumSpecies; i++ )
		{
			//double denMeta = Sp[i].BirthRate;
			double bio  = pow(Sp[i].BirthRate/a,1/dar);
		    if(bio>bioMax) 
		    	bio=bioMax;
		    else if(bio<bioMin) 
		    	bio=bioMin;
			Sp[i].BirthRate = bio;
		}
		privez=false;
	}

	int dx,dy;
	dx = data.getRows();
	dy = data.getCols();

	if( dx!=DimX || dy!=DimY )
		data.resize(DimX, DimY, 0.0);
		
	for(dy=0;dy<DimY; dy++)
		for(dx=0;dx<DimX; dx++)
		{

			int spc = C(dx,dy).Specie;
			if( spc>0 )
			{
				double bio = exp(mu+sigma*normran.dev())*Sp[spc].BirthRate;
	
				data(dx,dy) = bio;
				totBio += bio;
			}
		}
	return totBio;
}

int IPSNeutral::Reordering(simplmat <double> &newdata )
{
	int i,a,maxi;
	if( newdata.getRows()!=DimX || newdata.getCols()!=DimY)
		newdata.resize(DimX,DimY,0.0);
	else
		newdata.fill(0.0);
		
	double maxDen;
	
	double * den = new double[NumSpecies];

	for(i=0; i<NumSpecies; i++)
		den[i]=0;
	
	for(i=0; i<DimY; i++)
		for(int j=0;j<DimX;j++)
			{
			a = C(j,i).Specie;
			if( a>0 )
				den[ a-1 ]++;
			}
			
	int newSpecie=0;
	while(newSpecie<NumSpecies)
	{
		maxDen=maxi=0;
		for( i=0; i<NumSpecies; i++)
			if(den[i]>maxDen)
				{
					maxi=i+1;
					maxDen=den[i];
				}
		if (maxDen==0) break;
		
		//cout << maxi << "-" << maxDen << "\t";
		newSpecie++;
		for(i=0; i<DimY; i++)
			for(int j=0;j<DimX;j++)
				{
				a = C(j,i).Specie;
				if( a==maxi )
					newdata(j,i)=newSpecie;
				}
		den[maxi-1]=0;
	}
	return(1);
}

// Uses an previously calculated vector of densities
//
int IPSNeutral::Reordering(simplmat <double> &newdata, simplmat <double> &den )
{
	int i,a,maxi;
	if( newdata.getRows()!=DimX || newdata.getCols()!=DimY)
		newdata.resize(DimX,DimY,0.0);
	else
		newdata.fill(0.0);
		
	double maxDen;
	
	int newSpecie=0;
	while(newSpecie<NumSpecies)
	{
		maxDen=maxi=0;
		for( i=0; i<NumSpecies; i++)
			if(den(i)>maxDen)
				{
					maxi=i+1;
					maxDen=den(i);
				}
		if (maxDen==0) break;
		
		//cout << maxi << "-" << maxDen << "\t";
		newSpecie++;
		for(i=0; i<DimY; i++)
			for(int j=0;j<DimX;j++)
				{
				a = C(j,i).Specie;
				if( a==maxi )
					newdata(j,i)=newSpecie;
				}
		den(maxi-1)=0;
	}
	return(1);
}


// Convierte la matriz de especies solo con las que estan en species
// 
int IPSNeutral::Convert(simplmat <double> &data, const int * species )
{
	int dx,dy,sp=0;
	dx = data.getRows();
	dy = data.getCols();

	if( dx!=DimX || dy!=DimY )
		data.resize(DimX, DimY, 0.0);
		
	for(dy=0;dy<DimY; dy++)
		for(dx=0;dx<DimX; dx++)
		{
			sp=C(dx,dy).Specie;
			if( sp>0)
				if( species[sp-1] )
					data(dx,dy) = 1;
		}
	return 1;
}

// Convierte la matriz de especies solo con la "specie"
// specie [0 ... numSpecies]
// 
int IPSNeutral::Convert(simplmat <double> &data, const int specie )
{
	int dx,dy,sp=0;
	dx = data.getRows();
	dy = data.getCols();

	if( dx!=DimX || dy!=DimY )
		data.resize(DimX, DimY, 0.0);
		
	for(dy=0;dy<DimY; dy++)
		for(dx=0;dx<DimX; dx++)
		{
			sp=C(dx,dy).Specie-1;
			if( sp==specie)
				data(dx,dy) = 1;
		}
	return 1;
}

// Add a constant to a simplmat 
// 
// 
int IPSNeutral::AddConst(simplmat <double> &data, const double aa )
{
	int dx,dy,dimx,dimy,s=0;
	dimx = data.getRows();
	dimy = data.getCols();

	for(dy=0;dy<dimy; dy++)
		for(dx=0;dx<dimx; dx++)
		{
			data(dx,dy) = data(dx,dy) + aa;
			s += data(dx,dy);
		}
	return 1;
}

int IPSNeutral::PStats(simplmat <int> &data, const char * outFile, const char * ident,const string & type)
{
    hoshen_kopelman c;
    ofstream dout;
    bool privez=false;
    dout.open( outFile, ios::in );
    if( !dout )
            privez=true;

    dout.close();
    dout.clear();
    dout.open( outFile, ios::app );
    if( !dout )
    {
            cerr << "Cannot open output file.\n";
            return 0;
    }

    // With output=max multiClusters returns 2 pairs:
    // the 1st with maxCluster the 2nd with total number of clusters
    //
    vector<pair<int, unsigned int>> clusters = c.multiClusters(data,type);

    if(type=="max"){

	    if(privez){
	            dout << ident << "\tSpecies\tMaxClusterSize\tTotalClusters\tTotalSpecies\tSpanningSpecies\n";
	    }

	    //for(auto ites=clusters.begin(); ites!=clusters.end(); ++ites){
	    //dout << ident << "\t" << ites->first << "\t" << ites->second << endl;}

	    dout << ident << "\t" << clusters[0].first << "\t" << clusters[0].second 
	    	 << "\t" << clusters[1].second << "\t" << clusters[2].second << "\t" << clusters[3].first << endl;
	        
    }
    else
    {
	    if(privez){
	            dout << ident << "\tSpecies\tClusterSize\n";
	    }
	    for(auto ites=clusters.begin(); ites!=clusters.end(); ++ites){
	    dout << ident << "\t" << ites->first << "\t" << ites->second << endl;}
	    dout << ident << "\t" << -1 << "\t" << -1 << endl;
    }

	// return PatchStats(data,NumSpecies,outFile,ident);
	return 0;
}

int IPSNeutral::MFStats(simplmat <double> &data, simplmat <double> &q,
	int minBox, int maxBox, int deltaBox,const char * outFile,const char * ident, const char option)
{
	return 	MultifractalSBA(data, q,const_cast<char *>(outFile) ,minBox, maxBox, deltaBox, option,const_cast<char *>(ident));
}


