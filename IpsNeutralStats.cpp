#include "IpsNeutral.h"
#include "mf.h"

using namespace std;

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

int IPSNeutral::PStats(simplmat <double> &data, const char * outFile, const char * ident)
{
	// return PatchStats(data,NumSpecies,outFile,ident);
	return 0;
}

int IPSNeutral::MFStats(simplmat <double> &data, simplmat <double> &q,
	int minBox, int maxBox, int deltaBox,const char * outFile,const char * ident)
{
	return 	MultifractalSBA(data, q,const_cast<char *>(outFile) ,minBox, maxBox, deltaBox, 'S',const_cast<char *>(ident));
}


int IPSNeutral::MIStats(simplmat <double> &data, const char * outFile,const char * ident)
{
	// return MoranIRook(data, outFile, ident);
	return 0;
}

