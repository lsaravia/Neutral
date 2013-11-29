//
// Funciones Auxiliares
//
//#pragma implementation

#include <stdio.h>
#include <fstream>
#include <iomanip>
#include <sstream>
#include "IpsNeutral.h"
#include <ctype.h>
#include "grxkeys.h"
#include <time.h>
#include "bgi.hpp"
#include "RWFile.h"

#define uchar unsigned char

using namespace std;

int SpecieNeutral::Scan(char * buff)
	{
	int sp=0;
	istringstream ins(buff);
	ins >> sp >> BirthRate
				>> MortalityRate
				>> DispersalDistance
				>> ColonizationRate
				>> ReplacementRate;
	return sp;
	}


void IPSNeutral::Reset() // Corre nuevamente el modelo con parametros nuevos
	{
	C.fill( ActualCell.Elem(0) );
	A = 0;
	N = 1;
	T = 0;
	PrimeraEval=true;
	}

void IPSNeutral::ReRun()  // Corre nuevamente el modelo sin modificar parametros
	{
	C.fill( ActualCell.Elem(0) );
	A = 0;
	N = 1;
	T = 0;
	PrimeraEval=false;
	}

int IPSNeutral::ReadParms( long rndSeed, char * file )
	{
	ifstream in;
	char buff[255];

	int sp=0;


	in.open(file);
	if( !in )
		{
		cerr << "Cannot open Parms file.\n";
		return 1;
		}

	in >> DimX >> DimY;

	in >> NumSpecies;

	Init( NumSpecies, DimX, DimY, rndSeed);

	in.getline(buff,255);

	// Read only the neutral parameters
	in.getline(buff,254);
	sp = Sp[0].Scan(buff);
	
	for(int i=1;i<=NumSpecies;i++)
		{
		in.getline(buff,254);
		sp = Sp[i].Scan(buff);
		if( sp>NumSpecies || sp!=i )
			{
			cerr << "Error reading parameter file.\n";
			return 1;
			}
		}

	return 0;
	}


void IPSNeutral::ReadSetSeed(char * file)
	{
	FILE *in;
	char buff[255];

	int sp=0,age=0,cant=0,minx=0;


	if ((in = fopen(file, "rt")) == NULL)
		{
		fprintf(stderr, "Cannot open Parms file.\n");
		exit(1);
		}

	fgets(buff,80,in);
	do
		{

		if( sscanf(buff,"%i %i %i %i",&sp,&age,&cant,&minx) != 4 || sp>=NumSpecies )
			break;
		RandomSetSeed(sp,age,cant,minx);
		}
	while( fgets(buff,80,in)!=NULL );
	fclose(in);

	}

//  Read Seed file: Seed file must have SP as type.
//	mode == 1 -> Reads without model 
//	
int IPSNeutral::ReadSeed(char * fname,int mode)
	{
	FILE *in;
	char buff[256];
	int dx,dy;
	int spe,tipo=4,maxSpe=0;

	if ((in = fopen(fname, "rt")) == NULL)
		{
		fprintf(stderr, "Cannot open Seed file.\n");
		return 1;
		}
	fgets(buff,255,in);
	sscanf(buff,"%d %d",&dx,&dy);
	if( mode == 0 )
		{
		if( dx!=DimX || dy!=DimY )
			{
			fprintf(stderr, "Wrong dimension for Seed file.\n");
			exit(1);
			}
		}
	else
		{
		C.resize(dx,dy);
		DimX=dx;
		DimY=dy;
		N = A;
		Sp = NULL;
		}

	while( fgets(buff,4,in)!= NULL )
		{
		if( strncmp(buff,"SP",2)==0 )
			tipo = 0;
		else if(strncmp(buff,"AG",2)==0 )
			tipo = 1;
		else if(strncmp(buff,"BI",2)==0 )
			tipo = 3;
		else
			continue;

		for(dy=0;dy<DimY; dy++)
			for(dx=0;dx<DimX; dx++)
				{
				int ret = fscanf(in,"%i",&spe);
				if( ret == 0 || ret ==EOF )
					{
					cerr << "Seed File invalid.\n";
					exit(1);
					}

				switch (tipo) {
					case 0:
						C(dx,dy).Elem() = spe;
						if( spe>maxSpe )
							maxSpe = spe;
						break;

					default:
						cerr << "Seed File invalid.\n";
						exit(1);
					}

				}
		if( mode==0 )
			{
			if( maxSpe>NumSpecies )
				{
				fprintf(stderr, "Wrong Seed file.\n");
				exit(1);
				}
			}
		else
			{
			NumSpecies=maxSpe;
			}

		}

	fclose(in);
	return 0;

	}


int IPSNeutral::SaveSeed(const char * fname)
	{
	int i,j;
	ofstream sav( fname );
	if(!sav)
		{
		cerr << "Cannot open Seed file.\n";
		return 1;
		}

	sav << DimX << "\t" << DimY << endl;
	sav << "BI" << endl;
	for(i=0; i<DimY; i++)
		{
		for(j=0;j<DimX;j++)
			{
			sav<< setw(3) << int(C(j,i).Elem());
			}
		sav << endl;
		}
	sav << endl;

	return 0;
	}

void IPSNeutral::InitGraph(char * idrPal)
{
	IGraph(DimX,DimY,idrPal);
	char buff[20];

	for( int i=0; i<=NumSpecies; i++)
	{
		sprintf(buff,"%d",i);
		GLabel(buff,i);
		if(i>10) break;
	}
	
//  	BLabel();
}

// Agregar graficar la especie más abundante
//
void IPSNeutral::PrintGraph()
	{
	//int sa, baseName[9];
	ostringstream name;
	int sp=0;
    static bool privez=true;
    
	for(int i=0; i<DimX; i++)
		{
		for(int j=0;j<DimY;j++)
			{
			sp = C(i,j).Elem();
			PPix(i,j,sp);
			}
		}
    if( privez )
    {
		//while(GrKeyPressed()) GrKeyRead();
    	//getchar();
        privez = false;
    }
    
	if( GrKeyPressed() )
		{
		int sa=toupper(GrKeyRead());
		if( sa=='S')
			{
			EndGraph();
			cerr << "Tiempo :" << T << endl;
			cerr << "Ingrese nombre BASE de archivo : ";
//			cin.width(6);
//			cin >> baseName;
			name << "ipsNeutral" << T << ".sed" << ends;
			SaveSeed( name.str().c_str() );
			InitGraph();
			}
		else if( sa=='Q')
			exit(1);
			
		}
 
	}

void IPSNeutral::RandomSetSeed(int sp,unsigned age, int no, int minX)
	{
	int rx,ry;
	int i;
	for( i=0; i<no; i++)
		{
		while( 1 )
			{
			if( minX>0 )
				rx=Rand(minX);
			else
				rx=Rand(DimX-1);
			ry=Rand(DimY-1);
			if( C(rx,ry).Elem() == 0 )
				{
				C(rx,ry).Elem() = sp+1;
//				C(rx,ry).Age = age;
				break;
				}
			}
		}
	}

// DENTRO DE ESTA FUNCION CALCULA BIOMASA Y GUARDA y CALCULA ESPECTRO
//
int  IPSNeutral::PrintDensity(IPSParms p,const char *fname,const char *iname)
{
	ofstream dout;
	static int privez=0;
	double tot=0,totBio=0,totCells=DimX*DimY,diversity=0,freq=0,bioVol=0;
	simplmat <double> den(NumSpecies);
	//simplmat <int> calcDiv(NumSpecies);
	//double * den = new double[NumSpecies];
	int * calcDiv = new int[NumSpecies];
	int a,i;
	unsigned richness=0;

	
	if( fname!=NULL )
	{
		ostringstream name;
		name << fname << "Density.txt" << ends;
		dout.open( name.str().c_str(), ios::in );
		if( !dout )
			privez=1;

		dout.close();
        dout.clear();
		dout.open( name.str().c_str(), ios::app );
		if( !dout )
		{
			cerr << "Cannot open density file.\n";
			return 0;
		}
	}
	else
	{
		cerr << "File name cannot be NULL\n";
		return 0;		
	}

	if( privez )
	{
		privez=0;
		dout << iname <<"\tTime";
		for( a=0; a<NumSpecies; a++)
			{
			//dout.width(6);
			dout <<  "\t" << (a+1);
			}
		dout << "\tTot.Dens\tTot.Num\tRichness\tH\tRich>0.001\tH>0.001\tBiomass\tAvg.Bio" << endl;
	}

	if( iname==NULL )
		dout << T ;
	else
		dout << iname << "\t" << T ;
	
	// Calculates species densities
	//
	CalcDensity(den);
	for(i=0; i<NumSpecies; i++)
		calcDiv[i]=0;
	

	for( i=0; i<NumSpecies; i++)
	{
		freq = den(i)/totCells;
		dout << "\t" << freq;
		tot+= den(i);
		totBio+= freq;
		if(freq>0.0) 
		{
			diversity += freq * log(freq);
			richness++;
		}
		if( freq > 0.001 ) 				// Mark > 0.001 
			calcDiv[i] = true;
	}
	// If not saturated correct Shannon Diversity
	//
	if( tot<totCells)
	{
		diversity=0.0;
		for( i=0; i<NumSpecies; i++)
		{
			freq = den(i)/tot;
			if(freq>0.0) 
				diversity += freq * log(freq);
		}
	}
	
	dout << "\t" << totBio <<  "\t" << tot << "\t" << richness << "\t" << -1*diversity;
	
	// Calculates H diversity and richness for species with proportion > 0.001
	//
	totCells = diversity = richness = 0;
	for( i=0; i<NumSpecies; i++)
		if( calcDiv[i] )
			totCells+=den(i);
			
	for( i=0; i<NumSpecies; i++)
		if( calcDiv[i] )
		{
			freq = den(i)/totCells;
			diversity += freq * log(freq);
			richness++;
		}
			
	dout << "\t" << richness << "\t" << -1*diversity;

	if(p.bioCalc=='S')
	{
		ostringstream nam2;
		nam2 << iname << "\t" << T << ends;

		// Print the biomass spectrum
		bioVol = PrintDenBio(den, p.bioMax,p.bioMin,fname,nam2.str().c_str());   
	}
	dout << "\t" << bioVol << "\t" << bioVol/tot << endl;

	dout.close();
	delete []calcDiv;
	return tot;
};

// toma IPSparms para calcular correctamente biomasa y mfDIM
//
int  IPSNeutral::PrintPomac(IPSParms p, const char *fname,const char *iname)
{
	ofstream dout;
	static int privez=0;
	double tot=0,totBio=0,totCells=DimX*DimY,diversity=0,freq=0,maxf=0,bioVol=0;
	simplmat <double> den(NumSpecies);
	//double * den = new double[NumSpecies];
	int * calcDiv = new int[NumSpecies];
	int maxi=0,i;
	unsigned richness=0;


	if( fname!=NULL )
	{
		ostringstream name;
		name << fname << "Density.txt" << ends;
		dout.open( name.str().c_str(), ios::in );
		if( !dout )
			privez=1;

		dout.close();
        dout.clear();
		dout.open( name.str().c_str(), ios::app );
		if( !dout )
		{
			cerr << "Cannot open density file.\n";
			return 0;
		}
	}
	else
	{
		cerr << "File name cannot be NULL\n";
		return 0;		
	}

	if( privez )
	{
		privez=0;
		dout << iname <<"\tTime";
		dout << "\tTot.Dens\tTot.Num\tRichness\tH\tRich>0.001\tH>0.001\tBiomass\tAvg.Bio" << endl;
		}

	for(i=0; i<NumSpecies; i++)
		calcDiv[i]=0;
	
	// Calculate species' densities
	//
	CalcDensity(den);

	if( iname==NULL )
		dout << T ;
	else
		dout << iname << "\t" << T ;

	// Calculate Shannon Richness and the most abundant species
	//
	for( i=0; i<NumSpecies; i++)
	{
		freq = den(i)/totCells;
		tot+= den(i);
		totBio+= freq;
		if(freq>0.0) 
		{
			diversity += freq * log(freq);
			richness++;
			if( freq>maxf )
			{
				maxf=freq;
				maxi = i;
			}	
		}
		if( freq > 0.001 ) 				// Mark > 0.001 
			calcDiv[i] = true;
	}
	// If not saturated correct Shannon Diversity
	//
	if( tot<totCells)
	{
		diversity=0.0;
		for( i=0; i<NumSpecies; i++)
		{
			freq = den(i)/tot;
			if(freq>0.0) 
				diversity += freq * log(freq);
		}
	}

	dout << "\t" << totBio <<  "\t" << tot << "\t" << richness << "\t" << -1*diversity;
	
	// Calculates H diversity and richness for species with proportion > 0.001
	//
	totCells = diversity = richness = 0.0;
	for( i=0; i<NumSpecies; i++)
		if( calcDiv[i] )
			totCells+=den(i);
			
	for( i=0; i<NumSpecies; i++)
		if( calcDiv[i] )
		{
			freq = den(i)/totCells;
			diversity += freq * log(freq);
			richness++;
		}
			
	dout << "\t" << richness << "\t" << -1*diversity;

	ostringstream nam2;
	nam2 << iname << "\t" << T << ends;

	// Read q file for Multifractal spectrum
	//
	simplmat <double> dat;
	simplmat <double> q;
	RWFile file;
	if(!file.ReadSeed("q.sed", q))
		exit(1);

	// Calculates Dq converting to bioVolume =aN^(-4/3) 
	//
	if(p.bioCalc=='S')
	{
		ostringstream nam1;
		nam1 << fname << "mfBio.txt" << ends;
		bioVol = ConvertToBio(dat, den, p.bioMax,p.bioMin);
		// Print the biomass spectrum
		PrintDenBio(den, p.bioMax,p.bioMin,fname,nam2.str().c_str());   

		MFStats(dat,q,p.minBox,p.maxBox,p.deltaBox,nam1.str().c_str(),nam2.str().c_str());
	}
	dout << "\t" << bioVol << "\t" << bioVol/tot << endl;

	dout.close();
	
	// Calculates Dq reordering species from the most abundant with 1 
	//
	ostringstream nam3;
	nam3 << fname << "mfOrd.txt" << ends;

	if(Reordering(dat,den))
		MFStats(dat,q,p.minBox,p.maxBox,p.deltaBox,nam3.str().c_str(),nam2.str().c_str());
	
	

/*
	// Calculates Dq with q.sed for the most abundant specie
	//
	nam1 << fname << "mfCon.txt" << ends;

	Convert(dat,maxi);

	if( !MFStats(dat,q,4,512,20,nam1.str().c_str(),nam2.str().c_str()) )
	{
		Convert(dat,maxi);
		AddConst(dat,1.0);
		MFStats(dat,q,4,512,20,nam1.str().c_str(),nam2.str().c_str());
	}

	// Calculates infoDim of all species minus the most abundant
	//
	nam3 << fname << "mfSin.txt" << ends;
	for(i=0; i<NumSpecies; i++)
		calcDiv[i]=true;
	calcDiv[maxi]=false;
	dat.fill(0.0);
	Convert(dat,calcDiv);
	if(!MFStats(dat,q,4,512,20,nam3.str().c_str(),nam2.str().c_str()))
	{
		Convert(dat,calcDiv);
		AddConst(dat,1.0);
		MFStats(dat,q,4,512,20,nam3.str().c_str(),nam2.str().c_str());
	}
*/
	delete []calcDiv;					
	return tot;
	};

double IPSNeutral::PrintDenBio(simplmat <double> &den, float bioMax,float bioMin, const char * fname, const char * ident)
{
	static bool privez=false;
	ofstream dout;

	if( fname!=NULL )
	{
		ostringstream name;
		name << fname << "DenBio.txt" << ends;
		dout.open( name.str().c_str(), ios::in );
		if( !dout )
			privez=true;

		dout.close();
        dout.clear();
		dout.open( name.str().c_str(), ios::app );
		if( !dout )
		{
			cerr << "Cannot open density file.\n";
			return 0;
		}
	}
	else
	{
		cerr << "File name cannot be NULL\n";
		return 0;		
	}

	if( privez )
	{
		privez=false;
		dout << ident <<"\tTime";
		for( int i=0; i<NumSpecies; i++)
			{
			//dout.width(6);
			dout <<  "\t" << (i+1);
			}
		//dout << "\tTot.Dens\tTot.Num\tRichness\tH\tRich>0.001\tH>0.001" << endl;
		dout << "\tTot.Bio" << endl;
	}
	dout << ident;  // Parameters and time

/*
	for(int i=1; i<=NumSpecies; i++ )
	{	
		dout << "\t" << Sp[i].BirthRate;
	}
*/
	double dar=-4.0/3.0;		 //  inverse of Damuth Power exponent
	double a=0.0,totBio=0.0;  

	// set minimun value bioMin to a density of 90% of the total
	double maxN = DimX*DimY*0.9;
	// Calculate the constant 
	a = maxN/(pow(bioMin,(1/dar)));
	for(int i=0; i<NumSpecies; i++ )
	{	
		double bio  = pow(den(i)/a,dar);
	    if(bio>bioMax) 
	    	bio=bioMax;
	    else if(bio<bioMin) 
	    	bio=bioMin;

	    dout << "\t" << bio;
		totBio += bio*den(i);
	}
	dout << "\t" << totBio << endl;
	return(totBio);
}

int IPSNeutral::ReadIdrisi( char * fname, int mode )
	{
	ostringstream dname,iname;
	dname << fname << ".rdc" << ends;
	iname << fname << ".rst" << ends;

	char buff[255],type[10], *ptr;
	ifstream in;
	int cols=0,rows=0,dx,dy;

	in.open(dname.str().c_str());
	if( !in )
		{
		fprintf(stderr, "Cannot open doc file.\n");
		return 1;
		}
	while( !in.eof() )
		{
		in.read(buff,255);
		//in.width(20);
		//in >> buff
		ptr =strstr(buff,"columns");
		if( ptr!=NULL )
			cols = atoi(ptr+13);

		ptr=strstr(buff,"rows");
		if( ptr!=NULL )
			rows = atoi(ptr+13);

		ptr=strstr(buff,"data type");
		if( ptr!=NULL )
			{
			strncpy(type,ptr+14,10);
			type[9]='\0';
			}
		}

	in.close();
	if( mode == 0 )
		{
		if( cols!=DimX || rows!=DimY)
			return 1;
		}
	else
		{
		C.resize(cols,rows);
		DimX=cols;
		DimY=rows;
		N = A;
		}

	int spe=0,maxSpe=0;
	uchar echa=0;

#ifdef STREAMIO
	in.open(iname.str(), ios::binary | ios::in );

	if( !in )
		{
		cerr << "Cannot open img file.\n";
		return 1;
		}

	if( strstr( type,"integer")!=NULL || strstr( type,"byte")!=NULL)
		{
		for(dy=0;dy<DimY; dy++)
			for(dx=0;dx<DimX; dx++)
			{
			in.read(&echa,1);
			spe = echa;
			if( spe>maxSpe )
					maxSpe = spe;
			C(dx,dy).Elem(A) = spe;
			}
		}
	else
		{
		cerr << "Cannot read this file type\n";
		return 1;
		}
#else
	FILE * inn;
	if ((inn = fopen(iname.str().c_str(), "rb")) == NULL)
		{
		fprintf(stderr, "Cannot img file.\n");
		exit(1);
		}
	if( strstr( type,"integer")!=NULL || strstr( type,"byte")!=NULL)
		{
		for(dy=0;dy<DimY; dy++)
			for(dx=0;dx<DimX; dx++)
			{
			fread(&echa,1,1,inn);
			spe = echa;
			if( spe>maxSpe )
					maxSpe = spe;
			C(dx,dy).Elem() = spe;
			}
		}
	else
		{
		cerr << "Cannot read this file type\n";
		return 1;
		}


#endif

	NumSpecies = maxSpe;
	return 0;
	}


int IPSNeutral::ReadLineParms( const char * file)
{
	static ifstream inLineFile;
	string buff;
	static bool privez=true;
	if(privez)
	{
		inLineFile.open(file);
		if( !inLineFile )
		{
			cerr << "Cannot open Parms file" << endl;
			return 0;
		}
		privez=false;
		getline(inLineFile,buff);
	}
	
	getline(inLineFile,buff);
    if( inLineFile.eof() )
        return 0;
	
	istringstream ins(buff.c_str());
	
	ins >> Sp[0].BirthRate;
	ins >> Sp[0].MortalityRate;
	ins >> Sp[0].DispersalDistance;
	ins >> Sp[0].ColonizationRate;
	ins >> Sp[0].ReplacementRate;
	
	if(Sp[0].BirthRate<0 || Sp[0].DispersalDistance<0 || Sp[0].MortalityRate<0 || Sp[0].ColonizationRate<0 || Sp[0].ReplacementRate<0)
	{
		cerr << "Invalid parameter reading " << file << " file!!!!\b\b\b\b\b\b\b" << endl;
		exit(1);
	}

	
	return 1;
}

string IPSNeutral::PrintLineParms()
{
	ostringstream fb;
	fb << Sp[0].BirthRate << "\t";
	fb << Sp[0].MortalityRate << "\t";
	fb << Sp[0].DispersalDistance << "\t";
	fb << Sp[0].ColonizationRate << "\t";
	fb << Sp[0].ReplacementRate;
	
 	return fb.str();
}

