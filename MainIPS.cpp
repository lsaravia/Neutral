//
// CADis : CA Discreto con bordes continuos
//
#include <ctime>
#include <ctype.h>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include "IpsNeutral.h"
#include "bgi.hpp"
#include "RWFile.h"

using namespace std;


struct IPSParms
{
	unsigned long nRuns;
	unsigned long nEvals;
	long rndSeed;
	char gr;
	int grBio;
	unsigned long inter;
	unsigned long init;
	bool pomac;
	int modType;
	char de;
	char sa;
	char baseName[40];
	char idrPal[40];
	char mfDim;
	float minBox;
	float maxBox;
	float deltaBox;
	char patchStat;
	char moranI;
};

int ReadParms(char * pFile, IPSParms &p);
timespec diff(timespec start, timespec end);

int main(int argc, char * argv[])
{
	IPSParms p;
	
	IPSNeutral ca;
	string pvals; // Para construir string con valores de parametros
	
	if( argc > 1 )
		{
		ReadParms( argv[1], p );
		ca.ReadParms( p.rndSeed, argv[2] );
		pvals = ca.PrintLineParms();
		ca.InitParms(false);
		if(argc >3 )
			{
			if( strstr(argv[3],"sed")!=NULL )
				ca.ReadSeed(argv[3]);	// Lee mapa ASCII de las especies y edades
			else if( strstr(argv[3],"set")!=NULL )
				ca.ReadSetSeed(argv[3]);
			}
		}
	else
		{
		cerr << "Parameter files missing! (xx.par yy.inp zz.sed)" << endl;
		exit(1);
		}

	int r=1; // Contador de nro de corridas
	
	while(1)
	{
		if( p.pomac )
		{
			if( !ca.ReadLineParms( "pomac.lin" ) ) // si no hay mas lineas de parametros no hace mas simulaciones
				break;
			pvals = ca.PrintLineParms();
			ca.InitParms(p.pomac);
		}
			
		if( p.gr=='S' )
        {
			ca.InitGraph(p.idrPal);
        }
        else
			cerr << "Run " << r << endl;

		timespec t0,t1;
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
        
		for(int i=0; i<p.nEvals; i++)
		{
            if( i == 0 )
			{
       	        if( p.gr=='S' )
                	ca.PrintGraph();
                else
   	                cerr << "Initial Conditions \n";
                if( p.de=='S' && !p.pomac)
                    ca.PrintDensity( p.baseName, pvals.c_str() );
                }
			ca.Evaluate(p.modType);
			
			if( ((i+1) % p.inter)==0 || i==0 )
			{
       	        if( p.gr=='S' )
                	ca.PrintGraph();
                else
					cerr << "Time " << (i+1) << endl;
				
				if( (i+1)>=p.init)
				{	
					if( p.pomac )
						ca.PrintPomac( p.baseName, pvals.c_str() );
					else
					{
						if( p.de=='S')
							if( ca.PrintDensity( p.baseName, pvals.c_str() )== 0 )
								break;
						if( p.sa=='S' )
							{
							ostringstream name;
							name << p.baseName << "-" << (i+1) << ".sed" << ends;
							ca.SaveSeed( name.str().c_str() );
							}
						if( p.patchStat=='S')
							{
							ostringstream name,nam1;
							name << p.baseName << "Pat.txt" << ends;
							nam1 << argv[2] << "-" << (i+1) << ends;
							simplmat <double> dat;
							ca.Convert(dat);
							ca.PStats(dat,name.str().c_str(),nam1.str().c_str());
							}
						if( p.moranI=='S')
							{
							ostringstream name,nam1;
							name << p.baseName << "MI.txt" << ends;
							nam1 << argv[2] << "-"<< (i+1) << ends;
							simplmat <double> dat;
							ca.Convert(dat);
							ca.MIStats(dat,name.str().c_str(),nam1.str().c_str());
							}
						if( p.mfDim=='S')
							{
							RWFile file;

							ostringstream name,nam1;
							name << p.baseName << "mfOrd.txt" << ends;
							nam1 << argv[2] << "-" << (i+1) << ends;
							simplmat <double> dat;
							simplmat <double> q;
							if(!file.ReadSeed("q.sed", q))
								exit(1);

							ca.Convert(dat);
							if(ca.Reordering(dat))
								ca.MFStats(dat,q,p.minBox,p.maxBox,p.deltaBox,name.str().c_str(),nam1.str().c_str());
							}
					}
				}


			}
		}
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
		cout << "Run " << r << "     : " << diff(t0,t1).tv_sec*1000+diff(t0,t1).tv_nsec/1000000<< endl;
		// (t1.tv_sec * 1000 + t1.tv_nsec/1000.0) - (t0.tv_sec * 1000 + t0.tv_nsec/1000.0) << endl;

		r++;
		if( r>p.nRuns && !p.pomac )
			break;
				
		if( p.pomac )
			ca.Reset();
		else
			ca.ReRun();
			
		if(argc >3 )
		{
			if( strstr(argv[3],"sed")!=NULL )
				ca.ReadSeed(argv[3]);	// Lee mapa ASCII de las especies y edades
			else if( strstr(argv[3],"set")!=NULL )
				ca.ReadSetSeed(argv[3]);
		}
	}
	
	if( p.gr=='S' )
		EGraph();

	return 0;
}


int ReadParms(char * pFile, IPSParms &p)
{
	string buff;
	ifstream parms(pFile);
    if( !parms )
    	{
    	cerr << "Can't open parms file" << endl;
        exit(1);
        }
	parms >> buff;
	while( !parms.eof() )
	{
		if(buff=="nRuns") // numero de Corridas
    	{
    		parms >> p.nRuns;
    	}
    	else if(buff=="nEvals") // numero de evaluaciones (max time)
    	{
    		parms >> p.nEvals;
    	}
    	else if(buff=="rndSeed") // random number seed
    	{
    		parms >> p.rndSeed;
    	}
    	else if(buff=="gr") // Graficos? N/S
    	{
    		parms >> p.gr;
    		p.gr = toupper( p.gr );
    	}
    	else if(buff=="grBio") // Grafica: 0 = Species, 1=Biomasa
    	{
    		parms >> p.grBio;
    	}
	    else if(buff=="inter") // Intervalo de Graficacion/Muestreo
	    {
    		parms >> p.inter;
	    }
	    else if(buff=="init") // tiempo de Inicio de muestreo
	    {
    		parms >> p.init;
	    }
   	    else if(buff=="modType") // 1 Modelo sin saturación 2 Modelo saturado
	    {
    		parms >> p.modType;
	    }

   	    else if(buff=="pomac") // Lee parametros tipo pomac
	    {
    		parms >> p.pomac;
	    }

	    else if(buff=="de") // Calcula densidad por especie S/N
	    {
    		parms >> p.de;
			p.de = toupper( p.de );
	    }
	    else if(buff=="sa") // Salva seed en cada muestreo
	    {
    		parms >> p.sa;
			p.sa = toupper( p.sa );
	    }
	    else if(buff=="baseName") // Salva seed en cada muestreo
	    {
    		parms >> p.baseName;
	    }
	    else if(buff=="idrPal") // utiliza una palette de Idrisi
	    {
    		parms >> p.idrPal;
	    }
   	    else if(buff=="mfDim") // Calculates generalized dimension 	(S/N)
	    {
    		parms >> p.mfDim;
			p.mfDim = toupper( p.mfDim );
		}
   	    else if(buff=="minBox") // q inicial para genDim
	    {
    		parms >> p.minBox;
		}
   	    else if(buff=="maxBox") // q final para genDim
	    {
    		parms >> p.maxBox;
		}
   	    else if(buff=="deltaBox") // Delta q para genDim
	    {
    		parms >> p.deltaBox;
		}
   	    else if(buff=="patchStat") // Patch statistics
	    {
    		parms >> p.patchStat;
			p.patchStat = toupper( p.patchStat );
		}
   	    else if(buff=="moranI") // Patch statistics
	    {
    		parms >> p.moranI;
			p.patchStat = toupper( p.moranI );
		}
		else if( !buff.empty() )
		{
			cerr << "Error in input file, unrecognized token" << endl;
			exit(1);
		}
		parms >> buff;
    }
    return 1;
}

timespec diff(timespec start, timespec end)
{
	timespec temp;
	if ((end.tv_nsec-start.tv_nsec)<0) {
		temp.tv_sec = end.tv_sec-start.tv_sec-1;
		temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
	} else {
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	}
	return temp;
}
