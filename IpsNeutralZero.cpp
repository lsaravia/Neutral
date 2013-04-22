//	Clase Base : CABase
//
//	Clase Asociada : Specie
//
//
//#pragma implementation

#include "IpsNeutral.h"
#include "bgi.hpp"
#include "fortify.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <string>
#include <time.h>

using namespace std;

IPSNeutral::~IPSNeutral()
	{
	if(Sp!=NULL)
		delete Sp;
	}

void IPSNeutral::Init( unsigned numSp, int dimX, int dimY, int rr )
	{
	CABase::Init(numSp, dimX, dimY,rr);

	C.resize(dimX,dimY);

	//	ActualSp = 0;
    GlobalRate= 0;
    NumEvaluations = 0;
    
	// Usa los parametros desde 1 en adelante
	// The first species have the rates for the neutral models, 
	// after this only colonizationRate is taken as the metacommunity proportion.
	Sp	= new SpecieNeutral[NumSpecies+1];

	//for(int i=0; i<NumSpecies; i++)
	//	  Sp[i].Init();

	}

void IPSNeutral::Evaluate()
	{
	int x,y,i;
    static bool privez=true;

    if( privez )
    {
        for(i=1; i<=NumSpecies; i++ )
            GlobalRate += Sp[i].ColonizationRate; // Metacomunity proportion must be 1
        if( GlobalRate != 1 )
        {
			cerr << "Sum of proportion of species in metacommunity must be 1.\n";
			for(i=1; i<=NumSpecies; i++ )
				Sp[i].ColonizationRate /= GlobalRate;
		}
		
		if( Sp[0].ColonizationRate > 1)
			cerr << "La probabilidad de colonizar desde la metacomunidad no puede ser mayor que 1.\n";
			
		// Only two event birth or death make the total rate of change for each site
		GlobalRate = Sp[0].BirthRate+Sp[0].MortalityRate;

        Sp[0].BirthRate/=GlobalRate;
        Sp[0].MortalityRate/=GlobalRate;
//        Sp[0].BirthRate += Sp[0].MortalityRate; // Probabilidad acumulada
        Sp[0].ColonizationRate*= Sp[0].BirthRate; // Tasa de colonizacion desde la metacomunidad
        
        for(i=2; i<=NumSpecies; i++)
        {
            Sp[i].ColonizationRate += Sp[i-1].ColonizationRate;
        }
    	NumEvaluations = DimX*DimY*GlobalRate;
//  	max = DimX*DimY;
        privez=false;
    }
	for(i=0; i<NumEvaluations; i++)
		{
		x = Rand(DimX-1);
		y = Rand(DimY-1);
		EvalCell(x,y);
		}
	T++;

	}


void IPSNeutral::EvalCell(int x,int y)
{
    int actSp = C(x,y).Specie;
	double rnd = Rand();
    int dx,dy,dd,dis,x1,y1,i;

	// Colonization from the exterior
	//

	if(actSp==0)
	{
		if(rnd<Sp[0].ColonizationRate) // Colonize from metacommunity
		{
			rnd = Rand();
			for(i=1; i<=NumSpecies; i++)
				if(rnd<Sp[i].ColonizationRate)
				{
					C(x,y).Specie=i;
					break;
				}
        }
        else
        {
        //
        // Euclidean distance, Norm 2
        //
            dd=Sp[0].DispersalDistance;
    		dx= Rand(dd*2) - dd;
			dis=sqrt( static_cast<double>( dd * dd - dx * dx ));
			dy=Rand(dis*2) - dis ;
           	x1 = (x+ dx + DimX) % DimX;
        	y1 = (y+ dy + DimY) % DimY;
            C(x,y).Specie = C(x1,y1).Specie;
        }
	}
	else
	{
        if(rnd<Sp[0].MortalityRate)
		{
		// Si una especie muere es reemplazada por otra de la metacomunidad
		// o del entorno
			if(rnd<Sp[0].ColonizationRate) // Colonize from metacommunity
			{
				rnd = Rand();
				for(i=1; i<=NumSpecies; i++)
					if(rnd<Sp[i].ColonizationRate)
					{
						C(x,y).Specie=i;
						break;
					}
			}
			else
			{
			//
			// Euclidean distance, Norm 2
			//
				dd=Sp[0].DispersalDistance;
				dx= Rand(dd*2) - dd;
				dis=sqrt( static_cast<double>( dd * dd - dx * dx ));
				dy=Rand(dis*2) - dis ;
				x1 = (x+ dx + DimX) % DimX;
				y1 = (y+ dy + DimY) % DimY;
				C(x,y).Specie = C(x1,y1).Specie;
			}
		}
	}
}





