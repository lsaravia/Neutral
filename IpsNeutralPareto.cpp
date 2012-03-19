//	Clase Base : CABase
//
//	Clase Asociada : Specie
//
//
//#pragma implementation

#include "IpsNeutral.h"
#include "bgi.hpp"
//#include "fortify.h"
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
    PrimeraEval = true;
    
	// Usa los parametros desde 1 en adelante
	// The first species have the rates for the neutral models, 
	// after this only colonizationRate is taken as the metacommunity proportion.
	Sp	= new SpecieNeutral[NumSpecies+1];

	//for(int i=0; i<NumSpecies; i++)
	//	  Sp[i].Init();

	}

void IPSNeutral::Evaluate(int modType)
	{
	int x,y,i;
    
    if( PrimeraEval )
    {
		cerr << "ERROR modelo no inicializado" << endl;
		exit(1);
    }
    switch(modType){
		case 1: 					// Not saturated Neutral 
			for(i=0; i<NumEvaluations; i++)
				{
				x = Rand(DimX-1);
				y = Rand(DimY-1);
				EvalCell(x,y);
				}
			break;
		case 2:						// Saturated Neutral
			for(i=0; i<NumEvaluations; i++)
				{
				x = Rand(DimX-1);
				y = Rand(DimY-1);
				EvalCellZero(x,y);
				}
			break;				
		case 3:						// Not Saturated Hierarchical 
			for(i=0; i<NumEvaluations; i++)
				{
				x = Rand(DimX-1);
				y = Rand(DimY-1);
				EvalCellHierarchy(x,y);    
				}
			break;				
		case 4:						// Saturated Hierarchical 
		default:
			for(i=0; i<NumEvaluations; i++)
				{
				x = Rand(DimX-1);
				y = Rand(DimY-1);
				EvalCellZeroHierarchy(x,y);    
				}
		}
	T++;

	}

// Not saturated neutral model
// 
void IPSNeutral::EvalCell(int x,int y)
{
    int actSp = C(x,y).Specie;
	double rnd = Rand();
    int x1,y1,i;


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
	}
	else
	{
        if(rnd<Sp[0].MortalityRate)
            C(x,y).Specie=0;
		else 
        {
        //
        // Dispersal using Euclidean distance, Norm 2
        //
        EuclideanDispersal(x,y,x1,y1);
        	
		if( C(x1,y1).Specie == 0 )
			C(x1,y1).Specie = actSp;
		}
	}
}

// Zero-sum saturated neutral model
// 
void IPSNeutral::EvalCellZero(int x,int y)
{
    int actSp = C(x,y).Specie;
	double rnd = Rand();
    int x1,y1,i;

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
        	EuclideanDispersal(x,y,x1,y1);
			int dSp= C(x1,y1).Specie;
        	if( dSp>0 )
				C(x,y).Specie = dSp;
        } 
	}
	else
	{
        if(rnd<Sp[0].MortalityRate)
		{
			// If a species dies is replaced by another from metacomunity or neighborhood 
			// 
			rnd = Rand();
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
				EuclideanDispersal(x,y,x1,y1);
				int dSp= C(x1,y1).Specie;
				if( dSp>0 )
					C(x,y).Specie = dSp;
			}
		}
	}
}

// Not saturated hierarchical competition model
// 
void IPSNeutral::EvalCellHierarchy(int x,int y)
{
    int actSp = C(x,y).Specie;
	double rnd = Rand();
    int x1,y1,i;

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
	}
	else
	{
        if(rnd<Sp[0].MortalityRate)
            C(x,y).Specie=0;
        else 
        {
			rnd = Rand();
			if(rnd<Sp[0].ColonizationRate)
			{
				rnd = Rand();
				for(i=1; i<=NumSpecies; i++)
					if(rnd<Sp[i].ColonizationRate)
					{							// La reemplaza si la especie invasora es menor
						if( actSp>i) 			// Falta la probabilidad ro de reemplazo
							C(x,y).Specie=i;
						break;
					}
			}
			else
			{
				//
				// Dispersal using Euclidean distance, Norm 2
				//
				//EuclideanDispersal(x,y,x1,y1);
				//int dSp= C(x1,y1).Specie;
				//if( dSp>0 && dSp < actSp)
				//	C(x,y).Specie = dSp;

				EuclideanDispersal(x,y,x1,y1);
				int & dSp= C(x1,y1).Specie;
				if( dSp > actSp  || (dSp==0) ) 
					dSp = actSp;
			}
        }
	}
}

// Zero-sum saturated hierarchical model
// 
void IPSNeutral::EvalCellZeroHierarchy(int x,int y)
{
    int actSp = C(x,y).Specie;
	double rnd = Rand();
    int x1,y1,i;

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
/*        else
        {
        //
        // Euclidean distance, Norm 2
        //
        	EuclideanDispersal(x,y,x1,y1);
			int dSp= C(x1,y1).Specie;
			if( dSp>0 )
				C(x,y).Specie = dSp;
        } */
	}
	else
	{
        if(rnd<Sp[0].MortalityRate)
		{
		// Si una especie muere es reemplazada por otra de la metacomunidad
		// o del entorno
			rnd = Rand();
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
				EuclideanDispersal(x,y,x1,y1);
				int dSp = C(x1,y1).Specie;
				if( dSp>0 )
					C(x,y).Specie = dSp;
			}
		}
		else
		{
			
			rnd = Rand();
			if(rnd<Sp[0].ColonizationRate)
			{
				rnd = Rand();
				for(i=1; i<=NumSpecies; i++)
					if(rnd<Sp[i].ColonizationRate)
					{							// La reemplaza si la especie invasora es menor
						if( actSp>i) 			// Falta la probabilidad ro de reemplazo
							C(x,y).Specie=i;
						break;
					}
			}
			else
			{
				//
				// Dispersal using Euclidean distance, Norm 2
				//
				EuclideanDispersal(x,y,x1,y1);
				int & dSp= C(x1,y1).Specie;
				if( dSp > actSp  || (dSp==0) ) 
					dSp = actSp;
			}		
		}
	}
}


// Not saturated Density dependent neutral model 
// 
void IPSNeutral::EvalCellDensity(int x,int y)
{
    int actSp = C(x,y).Specie;
	double rnd = Rand();
    int x1,y1,i;


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
	}
	else
	{
		double rho=ParetoNeighborhood(x, y);
        if(rnd<Sp[0].MortalityRate)
            C(x,y).Specie=0;
		else 
        {
        //
        // Dispersal using Euclidean distance, Norm 2
        //
        EuclideanDispersal(x,y,x1,y1);
        	
		if( C(x1,y1).Specie == 0 )
			C(x1,y1).Specie = actSp;
		}
	}
}

double IPSNeutral::ParetoNeighborhood( int x, int y)
{
	int xf,yf,yi,dx,dy,dd;	
	dd=Sp[0].DispersalDistance;
	int k= Sp[0].ReplacementRate; // Es el parametro k de la funcion de Pareto 
	double sumTot=0,sumSp=0,dis;
	xf = x+dd;
	yf = y+dd;
	yi = y - dd;
	for(dx=x-dd; dx<=xf; dx++)
		for(dy=yi; dy<=yf; dy++)
			{
			x1 = (dx + DimX) % DimX;
			y1 = (dy + DimY) % DimY;
			dis = sqrt( static_cast<double>( (y-dy) * (y-dy) + (x-dx) * (x-dx) ));
			par = pow(1/dis, k);
			sumTot+= par;
			if( C(x1,y1).Specie>0)
				sumSp+=par;
			}
	return (sumSp/sumTot);
}



void IPSNeutral::InitParms(const bool pomac)
{
	int i=0;
	
	if(!pomac)
	{
		for(i=1; i<=NumSpecies; i++ )
		{
			GlobalRate += Sp[i].ColonizationRate; // Metacomunity proportion must be 1
		}
		if( !isEqual(GlobalRate,1) )
		{
			cerr << "Sum of proportion of species in metacommunity should be 1 : " << GlobalRate << endl;
			for(i=1; i<=NumSpecies; i++ )
				Sp[i].ColonizationRate /= GlobalRate;
		}

		for(i=2; i<=NumSpecies; i++)     // Probabilidad acumulada 
		{
			Sp[i].ColonizationRate += Sp[i-1].ColonizationRate;
		}

	}	
	
	if( Sp[0].ColonizationRate > 1)
		cerr << "The probability of colonization from the metacommunity can not exceed 1.\n";

	if( Sp[0].MortalityRate > Sp[0].BirthRate )
		cerr << "Mortality > Birth the process dies\n";

	// Only two event birth or death make the total rate of change for each site
	// GlobalRate = Sp[0].BirthRate > Sp[0].MortalityRate ? Sp[0].BirthRate : Sp[0].MortalityRate;
	GlobalRate = Sp[0].BirthRate + Sp[0].MortalityRate;

	Sp[0].BirthRate/=GlobalRate;
	Sp[0].MortalityRate/=GlobalRate;
//        Sp[0].BirthRate += Sp[0].MortalityRate; // Probabilidad acumulada
	Sp[0].ColonizationRate/= GlobalRate; // Tasa de colonizacion desde la metacomunidad
	
	NumEvaluations = DimX*DimY*GlobalRate;
	PrimeraEval=false;
	
}
