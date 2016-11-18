#ifndef __CABASE_HPP
#define __CABASE_HPP

#include <math.h>

//void IGraph(int dimX,int dimY);
//void PPix(int x,int y, int color, int type=0);
//void EGraph();
//void GLabel(char * buff,int color);


class CABase
	{
	protected:
	unsigned A; 	// actual
	unsigned N; 	// proximo

	int DimX,DimY;
	unsigned NumSpecies;
    unsigned T;
	int rndSeed;
	static double Pi2;
    unsigned char Initialized;


	// random number between 0 and 1 inclusive 
	virtual double Rand()=0;     // { return ranf(); };
	// random number between 0 and num inclusive 
	virtual int Rand(int num)=0;  //	 { return ignuin(0,num); };

	public:

	CABase( unsigned nroSp, int dimX, int dimY, int rr=0) { Init(nroSp,dimX,dimY,rr); };
	CABase();
	void Init(unsigned nroSp, int dimX, int dimY, int rr=0 );
	virtual ~CABase();

	virtual void EvalCell(int x,int y)=0;

	virtual void Evaluate();
	virtual void EmptyPrevious(){};

	virtual void PrintGraph();
	virtual void InitGraph();
	virtual void EndGraph();

	int  ReadParms( int &nEval, int &nRuns, char * file ){return 1;};
	virtual int  ReadSeed( char * file ){return 1;};
	virtual int  SaveSeed( char * fname ){return 1;};

	};

#endif
