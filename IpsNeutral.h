#ifndef __Neutral_H
#define __Neutral_H
//#include "fortify.h"
#include "cabase.hpp"
#include "smattpl.hpp"
#include <cmath>  /* for std::abs(double) */
#include "ran.h"


// Multiespecific Contact process
//

class CellNeutral
{
	public:

	int Specie;

	CellNeutral() { Specie=0; };
	int & Elem() { return Specie; };

	const CellNeutral & Elem( int nsp ) {
								Specie=nsp;
								return *this; };

	CellNeutral & operator=(const CellNeutral &src){ Specie = src.Specie;
								return *this;};

	CellNeutral(const CellNeutral &src){ Specie = src.Specie; };

};


class SpecieNeutral
{
	public:
	double BirthRate;
	double ColonizationRate;
	double MortalityRate;
	double ReplacementRate;
	double DispersalDistance;

	SpecieNeutral() { Init(); };
	virtual void Init(){
						BirthRate=0;
						MortalityRate=0;
						DispersalDistance=0;
						ReplacementRate=0;
						ColonizationRate=0;
						};

	virtual int Scan(char * buff);

};


class IPSNeutral : public CABase
{
	protected:
		SpecieNeutral * Sp;
		simplmat <CellNeutral> C;
//		int ActualSp;
		CellNeutral ActualCell ;
        double GlobalRate;
        int NumEvaluations;
        bool PrimeraEval; 			// primera evaluacion de una corrida
//        Ranfib ran;
        Ranf1 ran;
	public:

	IPSNeutral( unsigned numSp, int dimX, int dimY, int rr=0 ) : CABase(numSp,dimX,dimY,rr) {Init(numSp,dimX,dimY,rr);};
	IPSNeutral(){};
	~IPSNeutral();
	void Init(unsigned numSp, int dimX, int dimY, int rr=0 );
	int  ReadParms( long rndSeed, char * file );
	int  ReadLineParms( const char * file );
	string PrintLineParms();

	void ReadSetSeed( char * fname);
	void RandomSetSeed(int sp,unsigned age, int no, int minX);
	int Rand(int num);
	double Rand();

	void Evaluate(int modType);
	void EvalCell(int x,int y);
	void EvalCellZero(int x,int y);
	void EvalCellHierarchy(int x, int y);
	void EvalCellZeroHierarchy(int x, int y);
	void Reset();
	void ReRun();
	void InitParms(const bool pomac);
	void EuclideanDispersal(const int &x,const int &y,int &x1,int &y1);
	void PowerDispersal(const int &x,const int &y,int &x1,int &y1);
	void ExpDispersal(const int &x,const int &y,int &x1,int &y1);

	void PrintGraph();
   	void InitGraph(char * idrPal);
   	void InitGraph(){};
	int  PrintDensity(const char *fname=NULL,const char *iname=NULL);
	int  PrintPomac(const char *fname=NULL,const char *iname=NULL);

	int ReadSeed( char * fname, int mode);
	int ReadSeed( char * fname){return ReadSeed( fname,0 );};
	int ReadIdrisi( char * fname, int mode=0);

	int SaveSeed( const char * fname);

	int PStats(simplmat <double> &data, const char * outFile, const char * ident);

	int Convert(simplmat <double> &data);                                      // Convierte a simplmat
	int Convert(simplmat <double> &data, const int * species );
	int Convert(simplmat <double> &data, const int specie );
	int AddConst(simplmat <double> &data, const double aa );
	
	int MFStats(simplmat <double> &data, simplmat <double> &q,
						int minBox, int maxBox, int deltaBox,const char * outFile,const char * ident);
						
	int MIStats(simplmat <double> &data, const char * outFile, const char * ident);       // Moran's I Rook 
	
				
};

inline void IPSNeutral::EuclideanDispersal(const int &x,const int &y,int &x1,int &y1)
{
	int dx,dy,dd,dis;	
	dd=Sp[0].DispersalDistance;
	while(true)
	{
		dx= Rand(dd*2) - dd;
		dis=sqrt( static_cast<double>( dd * dd - dx * dx ));
		dy=Rand(dis*2) - dis ;
		if(dx!=0 || dy!=0)
			break;
	}
	x1 = (x+ dx + DimX) % DimX;
	y1 = (y+ dy + DimY) % DimY;
}

//   Generate a inverse power dispersal function with parameter alfa = DispersalDistance
//   and Xmin=1 using the inverse of the cumulative F = U ^ (-1/(alfa-1))
//   
inline void IPSNeutral::PowerDispersal(const int &x,const int &y,int &x1,int &y1)
{
	int dx,dy;	
	double dd=Sp[0].DispersalDistance,dis;
	while(true)
	{
		double ang;
		dis= pow(Rand(), -1/(dd-1));
		ang = Rand() * Pi2;
		dx = cos( ang ) * dis ;
		dy = sin( ang ) * dis ;
		if(dx!=0 || dy!=0)
			break;
	}
	x1 = (x+ dx + DimX) % DimX;
	y1 = (y+ dy + DimY) % DimY;
	if( x1<0 || y1<0 )
	{
		x1= (x1+DimX)% DimX;
		y1= (y1+DimY)% DimY;
	}
	//cout << dx <<"\t" << dy << endl;
}

//   Generate an exponential dispersal function with parameter beta = DispersalDistance
//   
inline void IPSNeutral::ExpDispersal(const int &x,const int &y,int &x1,int &y1)
{
	int dx,dy;	
	double dd=Sp[0].DispersalDistance,dis;
	while(true)
	{
		double ang;
		//do ang=Rand(); while(ang==0.0); 
		dis= -log(Rand())/dd;
		ang = Rand() * Pi2;
		dx = cos( ang ) * dis ;
		dy = sin( ang ) * dis ;
		if(dx!=0 || dy!=0)
			break;
	}
	x1 = (x+ dx + DimX) % DimX;
	y1 = (y+ dy + DimY) % DimY;
	if( x1<0 || y1<0 )
	{
		x1= (x1+DimX)% DimX;
		y1= (y1+DimY)% DimY;
	}
}


inline double IPSNeutral::Rand() { 
	//return ranf();
	return ran.doub(); 
	};

inline int IPSNeutral::Rand(int num) {
		return (ran.int64() % (num+1)); // between 0 and num inclusive 
		//return ignuin(0,num);
        };



bool isEqual(double x, double y);

inline bool isEqual(double x, double y)
{
   const double epsilon = 1e-5;
   return std::abs(x - y) <= epsilon * std::abs(x);
   // see Knuth section 4.2.2 pages 217-218
} 

#endif // __Neutral_H
