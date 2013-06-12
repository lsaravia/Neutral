#include <ctime>

typedef unsigned long long int Ullint;
typedef unsigned int Uint;

struct Ranf1 {
	// Fast randon generator based on numerical recipes, the period is 1.8 10^19 .
	Ullint v;
	Ranf1(Ullint j=0) : v(4101842887655102017LL) {
		if( j==0 ) j=time(NULL);
		v ^= j;
		v = int64();
	}

	inline Ullint int64() {
		v ^= v >> 21; v ^= v << 35; v ^= v >> 4;
		return v * 2685821657736338717LL;
	}

	inline double doub() { return 5.42101086242752217E-20 * int64(); }

	inline Uint int32() { return (Uint)int64(); }
};


struct Ranfib {
	// Knuth’s subtractive generator using only floating operations. 
	// and an integer fast generator with period 1.8 10^19 .
	double dtab[55], dd;
	int inext, inextp;
	Ullint v;

	Ranfib(Ullint j=0) : inext(0), inextp(31), v(4101842887655102017LL) 
	{
		if( j==0 ) j=time(NULL);
		// integer generator
		v ^= j;
		v = int64();
		// floating point Constructor.
		Ranf1 init(j);
		for (int k=0; k<55; k++) dtab[k] = init.doub();
	}

	double doub() {
		//Returns random double-precision floating value between 0. and 1.
		if (++inext == 55) inext = 0;
		if (++inextp == 55) inextp = 0;
		dd = dtab[inext] - dtab[inextp];
		if (dd < 0) dd += 1.0;
		return (dtab[inext] = dd);
	}

	inline Ullint int64() {
		v ^= v >> 21; v ^= v << 35; v ^= v >> 4;
		return v * 2685821657736338717LL;
	}

};
