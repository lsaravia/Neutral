
typedef unsigned long long int Ullint;
typedef unsigned int Uint;


struct Ran {
	Ullint u,v,w;
	Ran(Ullint j) : v(4101842887655102017LL), w(1) {
	u = j ^ v; int64();
	v = u; int64();
	w = v; int64();
	}

	//Return 64-bit random integer. See text for explanation of method.
	inline Ullint int64() {
	u = u * 2862933555777941757LL + 7046029254386353087LL;
	v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
	w = 4294957665U*(w & 0xffffffff) + (w >> 32);
	Ullint x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
	return (x + v) ^ w;
	}
	// Return random double-precision floating value in the range 0. to 1.
	inline double doub() { return 5.42101086242752217E-20 * int64(); }
	//Return 32-bit random integer.
	inline Uint int32() { return (Uint)int64(); }
};

struct Ranq1 {
	//Recommended generator for everyday use. The period is
	//1:8 10^19 .
	Ullint v;
	Ranq1(Ullint j) : v(4101842887655102017LL) {
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
	// Implements Knuth’s subtractive generator using only floating operations. See text for cautions.
	double dtab[55], dd;
	int inext, inextp;
	Ranfib(Ullint j) : inext(0), inextp(31) {
	//Constructor. Call with any integer seed. Uses Ranq1 to initialize.
	Ranq1 init(j);
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

	inline unsigned long int32()
	//Returns random 32-bit integer. Recommended only for testing purposes.
	{ return (unsigned long)(doub() * 4294967295.0);}
};
