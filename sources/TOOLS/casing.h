#ifndef CASING_H
#define CASING_H

#include "interface.h"


namespace contra
{

class Casing
{
public:
	Casing() : D(0.), L(0.), N(0), lambda_g(0.) {}
	Casing(SegmentData segmentData) :
		D(segmentData.D), L(segmentData.L), N(segmentData.N),
		lambda_g(segmentData.lambda_g) {}

	double get_D() { return D; }
	double get_L() { return L; }
	double get_N() { return N; }
	double get_lambda_g() { return lambda_g; }

private:
	double D;	// borehole diameter
	double L;	// length
	long N;	// number of nodes (discretization)
	double lambda_g;
};

}

#endif
