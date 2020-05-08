#ifndef SEGMENT_H
#define SEGMENT_H

#include <stdexcept>
#include <vector>
#include "casing.h"
#include "greeks.h"
#include "piping.h"
#include "configuration.h"
#include "stru3_matrix.h"
#include "interface.h"

namespace contra
{

class Segment{
public:
	Segment() {}
	Segment(SegmentData segmentData) : casing(segmentData),
		f1(stru3::DVec(casing.get_N())),
		f2(stru3::DVec(casing.get_N())),
		f3(stru3::DVec(casing.get_N()))
	{}

	Casing get_casing() { return casing; }
	Resistances set_resistances(Configuration* configuration);
	void set_functions(Piping* piping);

	void calculate_temperatures(Configuration* configuration);
	void configure(double* _T_in, double* _T_out, double* _T_s)
	{
		T_in = _T_in;
		T_out = _T_out;
		T_s = _T_s;
	}

	double get_f1(int n) { return f1[n]; }
	double get_f2(int n) { return f2[n]; }
	double get_f3(int n) { return f3[n]; }
	Greeks get_greeks() { return greeks; }
private:
	Casing casing;
	Greeks greeks;
	double* T_in;
	double* T_out;
	double* T_s;
	
	stru3::DVec f1, f2, f3;	// array

};

}
#endif
