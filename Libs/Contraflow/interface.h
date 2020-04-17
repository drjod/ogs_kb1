#ifndef INPUT_H_
#define INPUT_H_

#include "stru3_matrix.h"
#include <vector>

namespace contra
{

struct Resistances
{
	double R_0_Delta;
	double R_1_Delta;
};


struct SegmentData
{
	int N;
	double L;
	double D;
	double lambda_g;
};


struct FluidData
{
	double lambda;
	double mu;
	double c;
	double rho;
};


struct PipingData
{
	double d_0_i;
	double d_0_o;
	double d_1_i;
	double d_1_o;
	double w;  // shank space - 2U: across centre !
	double lambda_0;
	double lambda_1;
};

struct Result
{
	stru3::DVec T_in;
	stru3::DVec T_out;
	stru3::DVec T_fin;
	stru3::DVec T_fout;
	std::vector<Resistances> resistances_vec;
};

}

#endif
