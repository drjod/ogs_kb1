#ifndef PIPING_H
#define PIPING_H


#include "fluid.h"
#include "interface.h"
#include <stdexcept>

namespace contra
{

class Configuration;

class Piping
{
public:
	Piping() = default;
	Piping(PipingData pipingData, FluidData fluidData);
	void configure(int type);
	Configuration* get_configuration();
	Fluid get_fluid() { return fluid; }
	void set_flow(double L);
	void set_Q(double _Q) { Q = _Q; }

private:
	Fluid fluid;
	Configuration* configuration;

	double Q;	// flow rate
	double d_0_i;	// pipe 0 - inner radius
	double d_0_o;	// pipe 0 - outer radius
	double d_1_i;	// pipe 1 - inner radius
	double d_1_o;	// pipe 1 - outer radius
	double w;	// shank spacing - 0. for CX
	double lambda_0;
	double lambda_1;

	// secondary
	double d;	// d_1_i - d_0_o  - for CX
	double x;	// centre of filling
	double A_0; // area of pipe 0
	double A_1; // area of pipe 1

	double u_0;	// velocity - pipe 0
	double u_1;	// velocity - pipe 1
	double Re_0;	// Reynolds number - pipe 0
	double Re_1;	// Reynolds number - pipe 1
	double Nu_0;	// Nussel number - pipe 0
	double Nu_1;	// Nussel number - pipe 1

	//friend class Configuration;
	friend class Configuration_U;
	friend class Configuration_2U;
	friend class Configuration_CX;
};

}
#endif
