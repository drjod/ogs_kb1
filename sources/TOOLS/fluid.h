#ifndef FLUID_H
#define FLUID_H
#define _2_3 0.666667

#include <math.h>
#include "utilities.h"
#include "interface.h"

namespace contra
{

class Fluid
{
public:
	Fluid() = default;
	Fluid(FluidData fluidData) :
		lambda(fluidData.lambda), mu(fluidData.mu), c(fluidData.c), rho(fluidData.rho)
	{
		c_vol = c * rho;
		Pr = mu * c / lambda;
		DEBUG("Pr: " << Pr);
	} 

	double get_c_vol() { return c_vol; }
	double get_Pr() { return Pr; }
	double get_lambda() { return lambda; }
	double Reynolds(double u, double d) { return u * d * rho / mu; }
	double Nusselt_pipe(double Re, double L, double d) 
	{
		if(Re < 2300)
			return 4.365;
		else if(Re < 10000)
		{
			double gamma = (Re - 2300) / 7700;
			return (1 - gamma) * 4.365 + gamma * (1 + pow(d/L, _2_3)) * 38.8 * Pr / 
							(1 + .7880142765 * (pow(Pr, _2_3) - 1));
		}
		else
		{
			double xi_8 = pow(1.8 * log10(Re) - 1.5, -2.) / 8;
			DEBUG("xi: " << xi_8*8);
			DEBUG("d: " << d);
			return (1 + pow(d/L, _2_3)) * Re * Pr * xi_8 / (1 + 12.7 * sqrt(xi_8) * (pow(Pr, _2_3)-1));
		}
	}
	double Nusselt_ring(double Re, double L, double d_0_o, double d_1_i) 
	{
		const double d_oi = d_0_o / d_1_i;
		const double d_L = (d_1_i - d_0_o) / L;

		if(Re < 2300)
			return 3.66 + (4 - .102 /(d_oi + 0.02)) * pow(d_oi, .04) ;
		else if(Re < 10000)
		{
			double gamma = (Re - 2300) / 7700;
			return (1 - gamma) * ( 3.66 + (4 - .102 /(d_oi + 0.02)) * pow(d_oi, .04)  ) +
				gamma * ( (.86*pow(d_oi, .84) + (1 - 0.14*pow(d_oi, .6)))/ (1+d_oi)) * 
					(1 + pow(d_L, _2_3)) * 
					38.8 * Pr / (1 + .7880142765 * (pow(Pr, _2_3) - 1));
		}
		else
		{
			double xi_8 = pow(1.8 * log10(Re) - 1.5, -2.) / 8;
			DEBUG("xi: " << xi_8*8);
			return (1 + pow(d_L, _2_3)) * Re * Pr * xi_8 / (1 + 12.7 * sqrt(xi_8) * (pow(Pr, _2_3)-1))
				*  (.86*pow(d_oi, .84) + (1 - 0.14*pow(d_oi, .6)))/ (1+d_oi); 
		}
	}

private:
	double lambda;	// thermal conductivity
	double mu;	// viscosity
	double c;	// capacity
	double rho;	// density

	// secondary
	double c_vol;	// volumetric capacity
	double Pr;	// Prandl number
};

}
#endif
