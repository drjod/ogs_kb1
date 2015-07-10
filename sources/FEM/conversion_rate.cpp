
#include "conversion_rate.h"
#include <math.h>
#include <cmath>

//#define SIMPLE_KINETICS //wenn definiert, dann einfache Kinetik, sonst Schaube

#ifndef max
  #define max(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
  #define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

#define COMP_MOL_MASS_WATER 0.018016
#define COMP_MOL_MASS_N2   0.028013
#define COMP_MOL_MASS_O2 0.032

conversion_rate::conversion_rate(double T_solid, 
	                       double T_gas,  
						   double p_gas,
						   double x_reactive, 
						   double rho_s_initial,
						   double phi_S,
						   double delta_t,
						   FiniteElement::SolidReactiveSystem system)
:R(8.314510),p_eq(1.0) // , n_col(3)
{
	x = Eigen::VectorXd(1);
	update_param( T_solid, T_gas, p_gas, x_reactive, rho_s_initial, phi_S, delta_t, system);
	conversion_rate::rho_s_0 = rho_s_initial;
	if (system == FiniteElement::CaOH2){ //Definition auch in void CSolidProperties::SetSolidReactiveSystemProperties()
		rho_low = 1656.0;
		rho_up = 2200.0;
		reaction_enthalpy = -1.12e+05; //in J/mol; negative for exothermic composition reaction
		reaction_entropy = -143.5; //in J/mol K
		M_carrier = COMP_MOL_MASS_N2;
		M_react = COMP_MOL_MASS_WATER;
	}
	else if (system == FiniteElement::Mn3O4){//Definition auch in void CSolidProperties::SetSolidReactiveSystemProperties()
		rho_low = 4500.0;
		rho_up = 4860.0;
		reaction_enthalpy = -1.376e+05; //in J/mol; negative for exothermic composition reaction
		reaction_entropy = -114.1; //in J/mol K
		M_carrier = COMP_MOL_MASS_N2;
		M_react = COMP_MOL_MASS_O2;
	}

	tol_l = 1.0e-4;
	tol_u = 1.0 - tol_l;
	tol_rho = 0.1;
	
}

conversion_rate::~conversion_rate(void)
{
}

void conversion_rate::update_param(double T_solid, 
	                            double T_gas,  
								double p_gas,
								double x_reactive, 
								double rho_s_initial,
								double phi_S,
								double delta_t,
								FiniteElement::SolidReactiveSystem system)
{
	conversion_rate::T_s     = T_solid;
	conversion_rate::T       = T_gas; 
	conversion_rate::p_gas   = p_gas; // should be in unit bar
	conversion_rate::x_react   = x_reactive; 
	conversion_rate::rho_s   = rho_s_initial; 
	x(0)                  = rho_s_initial;
	conversion_rate::phi_solid = phi_S;
	conversion_rate::dt      = delta_t;
	conversion_rate::reaction_system = system;
}

void conversion_rate::calculate_qR()
{
	// step 1, calculate X_D and X_H
	X_D = (rho_s - rho_up - tol_rho)/(rho_low - rho_up - 2.0*tol_rho) ;
	X_D = (X_D < 0.5) ? max(tol_l,X_D) : min(X_D,tol_u); //constrain to interval [tol_l;tol_u]

	X_H = 1.0 - X_D;

	qR = 0.0;

	//Convert mass fraction into mole fraction
	double mol_frac_react;
	mol_frac_react = M_carrier*conversion_rate::x_react/(M_carrier*conversion_rate::x_react + M_react*(1.0-conversion_rate::x_react));
	
	////// step 2, calculate equilibrium
	T_s = conversion_rate::T_s;
	p_r_g = max(mol_frac_react * conversion_rate::p_gas, 1.0e-3); //avoid illdefined log

	// using the p_eq to calculate the T_eq - Clausius-Clapeyron
	T_eq = (reaction_enthalpy/R) / ((reaction_entropy/R) + log(p_r_g)); // unit of p in bar
	//Alternative: Use T_s as T_eq and calculate p_eq - for Schaube kinetics
	p_eq = exp((reaction_enthalpy/R)/T_s - (reaction_entropy/R));
	
	if (reaction_system == FiniteElement::CaOH2){ //Definition auch in void CSolidProperties::SetSolidReactiveSystemProperties()
		dXdt = Ca_hydration();
	}
	else if (reaction_system == FiniteElement::Mn3O4){//Definition auch in void CSolidProperties::SetSolidReactiveSystemProperties()
		dXdt = Mn_redox();
	}
    

	qR = (rho_up - rho_low) * dXdt;

}

void conversion_rate::set_rho_s(double new_rho_s)
{
	rho_s = new_rho_s; 
}

double conversion_rate::get_qR()
{
	return qR; 
}

void conversion_rate::get_x(Eigen::VectorXd& output_x)
{
	output_x = x; 
}


void conversion_rate::eval(double t, Eigen::VectorXd &y, Eigen::VectorXd &dydx)
{
	assert( y.size() == dydx.size() );

	this->set_rho_s( y(0) ); 
	this->calculate_qR();
	dydx(0) = this->get_qR();  

}

double conversion_rate::Ca_hydration()
{
	double dXdt;
		// step 3, calculate dX/dt
#ifdef SIMPLE_KINETICS
	if ( T_s < T_eq ) // hydration - simple model
#else
	if ( p_r_g > p_eq ) // hydration - Schaube model
#endif
	{
		//X_H = max(tol_l,X_H); //lower tolerance to avoid oscillations at onset of hydration reaction. Set here so that no residual reaction rate occurs at end of hydration.
#ifdef SIMPLE_KINETICS // this is from P. Schmidt
		dXdt = -1.0*(1.0-X_H) * (T_s - T_eq) / T_eq * 0.2 * conversion_rate::x_react;
#else //this is from Schaube
		if (X_H == tol_u || rho_s == rho_up)
			dXdt = 0.0;
		else if ( (T_eq-T_s) >= 50.0)
			dXdt = 13945.0 * exp(-89486.0/R/T_s) * pow(p_r_g/p_eq - 1.0,0.83) * 3.0 * (X_D) * pow(-1.0*log(X_D),0.666);
		else
			dXdt = 1.0004e-34 * exp(5.3332e4/T_s) * pow(p_r_g, 6.0) * (X_D);
#endif
	}
	else // dehydration
	{
		//X_D = max(tol_l,X_D); //lower tolerance to avoid oscillations at onset of dehydration reaction. Set here so that no residual reaction rate occurs at end of dehydration.
#ifdef SIMPLE_KINETICS // this is from P. Schmidt
		dXdt = -1.0* (1.0-X_D) * (T_s - T_eq) / T_eq * 0.05;
#else
		if (X_D == tol_u || rho_s == rho_low)
			dXdt = 0.0;
		else if (X_D < 0.2)
			dXdt = -1.9425e12 * exp( -1.8788e5/R/T_s )*pow(1.0-p_r_g/p_eq,3.0)*(X_H);
		else
			dXdt = -8.9588e9 * exp( -1.6262e5/R/T_s )*pow(1.0-p_r_g/p_eq,3.0)*2.0*pow(X_H, 0.5);
#endif
	}
	return dXdt;
}

double conversion_rate::Mn_redox()
{
	double dXdt, Avrami;

	if ( T_s < T_eq ) // oxidation
	{
		if (X_H == tol_u || rho_s == rho_up)
			dXdt = 0.0;
		else{
			Avrami = 282.277 - 0.912 * T_s + 9.949e-4 * T_s * T_s - 3.620e-7 * T_s * T_s * T_s;
			dXdt = 55271.0 * exp(-95.493e3 / (R * T_s)) * Avrami * (X_D) * pow(-1.0 * log(X_D),(Avrami - 1.0)/Avrami) * pow(1.0-T_s/T_eq,0.86);
		}
		if (p_r_g <= 1.0e-3)
			dXdt = 0.0;
	}
	else // reduction
	{
		dXdt = 0.0;
	}
	return dXdt;
}