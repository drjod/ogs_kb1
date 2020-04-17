#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include "piping.h"
#include "greeks.h"

namespace contra
{

class Configuration
{
public:
	Configuration(Piping* _piping) : piping(_piping), R_0_Delta(0.), R_1_Delta(0.), R_01_Delta(0.) {}
	virtual Resistances set_resistances(double D, double lambda_g) = 0;
	virtual void set_flow(double L) = 0;
	virtual Greeks set_greeks(Piping* piping) = 0;

	virtual void set_functions(double& f1, double& f2, double& f3, const double& gz, const Greeks& greeks) = 0;
	virtual double F4(const double &z, const double &a, const double &b, const Greeks& greeks) = 0;
        virtual double F5(const double &z, const double &a, const double &b, const Greeks& greeks) = 0;

	double get_R_0_Delta() { return R_0_Delta; }
	double get_R_1_Delta() { return R_1_Delta; }
	double get_R_12_Delta() { return R_01_Delta; }

protected:
	Piping* piping;
	double R_0_Delta;
	double R_1_Delta;
	double R_01_Delta;

private:

};

class Configuration_U : public Configuration
{
public:
	Configuration_U(Piping* _piping);
	//void set_resistances_pipe(){}
	Resistances set_resistances(double D, double lambda_g);

	void set_flow(double L);
	Greeks set_greeks(Piping* piping);
	void set_functions(double& f1, double& f2, double& f3, const double& gz, const Greeks& greeks);
	double F4(const double &z, const double &a, const double &b, const Greeks& greeks);
        double F5(const double &z, const double &a, const double &b, const Greeks& greeks);
private:
	double R_adv;	// advective resistance
	double R_con_a;
	double R_con_b;

	double R_gs;
	double R_fg;
	double R_gg;


};

class Configuration_2U : public Configuration
{
public:
	Configuration_2U(Piping* _piping);
	//void set_resistances_pipe(){}
	Resistances set_resistances(double D, double lambda_g);

	void set_flow(double L);
	Greeks set_greeks(Piping* piping);
	void set_functions(double& f1, double& f2, double& f3, const double& gz, const Greeks& greeks);
	double F4(const double &z, const double &a, const double &b, const Greeks& greeks);
        double F5(const double &z, const double &a, const double &b, const Greeks& greeks);
private:
	double R_gs;
	double R_con_a;
	double R_con_b;
	double R_fg;
	double R_gg_1;
	double R_gg_2;
	double R_adv;	// advective resistance
};

class Configuration_CX : public Configuration
{
public:
	Configuration_CX(Piping* _piping);

	Resistances set_resistances(double D, double lambda_g); 

	void set_flow(double L);
	Greeks set_greeks(Piping* piping);
	void set_functions(double& f1, double& f2, double& f3, const double& gz, const Greeks& greeks);
	double F4(const double &z, const double &a, const double &b, const Greeks& greeks);
        double F5(const double &z, const double &a, const double &b, const Greeks& greeks);
private:
	double R_gs;
	double R_con_0_a; 
	double R_con_1_a; 
	double R_con_b;
	double R_ff;
	double R_fg;

	double R_adv_0;	// advective resistance - pipe 0
	double R_adv_1;	// advective resistance - pipe 1
	double R_adv_2;	// advective resistance
};

}

#endif
