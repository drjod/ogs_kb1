#include <vector>
//using namespace std; 
class Fluid
{
private:
public:
	Fluid(void);
	~Fluid(void);


/* Data */
	static double TT, PP;

	//The Properties of Gases and Liquids, Fifth Edition   Bruce E. Poling, John M. Prausnitz, John P. O’Connell
	typedef struct
	{
		std::string name;
		std::string formula;
		double Tc; //critical temperature, K
		double Vc; //critical volume, cm3/mol
		double Pc; //critical pressure, bar
		double w;  //omega, acentric factor (See Chap. 2)
		double u;  //mu, dipole moment, D(debye)
		double k;  //association factor Eq. (9-4.11) 
		double M;  //molecular weight, g/mol
		double x;  //mole fraction
	}compound_properties; // eq 9-5.24 to eq 9-5.44
	static std::vector<compound_properties> Component;



/* Methods */
	static double co2_viscosity (double rho, double T);
	static void initial(void);
	static double viscosity_Chung(double T, std::vector<std::string> component_formula, std::vector<double> mole_amount); //low pressure condition
	static double viscosity_Chung_p(double T, double P, double rho, std::vector<std::string> component_formula, std::vector<double> mole_amount); //high pressure condition
	static double viscosity_Chung_CO2(double T, double P);
	static double viscosity_Chung_H2O(double T, double P);
	static double viscosity_TR(double T, double P, std::vector<std::string> component_formula, std::vector<double> mole_amount);

	static void fluid_properties_calc(double T, double P, double mCO2, double &rho_g, double &rho_l, double &solu_CO2, double &eta_g, double &eta_l);
	static void fluid_properties_calc_CH4(double T, double P, double mCH4, double &rho_g, double &rho_l, double &solu_CH4, double &eta_g, double &eta_l);

	static void fluid_properties_calc_H2(double T, double P, double mH2,  double &rho_g, double &rho_l, double &solu_H2, double &eta_g, double &eta_l);

	static double Gex_RT(double T, double P, double m);

	static double ion_limiting_conductivity(double T, double P, std::string ion);

	static double viscosity_XCl(double T, double P, double m);//to test the method

	static double viscosity_CO2_XCl(double T, double P, double mCO2, double mLiCl, double mNaCl, double mKCl);	



	static double Jones_Dole_viscosity(double T, double P, double mCO2, double mLiCl, double mNaCl, double mKCl);

	static double molarity_to_molality_NaCl(double T, double P, double c);

	static double molefraction_to_molality_NaCl(double x);

	static double molality_to_molarity_CO2(double T, double P, double mNaCl, double mCO2);

	static int viscosity_interface();
	static void entrance(void);
};

class NIST_H2
{
private:
public:
	NIST_H2(void);
	~NIST_H2(void);
/* Data */
	static int f;//f==0 normal hydrogen, f==1 parahydrogen, f==2 orthohydrogen
	static double TT, PP, DD, dP; //density unit mol dm-3
	static double R, Rm;
	static double M, Tc, Pc, Dc;
/* Methods */
	static void ReferenceConstants(void);
	static double alpha_0(double tau, double delta);
	static double alpha_r(double tau, double delta);
	static double dX_alpha_0(double tau, double delta);
	static double dY_alpha_0(double tau, double delta);
	static double dX_alpha_r(double tau, double delta);
	static double dY_alpha_r(double tau, double delta);
	static double pressure(double T, double D);
	static double dpressure(double D);
	static double density(double T, double P);
	static double volume(double T, double P);
	static double compressibility_factor(double T, double P);
	static double entropy(double T, double P);	
	static double interal_energy(double T, double P);
	static double gibbs_energy(double T, double P);
	static double helmholtz_energy(double T, double P);
	static double enthalpy(double T, double P);
	static double isochoric_heat_capacity(double T, double P);
	static double isobaric_heat_capacity(double T, double P);
	static double dQ(double dT);
	static double adiabatic_temperature(double T, double P, double P1, int n);
	static double adiabatic_temperature_id(double T, double P, double P1); //calculate as ideal gas
	static double lnphi(double T, double P);
	static double dY_alphar(double tau, double delta);
	static double uJT(double T, double P);
	static double dH(double dT);
	static double isenthalpic_temperature(double T, double P, double P1, int n);
};

//Roland Span, et al, 2000, JPCRD v29,1361
class NIST_N2
{
private:
public:
	NIST_N2(void);
	~NIST_N2(void);
/* Data */
	static double TT, PP, DD, dP; //density unit mol dm-3
	static double R, Rm;
	static double M, Tc, Pc, Dc;
/* Methods */
	static void ReferenceConstants(void);
	static double alpha_0(double tau, double delta);
	static double alpha_r(double tau, double delta);
	static double dX_alpha_0(double tau, double delta);
	static double dY_alpha_0(double tau, double delta);
	static double dX_alpha_r(double tau, double delta);
	static double dY_alpha_r(double tau, double delta);
	static double pressure(double T, double D);
	static double dpressure(double D);
	static double density(double T, double P);
	static double volume(double T, double P);
	static double compressibility_factor(double T, double P);
	static double entropy(double T, double P);	
	static double interal_energy(double T, double P);
	static double gibbs_energy(double T, double P);
	static double helmholtz_energy(double T, double P);
	static double enthalpy(double T, double P);
	static double isochoric_heat_capacity(double T, double P);
	static double isobaric_heat_capacity(double T, double P);
	static double dQ(double dT);
	static double adiabatic_temperature(double T, double P, double P1, int n);
	static double adiabatic_temperature_id(double T, double P, double P1); //calculate as ideal gas
	static double uJT(double T, double P);
	static double dH(double dT);
	static double isenthalpic_temperature(double T, double P, double P1, int n);

};

//U.Setzmann, W. Wagner, JPCRD, v20,1061, 1991
class NIST_CH4
{
private:
public:
	NIST_CH4(void);
	~NIST_CH4(void);
/* Data */
	static double TT, PP, DD, dP; //density unit mol dm-3
	static double R, Rm;
	static double M, Tc, Pc, Dc;
/* Methods */
	static void ReferenceConstants(void);
	static double alpha_0(double tau, double delta);
	static double alpha_r(double tau, double delta);
	static double dX_alpha_0(double tau, double delta);
	static double dY_alpha_0(double tau, double delta);
	static double dX_alpha_r(double tau, double delta);
	static double dY_alpha_r(double tau, double delta);
	static double pressure(double T, double D);
	static double dpressure(double D);
	static double density(double T, double P);
	static double volume(double T, double P);
	static double compressibility_factor(double T, double P);
	static double entropy(double T, double P);	
	static double interal_energy(double T, double P);
	static double gibbs_energy(double T, double P);
	static double helmholtz_energy(double T, double P);
	static double enthalpy(double T, double P);
	static double isochoric_heat_capacity(double T, double P);
	static double isobaric_heat_capacity(double T, double P);
	static double dQ(double dT);
	static double adiabatic_temperature(double T, double P, double P1, int n);
	static double adiabatic_temperature_id(double T, double P, double P1); //calculate as ideal gas
	static double uJT(double T, double P);
	static double dH(double dT);
	static double isenthalpic_temperature(double T, double P, double P1, int n);
};

//R. Span, W. Wagner, JPCRD, v25, 1509, 1996
class NIST_CO2
{
private:
public:
	NIST_CO2(void);
	~NIST_CO2(void);
/* Data */
	static double TT, PP, DD, dP; //density unit mol dm-3
	static double R, Rm;
	static double M, Tc, Pc, Dc;
/* Methods */
	static void ReferenceConstants(void);
	static double alpha_0(double tau, double delta);
	static double alpha_r(double tau, double delta);
	static double dX_alpha_0(double tau, double delta);
	static double dY_alpha_0(double tau, double delta);
	static double dX_alpha_r(double tau, double delta);
	static double dY_alpha_r(double tau, double delta);
	static double pressure(double T, double D);
	static double dpressure(double D);
	static double density(double T, double P);
	static double volume(double T, double P);
	static double compressibility_factor(double T, double P);
	static double entropy(double T, double P);	
	static double interal_energy(double T, double P);
	static double gibbs_energy(double T, double P);
	static double helmholtz_energy(double T, double P);
	static double enthalpy(double T, double P);
	static double isochoric_heat_capacity(double T, double P);
	static double isobaric_heat_capacity(double T, double P);
	static double dQ(double dT);
	static double adiabatic_temperature(double T, double P, double P1, int n);
	static double adiabatic_temperature_id(double T, double P, double P1); //calculate as ideal gas
	static double vapor_pressure(double T);
	static double saturated_liquid_density(double T);
	static double saturated_vapor_density(double T);
	static double uJT(double T, double P);
	static double dH(double dT);
	static double isenthalpic_temperature(double T, double P, double P1, int n);
};

//R. Schmidt, W. Wagner, 1985 Fluid Phase Equilibria, 19:175-200
class NIST_O2
{
private:
public:
	NIST_O2(void);
	~NIST_O2(void);
/* Data */
	static double TT, PP, DD, dP; //density unit mol dm-3
	static double R, Rm;
	static double M, Tc, Pc, Dc;
/* Methods */
	static void ReferenceConstants(void);
	static double alpha_0(double tau, double delta);
	static double alpha_r(double tau, double delta);
	static double dX_alpha_0(double tau, double delta);
	static double dY_alpha_0(double tau, double delta);
	static double dX_alpha_r(double tau, double delta);
	static double dY_alpha_r(double tau, double delta);
	static double pressure(double T, double D);
	static double dpressure(double D);
	static double density(double T, double P);
	static double volume(double T, double P);
	static double compressibility_factor(double T, double P);
	static double entropy(double T, double P);	
	static double interal_energy(double T, double P);
	static double gibbs_energy(double T, double P);
	static double helmholtz_energy(double T, double P);
	static double enthalpy(double T, double P);
	static double isochoric_heat_capacity(double T, double P);
	static double isobaric_heat_capacity(double T, double P);
	static double dQ(double dT);
	static double adiabatic_temperature(double T, double P, double P1, int n);
	static double adiabatic_temperature_id(double T, double P, double P1); //calculate as ideal gas
	static double uJT(double T, double P);
	static double dH(double dT);
	static double isenthalpic_temperature(double T, double P, double P1, int n);
};

class NIST_M
{
private:
public:
	NIST_M(void);
	~NIST_M(void);
/* Data */
	static double TT, PP, dP, VV; //density unit mol dm-3
	static double mH2, mN2, mO2, mCH4, mCO2;
/* Methods */
	//total enthalpy, unit kJ
	static double enthalpy(double T, double P, double mH2, double mN2, double mO2, double mCH4, double mCO2); 
	//total volume, unit dm^3
	static double volume(double T, double P, double mH2, double mN2, double mO2, double mCH4, double mCO2);
	//total interal energy, unit kJ
	static double interal_energy(double T, double P, double mH2, double mN2, double mO2, double mCH4, double mCO2);
	//total gibbs energy, unit kJ
	static double gibbs_energy(double T, double P, double mH2, double mN2, double mO2, double mCH4, double mCO2);
	//total helmholtz energy, unit kJ
	static double helmholtz_energy(double T, double P, double mH2, double mN2, double mO2, double mCH4, double mCO2);
	//unit g/dm^3
	static double density(double T, double P, double mH2, double mN2, double mO2, double mCH4, double mCO2);
	//unit J
	static double entropy(double T, /*double P,*/ double mH2, double mN2, double mO2, double mCH4, double mCO2);
	//total heat capacity, J/K
	static double isochoric_heat_capacity(double T, double P, double mH2, double mN2, double mO2, double mCH4, double mCO2);
	//total heat capacity, J/K
	static double isobaric_heat_capacity(double T, double P, double mH2, double mN2, double mO2, double mCH4, double mCO2);
	static double dQ(double dT);
	static double adiabatic_temperature(double T0, double P0, double P1, int n, double mH2, double mN2, double mO2, double mCH4, double mCO2);
	static double uJT(double T, double P, double mH2, double mN2, double mO2, double mCH4, double mCO2);
	static double dH(double dT);
	static double isenthalpic_temperature(double T0, double P0, double P1, int n, double mH2, double mN2, double mO2, double mCH4, double mCO2);
};
