#include <iostream>
#include <math.h>
#include <vector>
#include "eqbrm.h"
#include "Fluid.h"
#include "VLE.h"
#include "Density.h"
#include "PITZcalc.h"
#include "IAPWS-IF97.h"
#include "IO.h"
#include "DiffOH.h"
#include "NR.h"
//using namespace std;


//J.W.Leachmann, J. Phys. Chem. Ref. Data, Vol.38, No.3, Page, 721, 2009,
int NIST_H2::f;
double NIST_H2::TT, NIST_H2::PP, NIST_H2::DD, NIST_H2::dP;
double NIST_H2::R, NIST_H2::Rm;
double NIST_H2::M, NIST_H2::Tc, NIST_H2::Pc, NIST_H2::Dc;

//-----------------------------
NIST_H2::NIST_H2(void){}
NIST_H2::~NIST_H2(void){}

void NIST_H2::ReferenceConstants(void){
	Rm = 8.31451   ;//kJ kmol-1 K-1
	M  = 2.02      ;//kg kmol-1 (molecular weight)
	Tc = 33.145    ;//K         (Table 1)
	Pc = 1.2964    ;//MPa       (Table 1)
	Dc = 15.508    ;//mol dm-3  (Table 1)
}

double NIST_H2::alpha_0(double tau, double delta){
		
	double res;

	if(NIST_H2::f==1){
		const double a[10] = {0,-1.4485891134,1.884521239,4.30256,13.0289,-47.7365,50.0013,-18.6261,0.993973,0.536078};
		const double b[10] = {0,0,0,-15.1496751472,-25.0925982148,-29.4735563787,-35.4059141417,-40.724998482,-163.7925799988,-309.2173173842};

		res=log(delta)+1.5*log(tau)+a[1]+a[2]*tau; //Eq.31
		for(int i=3;i<8;i++)
			res+=a[i]*log(1-exp(b[i]*tau));//Eq.31
	}
	else if(NIST_H2::f==2){
		const double a[10] = {0,-1.4675442336,1.884506886,2.54151,-2.3661,1.00365,1.22447,0,0,0};
		const double b[10] = {0,0,0,-25.7676098736,-43.4677904877,-66.0445514750,-209.7531607465,0,0,0};

		res=log(delta)+1.5*log(tau)+a[1]+a[2]*tau; //Eq.31
		for(int i=3;i<8;i++)
		res+=a[i]*log(1-exp(b[i]*tau));//Eq.31
	}
	else
	{
		const double a[10] = {0,-1.4579856475,1.888076782,1.616,-0.4117,-0.792,0.758,1.217,0,0}; //J.W.Leachmann,2009(Table 4)
		const double b[10] = {0,0,0,-16.0205159149,-22.6580178006,-60.0090511389,-74.9434303817,-206.9392065168,0,0}; //J.W.Leachmann,2009(Table 4)

		res=log(delta)+1.5*log(tau)+a[1]+a[2]*tau; //Eq.31
		for(int i=3;i<8;i++)
			res+=a[i]*log(1-exp(b[i]*tau));//Eq.31
	}
	return res;
}

double NIST_H2::alpha_r(double tau, double delta){
	double res = 0.;
	//J.W.Leachmann,2009, Table 5 and Table 6
	if(NIST_H2::f==1){
		const double N[15] = {0,-7.33375,0.01,2.60375,4.66279,0.68239,-1.47078,0.135801,-1.05327,0.328239,-0.0577833,0.0449743,0.0703464,-0.0401766,0.11951};
		const double t[15] = {0,0.6855,1,1,0.489,0.774,1.133,1.386,1.619,1.162,3.96,5.276,0.99,6.791,3.19};
		const int d[15] = {0,1,4,1,1,2,2,3,1,3,2,1,3,1,1};
		const int p[10]	= {0,0,0,0,0,0,0,0,1,1};
		const double phi[15] = {0,0,0,0,0,0,0,0,0,0,-1.7437,-0.5516,-0.0634,-2.1341,-1.777};
		const double beta[15] = {0,0,0,0,0,0,0,0,0,0,-0.194,-0.2019,-0.0301,-0.2383,-0.3253};
		const double gamma[15] = {0,0,0,0,0,0,0,0,0,0,0.8048,1.5248,0.6648,0.6832,1.493};
		const double D[15] = {0,0,0,0,0,0,0,0,0,0,1.5487,0.1785,1.28,0.6319,1.7104};

		for(int i=1;i<8;i++)
			res+=N[i]*pow(delta,d[i])*pow(tau,t[i]);	//Eq.32
		for(int i=8;i<10;i++)
			res+=N[i]*pow(delta,d[i])*pow(tau,t[i])*exp(-pow(delta,p[i]));//Eq.32
		for(int i=10;i<11;i++)
			res+=N[i]*pow(delta,d[i])*pow(tau,t[i])*exp(phi[i]*(delta-D[i])*(delta-D[i])+beta[i]*(tau-gamma[i])*(tau-gamma[i]));//Eq.32
	}
	else if(NIST_H2::f==2){
		const double N[15] = {0,-6.83148,0.01,2.11505,4.38353,0.211292,-1.00939,0.142086,-0.87696,0.804927,-0.710775,0.0639688,0.0710858,-0.087654,0.647088};
		const double t[15] = {0,0.7333,1,1.1372,0.5136,0.5638,1.6248,1.829,2.404,2.105,4.1,7.658,1.259,7.589,3.946};
		const int d[15] = {0,1,4,1,1,2,2,3,1,3,2,1,3,1,1};
		const int p[10]	= {0,0,0,0,0,0,0,0,1,1};
		const double phi[15] = {0,0,0,0,0,0,0,0,0,0,-1.169,-0.894,-0.04,-2.072,-1.306};
		const double beta[15] = {0,0,0,0,0,0,0,0,0,0,-0.4555,-0.4046,-0.0869,-0.4415,-0.5743};
		const double gamma[15] = {0,0,0,0,0,0,0,0,0,0,1.5444,0.6627,0.763,0.6587,1.4327};
		const double D[15] = {0,0,0,0,0,0,0,0,0,0,0.6366,0.3876,0.9437,0.3976,0.9626};

		for(int i=1;i<8;i++)
			res+=N[i]*pow(delta,d[i])*pow(tau,t[i]);	//Eq.32
		for(int i=8;i<10;i++)
			res+=N[i]*pow(delta,d[i])*pow(tau,t[i])*exp(-pow(delta,p[i]));//Eq.32
		for(int i=10;i<11;i++)
			res+=N[i]*pow(delta,d[i])*pow(tau,t[i])*exp(phi[i]*(delta-D[i])*(delta-D[i])+beta[i]*(tau-gamma[i])*(tau-gamma[i]));//Eq.32
	}
	else
	{
		const double N[15] = {0,-6.93643,0.01,2.1101,4.52059,0.732564,-1.34086,0.130985,-0.777414,0.351944,-0.0211716,0.0226312,0.032187,-0.0231752,0.0557346};
		const double t[15] = {0,0.6844,1,0.989,0.489,0.803,1.1444,1.409,1.754,1.311,4.187,5.646,0.791,7.249,2.986};
		const int d[15] = {0, 1, 4, 1, 1, 2, 2, 3, 1, 3, 2, 1, 3, 1, 1};
		const int p[10]	= {0, 0, 0, 0, 0, 0, 0, 0, 1, 1};
		const double phi[15] = {0,0,0,0,0,0,0,0,0,0,-1.685,-0.489,-0.103,-2.506,-1.607};
		const double beta[15] = {0,0,0,0,0,0,0,0,0,0,-0.171,-0.2245,-0.1304,-0.2785,-0.3967};
		const double gamma[15] = {0,0,0,0,0,0,0,0,0,0,0.7164,1.3444,1.4517,0.7204,1.5445};
		const double D[15] = {0,0,0,0,0,0,0,0,0,0,1.506,0.156,1.736,0.67,1.662};
	
		for(int i=1;i<8;i++)
			res+=N[i]*pow(delta,d[i])*pow(tau,t[i]);	//Eq.32
		for(int i=8;i<10;i++)
			res+=N[i]*pow(delta,d[i])*pow(tau,t[i])*exp(-pow(delta,p[i]));//Eq.32
		for(int i=10;i<11;i++)
			res+=N[i]*pow(delta,d[i])*pow(tau,t[i])*exp(phi[i]*(delta-D[i])*(delta-D[i])+beta[i]*(tau-gamma[i])*(tau-gamma[i]));//Eq.32
	}
	return res;
}


double NIST_H2::dX_alpha_0(double tau, double delta){ return NR::dfridrX(alpha_0,tau,delta);}
double NIST_H2::dY_alpha_0(double tau, double delta){ return NR::dfridrY(alpha_0,tau,delta);}
double NIST_H2::dX_alpha_r(double tau, double delta){ return NR::dfridrX(alpha_r,tau,delta);}
double NIST_H2::dY_alpha_r(double tau, double delta){ return NR::dfridrY(alpha_r,tau,delta);}



double NIST_H2::dY_alphar(double tau, double delta){
	const double N[15] = {0,-6.93643,0.01,2.1101,4.52059,0.732564,-1.34086,0.130985,-0.777414,0.351944,-0.0211716,0.0226312,0.032187,-0.0231752,0.0557346};
	const double t[15] = {0,0.6844,1,0.989,0.489,0.803,1.1444,1.409,1.754,1.311,4.187,5.646,0.791,7.249,2.986};
	const int d[15] = {0, 1, 4, 1, 1, 2, 2, 3, 1, 3, 2, 1, 3, 1, 1};
	const int p[10]	= {0, 0, 0, 0, 0, 0, 0, 0, 1, 1};
	const double phi[15] = {0,0,0,0,0,0,0,0,0,0,-1.685,-0.489,-0.103,-2.506,-1.607};
	const double beta[15] = {0,0,0,0,0,0,0,0,0,0,-0.171,-0.2245,-0.1304,-0.2785,-0.3967};
	const double gamma[15] = {0,0,0,0,0,0,0,0,0,0,0.7164,1.3444,1.4517,0.7204,1.5445};
	const double D[15] = {0,0,0,0,0,0,0,0,0,0,1.506,0.156,1.736,0.67,1.662};
	double res;
	int i;
	res=0.0;
	for(i=1;i<8;i++)
		res+=d[i]*N[i]*pow(delta,d[i])*pow(tau,t[i]);	//Eq.32
	for(i=8;i<10;i++)
		res+=(d[i]-p[i]*pow(delta,p[i]))*N[i]*pow(delta,d[i])*pow(tau,t[i])*exp(-pow(delta,p[i]));//Eq.32
	for(i=10;i<11;i++)
		res+=(d[i]+2*delta*phi[i]*(delta-D[i]))*N[i]*pow(delta,d[i])*pow(tau,t[i])*exp(phi[i]*(delta-D[i])*(delta-D[i])+beta[i]*(tau-gamma[i])*(tau-gamma[i]));//Eq.32
	return res/delta;
}


double NIST_H2::pressure(double T, double D){
	NIST_H2::ReferenceConstants();
	double tau,delta;
	delta=D/Dc;
	tau=Tc/T;
	return 1.0e-3*D*Rm*T*(1+delta*NR::dfridrY(alpha_r,tau,delta)); 
	//Eq. 56  MPa unit of pressure J.Phys.Chem.Ref.Data, vol.29, No.6,2000
}

double NIST_H2::dpressure(double D){ 
	return NIST_H2::PP-NIST_H2::pressure(NIST_H2::TT,D);
}

double NIST_H2::density(double T, double P){
	double Ttemp=NIST_H2::TT, Ptemp=NIST_H2::PP,res;
	NIST_H2::TT=T;
	NIST_H2::PP=P;
	res=NR::zbrent(dpressure, 0.001, 100.0, 1.0e-8);
	NIST_H2::TT = Ttemp;
	NIST_H2::PP = Ptemp;
	return res;
}

double NIST_H2::volume(double T, double P){
	return 1.0/density(T,P);
}

double NIST_H2::compressibility_factor(double T, double P){
	double tau,delta;
	delta= NIST_H2::density(T,P)/Dc;
	tau=Tc/T;
	return 1.0+delta*NR::dfridrY(alpha_r,tau,delta); //Eq.57
}

double NIST_H2::entropy(double T, double P){
	double tau,delta;
	delta= NIST_H2::density(T,P)/Dc;
	tau=Tc/T;
	return Rm*(tau*(NR::dfridrX(alpha_0,tau,delta)+NR::dfridrX(alpha_r,tau,delta))-alpha_0(tau,delta)-alpha_r(tau,delta));//Eq.60
}

double NIST_H2::interal_energy(double T, double P){
	double tau,delta;
	delta= NIST_H2::density(T,P)/Dc;
	tau=Tc/T;
	return 1.0e-3*Rm*T*tau*(NR::dfridrX(alpha_0,tau,delta)+NR::dfridrX(alpha_r,tau,delta)); //Eq.58
}

double NIST_H2::gibbs_energy(double T, double P){
	double tau,delta;
	delta= NIST_H2::density(T,P)/Dc;
	tau=Tc/T;
	return 1.0e-3*Rm*T*(1.0+alpha_0(tau,delta)+alpha_r(tau,delta)+delta*(NR::dfridrY(alpha_r,tau,delta)));//Eq.61
}

double NIST_H2::helmholtz_energy(double T, double P){
	double tau,delta;
	delta= NIST_H2::density(T,P)/Dc;
	tau=Tc/T;
	return 1.0e-3*Rm*T*(alpha_0(tau,delta)+alpha_r(tau,delta));
}

double NIST_H2::enthalpy(double T, double P){
	double tau,delta,dr0,dr1,dr2;
	delta= NIST_H2::density(T,P)/Dc;
	tau=Tc/T;
	dr0=NR::dfridrX(alpha_0,tau,delta);
	dr1=NR::dfridrX(alpha_r,tau,delta);
	dr2=NR::dfridrY(alpha_r,tau,delta);
	return 1.0e-3*Rm*T*(tau*(dr0+dr1)+delta*dr2+1);//Eq.59 unit kJ/mol for enthalpy of H2
}

double NIST_H2::isochoric_heat_capacity(double T, double P){
	double tau,delta;
	delta= NIST_H2::density(T,P)/Dc;
	tau=Tc/T;
	return -Rm*tau*tau*(NR::dfridrX(dX_alpha_0,tau,delta)+NR::dfridrX(dX_alpha_r,tau,delta));//Eq.62
}

double NIST_H2::isobaric_heat_capacity(double T, double P){
	double tau,delta,Cv,term1,term2;
	delta= NIST_H2::density(T,P)/Dc;
	tau=Tc/T;
	Cv=NIST_H2::isochoric_heat_capacity(T,P);
	term1=1.0+delta*NR::dfridrY(alpha_r,tau,delta)-delta*tau*NR::dfridrY(dX_alpha_r,tau,delta);
	term1=term1*term1;
	term2=1.0+2.0*delta*NR::dfridrY(alpha_r,tau,delta)+delta*delta*NR::dfridrY(dY_alpha_r,tau,delta);
	return Cv+Rm*term1/term2;//Eq.63
}

double NIST_H2::dQ(double dT){
	double H1,H2,T1,P1,T2,P2,dH,res;
	T1=NIST_H2::TT;
	T2=NIST_H2::TT+dT;
	P1=NIST_H2::PP;
	P2=NIST_H2::PP+NIST_H2::dP;
	H1=NIST_H2::enthalpy(T1,P1) ;
	H2=NIST_H2::enthalpy(T2,P2);
	dH=H2-H1;
	res=dH-NIST_H2::dP/NIST_H2::DD;
	return res;
}

double NIST_H2::dH(double dT){
	double H1, H2, T1, P1, T2, P2, dH;
	T1 = NIST_H2::TT;
	T2 = NIST_H2::TT + dT;
	P1 = NIST_H2::PP;
	P2 = NIST_H2::PP + NIST_H2::dP;
	H1 = NIST_H2::enthalpy(T1, P1);
	H2 = NIST_H2::enthalpy(T2, P2);
	dH = H2 - H1;
	return dH;
}

double NIST_H2::adiabatic_temperature(double T, double P, double P1, int n){
	int i;
	double dT=0.,Tx,Px;
	NIST_H2::TT=T;
	NIST_H2::PP=P;
	NIST_H2::dP=(P1-P)/n;
	NIST_H2::DD=NIST_H2::density(T,P);
	Tx=NIST_H2::TT;
	Px=NIST_H2::PP;
	if(fabs(P1-P)>1.0e-6)
		for(i=0;i<n;i++){
			if(P1>P) dT=NR::zbrent(dQ, 1.0e-6, 20.0, 1.0e-8);
			if(P1<P) dT=NR::zbrent(dQ, -1.0e-6, -20.0, 1.0e-8);
			Tx+=dT;
			Px+=NIST_H2::dP;
			NIST_H2::TT=Tx;
			NIST_H2::PP=Px;
			NIST_H2::DD=NIST_H2::density(NIST_H2::TT,NIST_H2::PP);
		}
	return Tx;	
}

double NIST_H2::isenthalpic_temperature(double T, double P, double P1, int n){
	int i;
	double dT=0., Tx, Px;
	double uJT = NIST_H2::uJT(T, P);
	NIST_H2::TT = T;
	NIST_H2::PP = P;
	NIST_H2::dP = (P1 - P) / n;
	Tx = NIST_H2::TT;
	Px = NIST_H2::PP;
	if (fabs(P1 - P)>1.0e-6)
	for (i = 0; i<n; i++){
		if ((P1 - P)*uJT>0)
			dT = NR::zbrent(dH, 1.0e-6, 5.0, 1.0e-8);
		if ((P1 - P)*uJT<0)
			dT = NR::zbrent(dH, -1.0e-6, -5.0, 1.0e-8);
		Tx += dT;
		Px += NIST_H2::dP;
		NIST_H2::TT = Tx;
		NIST_H2::PP = Px;
	}
	return Tx;
}

double NIST_H2::adiabatic_temperature_id(double T, double P, double P1){
	double gamma,V,V1,Cv,Cp;
	Cp=28.824;//J K-1 mol-1
	Cv=Cp-Rm;//Chapter.2 Atkins' Physical Chemistry, ninth edition
	gamma=Cp/Cv;
	V=Rm*T/P;//ideal gas
	V1=V*pow(P/P1,1/gamma); //P*V^gamma = P1*V1^gamma Eq(2.29) =>  V1/V=(P/P1)^(1/gamma)
	return P1*V1/Rm;
}

double NIST_H2::lnphi(double T, double P){
	double tau,delta,res;
	delta= NIST_H2::density(T,P)/Dc;
	tau=Tc/T;
	res=NIST_H2::alpha_r(tau,delta)+delta*NIST_H2::dY_alphar(tau,delta)-log(delta*NIST_H2::dY_alphar(tau,delta)+1);
	return res;
}

double NIST_H2::uJT(double T, double P){
	double Cp,uT,res;
	Cp = NIST_H2::isobaric_heat_capacity(T, P);
	uT = NR::dfridrY(NIST_H2::enthalpy, T, P);
	res = -1.0e3*uT / Cp; //unit convert to K MPa-1
	return res;
}


//EoS for N2, R. Span et al, J. Phys. Chem. Ref. Data, vol.29, No.6, page 1361, 2000
double NIST_N2::TT, NIST_N2::PP, NIST_N2::DD, NIST_N2::dP;
double NIST_N2::R, NIST_N2::Rm;
double NIST_N2::M, NIST_N2::Tc, NIST_N2::Pc, NIST_N2::Dc;

void NIST_N2::ReferenceConstants(void){
	Rm = 8.31451   ;//kJ kmol-1 K-1 (molar gas constant)
	M  = 28.01348  ;//kg kmol-1 (molecular weight)
	Tc = 126.192   ;//K
	Pc = 3.3958    ;//MPa
	Dc = 11.1839   ;//mol dm-3
}

double NIST_N2::alpha_0(double tau, double delta){
	const double a[9] = {0,2.5,-12.76952708,-0.00784163,-1.9348193e-4,-1.247742e-5,6.678326e-8,1.012941,26.65788};
	return log(delta)+a[1]*log(tau)+a[2]+a[3]*tau+a[4]/tau+a[5]/tau/tau+a[6]/tau/tau/tau+a[7]*log(1-exp(-a[8]*tau));
}

double NIST_N2::alpha_r(double tau, double delta){
	//Table 17, 18
	const double N[37]
	={0,0.924803575275,-0.492448489428,0.661883336938,-0.192902649201e1,-0.622469309629e-1,0.349943957581,0.564857472498,
	   -0.161720005987e1,-0.481395031883,0.421150636384,-0.161962230825e-1,0.172100994165,0.735448924933e-2,0.168077305479e-1,
	   -0.107626664179e-2,-0.137318088513e-1,0.635466899859e-3,0.304432279419e-2,-0.435762336045e-1,-0.723174889316e-1,
		0.389644315272e-1,-0.212201363910e-1,0.408822981509e-2,-0.551990017984e-4,-0.462016716479e-1,-0.300311716011e-2,
		0.368825891208e-1,-0.255856846220e-2,0.896915264558e-2,-0.441513370350e-2,0.133722924858e-2,0.264832491957e-3,
		0.196688194015e2,-0.209115600730e2,0.167788306989e-1,0.262767566274e4};
	const double i[37]
	={0,1.0,1.0,2.0,2.0,3.0,3.0,1.0,1.0,1.0,3.0,3.0,4.0,6.0,6.0,7.0,7.0,8.0,8.0,1.0,2.0,3.0,4.0,5.0,8.0,4.0,5.0,5.0,8.0,3.0,5.0,6.0,9.0,1.0,1.0,3.0,2.0};
	const double j[37]
	={0,0.25,0.875,0.5,0.875,0.375,0.75,0.5,0.75,2,1.25,3.5,1,0.5,3,0,2.75,0.75,2.5,4,6,6,3,3,6,16,11,15,12,12,7,4,16,0,1,2,3};
	const double l[37]
	={0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,4,4,4,4,2,2,2,2};
	const double phi[37]
	={0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,20,20, 15,25};
	const double beta[37]
	={0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,325,325, 300,275};
	const double gamma[37]
	={0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,1.16,1.16, 1.13,1.25};
	double res=0;
	int k;
	for(k=1;k<=6;k++)
		res+=N[k]*pow(delta,i[k])*pow(tau,j[k]);
	for(k=7;k<=32;k++)
		res+=N[k]*pow(delta,i[k])*pow(tau,j[k])*exp(-pow(delta,l[k]));
	for(k=33;k<=36;k++)
		res+=N[k]*pow(delta,i[k])*pow(tau,j[k])*exp(-phi[k]*(delta-1)*(delta-1)-beta[k]*(tau-gamma[k])*(tau-gamma[k]));
	return res;
}

double NIST_N2::dX_alpha_0(double tau, double delta){ return NR::dfridrX(alpha_0,tau,delta);}
double NIST_N2::dY_alpha_0(double tau, double delta){ return NR::dfridrY(alpha_0,tau,delta);}
double NIST_N2::dX_alpha_r(double tau, double delta){ return NR::dfridrX(alpha_r,tau,delta);}
double NIST_N2::dY_alpha_r(double tau, double delta){ return NR::dfridrY(alpha_r,tau,delta);}


double NIST_N2::pressure(double T, double D){
	NIST_N2::ReferenceConstants();
	double tau,delta;
	delta=D/Dc;
	tau=Tc/T;
	return 1.0e-3*D*Rm*T*(1+delta*NR::dfridrY(alpha_r,tau,delta)); 
	//Eq. 56  MPa unit of pressure J.Phys.Chem.Ref.Data, vol.29, No.6,2000
}

double NIST_N2::dpressure(double D){ 
	return NIST_N2::PP-NIST_N2::pressure(NIST_N2::TT,D);
}

double NIST_N2::density(double T, double P){
	double Ttemp = NIST_N2::TT, Ptemp = NIST_N2::PP, res;
	NIST_N2::TT=T;
	NIST_N2::PP=P;
	res= NR::zbrent(dpressure, 0.001, 100.0, 1.0e-8);
	NIST_N2::TT = Ttemp;
	NIST_N2::PP = Ptemp;
	return res;
}

double NIST_N2::compressibility_factor(double T, double P){
	double tau,delta;
	delta= NIST_N2::density(T,P)/Dc;
	tau=Tc/T;
	return 1.0+delta*NR::dfridrY(alpha_r,tau,delta); //Eq.57
}

double NIST_N2::volume(double T, double P){
	return 1.0/density(T,P);
}

double NIST_N2::interal_energy(double T, double P){
	double tau,delta;
	delta= NIST_N2::density(T,P)/Dc;
	tau=Tc/T;
	return 1.0e-3*Rm*T*tau*(NR::dfridrX(alpha_0,tau,delta)+NR::dfridrX(alpha_r,tau,delta)); //Eq.58
}

double NIST_N2::enthalpy(double T, double P){
	double tau,delta,dr0,dr1,dr2;
	delta= NIST_N2::density(T,P)/Dc;
	tau=Tc/T;
	dr0=NR::dfridrX(alpha_0,tau,delta);
	dr1=NR::dfridrX(alpha_r,tau,delta);
	dr2=NR::dfridrY(alpha_r,tau,delta);
	return 1.0e-3*Rm*T*(tau*(dr0+dr1)+delta*dr2+1);//Eq.59 unit kJ/mol for enthalpy of H2
}

double NIST_N2::entropy(double T, double P){
	double tau,delta;
	delta= NIST_N2::density(T,P)/Dc;
	tau=Tc/T;
	return Rm*(tau*(NR::dfridrX(alpha_0,tau,delta)+NR::dfridrX(alpha_r,tau,delta))-alpha_0(tau,delta)-alpha_r(tau,delta));//Eq.60
}

double NIST_N2::gibbs_energy(double T, double P){
	double tau,delta;
	delta= NIST_N2::density(T,P)/Dc;
	tau=Tc/T;
	return 1.0e-3*Rm*T*(1.0+alpha_0(tau,delta)+alpha_r(tau,delta)+delta*(NR::dfridrY(alpha_r,tau,delta)));//Eq.61
}

double NIST_N2::helmholtz_energy(double T, double P){
	double tau,delta;
	delta= NIST_N2::density(T,P)/Dc;
	tau=Tc/T;
	return 1.0e-3*Rm*T*(alpha_0(tau,delta)+alpha_r(tau,delta));
}


double NIST_N2::isochoric_heat_capacity(double T, double P){
	double tau,delta;
	delta= NIST_N2::density(T,P)/Dc;
	tau=Tc/T;
	return -Rm*tau*tau*(NR::dfridrX(dX_alpha_0,tau,delta)+NR::dfridrX(dX_alpha_r,tau,delta));//Eq.62
}

double NIST_N2::isobaric_heat_capacity(double T, double P){
	double tau,delta,Cv,term1,term2;
	delta= NIST_N2::density(T,P)/Dc;
	tau=Tc/T;
	Cv=NIST_N2::isochoric_heat_capacity(T,P);
	term1=1.0+delta*NR::dfridrY(alpha_r,tau,delta)-delta*tau*NR::dfridrY(dX_alpha_r,tau,delta);
	term1=term1*term1;
	term2=1.0+2.0*delta*NR::dfridrY(alpha_r,tau,delta)+delta*delta*NR::dfridrY(dY_alpha_r,tau,delta);
	return Cv+Rm*term1/term2;//Eq.63
}

double NIST_N2::dQ(double dT){
	double H1,H2,T1,P1,T2,P2,dH,res;
	T1=NIST_N2::TT;
	T2=NIST_N2::TT+dT;
	P1=NIST_N2::PP;
	P2=NIST_N2::PP+NIST_N2::dP;
	H1=NIST_N2::enthalpy(T1,P1) ;
	H2=NIST_N2::enthalpy(T2,P2);
	dH=H2-H1;
	res=dH-NIST_N2::dP/NIST_N2::DD;
	return res;
}

double NIST_N2::adiabatic_temperature(double T, double P, double P1, int n){
	int i;
	double dT=0.,Tx,Px;
	NIST_N2::TT=T;
	NIST_N2::PP=P;
	NIST_N2::dP=(P1-P)/n;
	NIST_N2::DD=NIST_N2::density(T,P);
	Tx=NIST_N2::TT;
	Px=NIST_N2::PP;
	if(fabs(P1-P)>1.0e-6)
		for(i=0;i<n;i++){
			if(P1>P) dT=NR::zbrent(dQ, 1.0e-6, 20.0, 1.0e-8);
			if(P1<P) dT=NR::zbrent(dQ, -1.0e-6, -20.0, 1.0e-8);
			Tx+=dT;
			Px+=NIST_N2::dP;
			NIST_N2::TT=Tx;
			NIST_N2::PP=Px;
			NIST_N2::DD=NIST_N2::density(NIST_N2::TT,NIST_N2::PP);
		}
	return Tx;	
}

double NIST_N2::adiabatic_temperature_id(double T, double P, double P1){
	double gamma,V,V1,Cv,Cp;
	Cp=29.125;//J K-1 mol-1, Table 2.8
	Cv=Cp-Rm;//Chapter.2 Atkins' Physical Chemistry, ninth edition
	gamma=Cp/Cv;
	V=Rm*T/P;//ideal gas
	V1=V*pow(P/P1,1/gamma); //P*V^gamma = P1*V1^gamma Eq(2.29) =>  V1/V=(P/P1)^(1/gamma)
	return P1*V1/Rm;
}


double NIST_N2::isenthalpic_temperature(double T, double P, double P1, int n){
	int i;
	double dT=0., Tx, Px;
	double uJT = NIST_N2::uJT(T, P);
	NIST_N2::TT = T;
	NIST_N2::PP = P;
	NIST_N2::dP = (P1 - P) / n;
	Tx = NIST_N2::TT;
	Px = NIST_N2::PP;
	if (fabs(P1 - P)>1.0e-6)
	for (i = 0; i<n; i++){
		if ((P1 - P)*uJT>0)
			dT = NR::zbrent(dH, 1.0e-6, 5.0, 1.0e-8);
		if ((P1 - P)*uJT<0)
			dT = NR::zbrent(dH, -1.0e-6, -5.0, 1.0e-8);
		Tx += dT;
		Px += NIST_N2::dP;
		NIST_N2::TT = Tx;
		NIST_N2::PP = Px;
	}
	return Tx;
}


double NIST_N2::dH(double dT){
	double H1, H2, T1, P1, T2, P2, dH;
	T1 = NIST_N2::TT;
	T2 = NIST_N2::TT + dT;
	P1 = NIST_N2::PP;
	P2 = NIST_N2::PP + NIST_N2::dP;
	H1 = NIST_N2::enthalpy(T1, P1);
	H2 = NIST_N2::enthalpy(T2, P2);
	dH = H2 - H1;
	return dH;
}
double NIST_N2::uJT(double T, double P){
	double Cp, uT, res;
	Cp = NIST_N2::isobaric_heat_capacity(T, P);
	uT = NR::dfridrY(NIST_N2::enthalpy, T, P);
	res = -1.0e3*uT / Cp; //unit convert to K MPa-1
	return res;
}






//EoS for CH4, U. Setzmann, W. Wagner, J. Phys. Chem. Ref. Data, vol.20, No.6, page 1061, 1991
double NIST_CH4::TT, NIST_CH4::PP, NIST_CH4::DD, NIST_CH4::dP;
double NIST_CH4::R, NIST_CH4::Rm;
double NIST_CH4::M, NIST_CH4::Tc, NIST_CH4::Pc, NIST_CH4::Dc;

void NIST_CH4::ReferenceConstants(void){ //Table 35
	R  = 0.5182705 ;//kJ kg-1 K-1
	Rm = 8.31451   ;//kJ kmol-1 K-1 (molar gas constant)
	M  = 16.0428   ;//kg kmol-1 (molecular weight)
	Tc = 190.564   ;//K
	Pc = 4.5922    ;//MPa
	Dc = 162.66    ;//kg m-3
	Dc = Dc / M;    //mol dm-3
}

double NIST_CH4::alpha_0(double tau, double delta){
	const double a[9] = {0,9.91243972,-6.33270087,3.0016,0.008449,4.6942,3.4865,1.6572,1.4115};
	const double d[9] = {0,0,0,0,3.4004324,10.26951575,20.43932747,29.93744884,79.13351945};
	int i;
	double res;
	res=log(delta)+a[1]+a[2]*tau+a[3]*log(tau);
	for(i=4;i<=8;i++) 
		res+=a[i]*log(1.0-exp(-d[i]*tau));
	return res;
}

double NIST_CH4::alpha_r(double tau, double delta){
	//Table 17, 18
	const double N[41]={
		 0.0            ,
		 0.4367901028e-1,
		 0.6709236199   ,
		-0.1765577859e1 ,
		 0.8582330241   ,
		-0.1206513052e1 ,
		 0.5120467220   ,
		-0.4000010791e-3,
		-0.1247842423e-1,
		 0.3100269701e-1,
		 0.1754748522e-2,
		-0.3171921605e-5,
		-0.2240346840e-5,
		 0.2947056156e-6,
		 0.1830487909   ,
		 0.1511883679   ,
		-0.4289363877   ,
		 0.6894002446e-1,
		-0.1408313996e-1,
		-0.3063054830e-1,
		-0.2969906708e-1,
		-0.1932040831e-1,
		-0.1105739959   ,
		 0.9952548995e-1,
		 0.8548437825e-2,
		-0.6150555662e-1,
		-0.4291792423e-1,
		-0.1813207290e-1,
		 0.3445904760e-1,
		-0.2385919450e-2,
		-0.1159094939e-1,
		 0.6641693602e-1,
		-0.2371549590e-1,
		-0.3961624905e-1,
		-0.1387292044e-1,
		 0.3389489599e-1,
		-0.2927378753e-2,
		 0.9324799946e-4,
		-0.6287171518e1 ,
		 0.1271069467e2 ,
		-0.6423953466e1 };
	const double c[41]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,2,2,2,2,2,3,3,3,3,4,4,4,4,4,4,4};
	const double d[41]={0,1,1,1,2,2,2,2,3,4,4,8,9,10,1,1,1,2,4,5,6,1,2,3,4,4,3,5,5,8,2,3,4,4,4,5,6,2,0,0,0};
	const double t[41]={0.0,-0.5,0.5,1.0,0.5,1.0,1.5,4.5,0.0,1.0,3.0,1.0,3.0,3.0,0.0,1.0,2.0,0.0,0.0,2.0,2.0,5.0,
		5.0,5.0,2.0,4.0,12.0,8.0,10.0,10.0,10.0,14.0,12.0,18.0,22.0,18.0,14.0,2.0,0.0,1.0,2.0};
	const double phi[41]={0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,20,40,40,40};
	const double beta[41]={0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,200,250,250,250};
	const double gamma[41]={0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,1.07,1.11,1.11,1.11};	
	
	double res=0;
	int i;
	for(i=1;i<=13;i++)
		res+=N[i]*pow(delta,d[i])*pow(tau,t[i]);
	for(i=14;i<=36;i++)
		res+=N[i]*pow(delta,d[i])*pow(tau,t[i])*exp(-pow(delta,c[i]));
	for(i=37;i<=40;i++)
		res+=N[i]*pow(delta,d[i])*pow(tau,t[i])*exp(-phi[i]*(delta-1.0)*(delta-1.0)-beta[i]*(tau-gamma[i])*(tau-gamma[i]));
	return res;
}


double NIST_CH4::dX_alpha_0(double tau, double delta){ return NR::dfridrX(alpha_0,tau,delta);}
double NIST_CH4::dY_alpha_0(double tau, double delta){ return NR::dfridrY(alpha_0,tau,delta);}
double NIST_CH4::dX_alpha_r(double tau, double delta){ return NR::dfridrX(alpha_r,tau,delta);}
double NIST_CH4::dY_alpha_r(double tau, double delta){ return NR::dfridrY(alpha_r,tau,delta);}


double NIST_CH4::pressure(double T, double D){//density in mol/l
	NIST_CH4::ReferenceConstants();
	double tau,delta;
	delta=D/Dc;
	tau=Tc/T;
	return 1.0e-3*D*M*R*T*(1+delta*NR::dfridrY(alpha_r,tau,delta)); 
	//Eq. 56  MPa unit of pressure J.Phys.Chem.Ref.Data, vol.29, No.6,2000
}

double NIST_CH4::dpressure(double D){ 
	return NIST_CH4::PP-NIST_CH4::pressure(NIST_CH4::TT,D);
}

double NIST_CH4::density(double T, double P){
	double Ttemp = NIST_CH4::TT, Ptemp = NIST_CH4::PP, res;	
	NIST_CH4::TT=T;
	NIST_CH4::PP=P;	
	res= NR::zbrent(dpressure, 0.0005, 100.0, 1.0e-8); // mol/l
	NIST_CH4::TT = Ttemp;
	NIST_CH4::PP = Ptemp;
	return res;
}

double NIST_CH4::compressibility_factor(double T, double P){
	double tau,delta;
	delta= NIST_CH4::density(T,P)/Dc;	
	tau=Tc/T;
	return 1.0+delta*NR::dfridrY(alpha_r,tau,delta); //Eq.57
}

double NIST_CH4::volume(double T, double P){
	return 1.0/density(T,P);
}

double NIST_CH4::interal_energy(double T, double P){
	double tau,delta;
	delta= NIST_CH4::density(T,P)/Dc;	
	tau=Tc/T;
	return R*T*tau*(NR::dfridrX(alpha_0,tau,delta)+NR::dfridrX(alpha_r,tau,delta)) *M/1000 + 14.614;//ref to abs zero K
}

double NIST_CH4::enthalpy(double T, double P){
	double tau,delta,dr0,dr1,dr2;
	delta= NIST_CH4::density(T,P)/Dc;	
	tau=Tc/T;
	dr0=NR::dfridrX(alpha_0,tau,delta);
	dr1=NR::dfridrX(alpha_r,tau,delta);
	dr2=NR::dfridrY(alpha_r,tau,delta);
	return R*T*(tau*(dr0+dr1)+delta*dr2+1) *M/1000 + 14.614;
}

double NIST_CH4::entropy(double T, double P){
	double tau,delta;
	delta= NIST_CH4::density(T,P)/Dc;	
	tau=Tc/T;
	return M*R*(tau*(NR::dfridrX(alpha_0,tau,delta)+NR::dfridrX(alpha_r,tau,delta))-alpha_0(tau,delta)-alpha_r(tau,delta)) + 107.1166;//Eq.60 //ref to abs zero K
}

double NIST_CH4::gibbs_energy(double T, double P){
	double tau,delta;
	delta= NIST_CH4::density(T,P)/Dc;	
	tau=Tc/T;
	return R*T*(1.0+alpha_0(tau,delta)+alpha_r(tau,delta)+delta*(NR::dfridrY(alpha_r,tau,delta))) *M/1000 -17.3228;
}

double NIST_CH4::helmholtz_energy(double T, double P){
	double tau,delta;
	delta= NIST_CH4::density(T,P)/Dc;	
	tau=Tc/T;
	return R*T*(alpha_0(tau,delta)+alpha_r(tau,delta)) *M/1000 -17.3228;
}


double NIST_CH4::isochoric_heat_capacity(double T, double P){
	double tau,delta;
	delta= NIST_CH4::density(T,P)/Dc;	
	tau=Tc/T;
	return -M*R*tau*tau*(NR::dfridrX(dX_alpha_0,tau,delta)+NR::dfridrX(dX_alpha_r,tau,delta));//Eq.62
}

double NIST_CH4::isobaric_heat_capacity(double T, double P){
	double tau,delta,Cv,term1,term2;
	delta= NIST_CH4::density(T,P)/Dc;	
	tau=Tc/T;
	Cv=NIST_CH4::isochoric_heat_capacity(T,P);
	term1=1.0+delta*NR::dfridrY(alpha_r,tau,delta)-delta*tau*NR::dfridrY(dX_alpha_r,tau,delta);
	term1=term1*term1;
	term2=1.0+2.0*delta*NR::dfridrY(alpha_r,tau,delta)+delta*delta*NR::dfridrY(dY_alpha_r,tau,delta);
	return Cv+M*R*term1/term2;//Eq.63
}

double NIST_CH4::dQ(double dT){
	double H1,H2,T1,P1,T2,P2,dH,res;
	T1=NIST_CH4::TT;
	T2=NIST_CH4::TT+dT;
	P1=NIST_CH4::PP;
	P2=NIST_CH4::PP+NIST_CH4::dP;
	H1=NIST_CH4::enthalpy(T1,P1) ;
	H2=NIST_CH4::enthalpy(T2,P2);
	dH=H2-H1;
	res=dH-NIST_CH4::dP/NIST_CH4::DD;
	return res;
}

double NIST_CH4::adiabatic_temperature(double T, double P, double P1, int n){
	int i;
	double dT=0.,Tx,Px;
	NIST_CH4::TT=T;
	NIST_CH4::PP=P;
	NIST_CH4::dP=(P1-P)/n;
	NIST_CH4::DD=NIST_CH4::density(T,P);
	Tx=NIST_CH4::TT;
	Px=NIST_CH4::PP;
	if(fabs(P1-P)>1.0e-6)
		for(i=0;i<n;i++){
			if(P1>P) dT=NR::zbrent(dQ, 1.0e-6, 20.0, 1.0e-8);
			if(P1<P) dT=NR::zbrent(dQ, -1.0e-6, -20.0, 1.0e-8);
			Tx+=dT;
			Px+=NIST_CH4::dP;
			NIST_CH4::TT=Tx;
			NIST_CH4::PP=Px;
			NIST_CH4::DD=NIST_CH4::density(NIST_CH4::TT,NIST_CH4::PP);
		}
	return Tx;	
}

double NIST_CH4::adiabatic_temperature_id(double T, double P, double P1){
	double gamma,V,V1,Cv,Cp;
	Cp=35.31;//J K-1 mol-1, Table 2.6
	Cv=Cp-Rm;//Chapter.2 Atkins' Physical Chemistry, ninth edition
	gamma=Cp/Cv;
	V=Rm*T/P;//ideal gas
	V1=V*pow(P/P1,1/gamma); //P*V^gamma = P1*V1^gamma Eq(2.29) =>  V1/V=(P/P1)^(1/gamma)
	return P1*V1/Rm;
}


double NIST_CH4::isenthalpic_temperature(double T, double P, double P1, int n){
	int i;
	double dT=0., Tx, Px;
	double uJT = NIST_CH4::uJT(T, P);
	NIST_CH4::TT = T;
	NIST_CH4::PP = P;
	NIST_CH4::dP = (P1 - P) / n;
	Tx = NIST_CH4::TT;
	Px = NIST_CH4::PP;
	if (fabs(P1 - P)>1.0e-6)
	for (i = 0; i<n; i++){
		if ((P1 - P)*uJT>0)
			dT = NR::zbrent(dH, 1.0e-6, 5.0, 1.0e-8);
		if ((P1 - P)*uJT<0)
			dT = NR::zbrent(dH, -1.0e-6, -5.0, 1.0e-8);
		Tx += dT;
		Px += NIST_CH4::dP;
		NIST_CH4::TT = Tx;
		NIST_CH4::PP = Px;
	}
	return Tx;
}


double NIST_CH4::dH(double dT){
	double H1, H2, T1, P1, T2, P2, dH;
	T1 = NIST_CH4::TT;
	T2 = NIST_CH4::TT + dT;
	P1 = NIST_CH4::PP;
	P2 = NIST_CH4::PP + NIST_CH4::dP;
	H1 = NIST_CH4::enthalpy(T1, P1);
	H2 = NIST_CH4::enthalpy(T2, P2);
	dH = H2 - H1;
	return dH;
}
double NIST_CH4::uJT(double T, double P){
	double Cp, uT, res;
	Cp = NIST_CH4::isobaric_heat_capacity(T, P);
	uT = NR::dfridrY(NIST_CH4::enthalpy, T, P);
	res = -1.0e3*uT / Cp; //unit convert to K MPa-1
	return res;
}






//EoS for CO2, R. Span, W. Wagner, JPCRD, v25, 1509, 1996
double NIST_CO2::TT, NIST_CO2::PP, NIST_CO2::DD, NIST_CO2::dP;
double NIST_CO2::R, NIST_CO2::Rm;
double NIST_CO2::M, NIST_CO2::Tc, NIST_CO2::Pc, NIST_CO2::Dc;

void NIST_CO2::ReferenceConstants(void){ //Table 35
	R  = 0.1889241 ;//kJ kg-1 K-1
	Rm = 8.31451   ;//kJ kmol-1 K-1 (molar gas constant)
	M  = 44.0098   ;//kg kmol-1 (molecular weight)
	Tc = 304.1282  ;//K
	Pc = 7.3773    ;//MPa
	Dc = 467.6     ;//kg m-3
	Dc = Dc / M;    //mol dm-3
}

double NIST_CO2::alpha_0(double tau, double delta){
	const double a[9] = {0,8.37304456,-3.70454304,2.50000000,1.99427042,0.62105248,0.41195293,1.04028922,0.08327678};
	const double d[9] = {0,0,0,0,3.15163,6.11190,6.77708,11.32384,27.08792};
	int i;
	double res;
	res=log(delta)+a[1]+a[2]*tau+a[3]*log(tau);
	for(i=4;i<=8;i++) 
		res+=a[i]*log(1.0-exp(-d[i]*tau));
	return res;
}


double NIST_CO2::alpha_r(double tau, double delta){
	//Table 17, 18
	const double N[43]={
		 0.0                ,
		 0.38856823203161   ,
		 0.29385475942740e1 ,
		-0.55867188534934e1 ,
		-0.76753199592477   ,
		 0.31729005580416   ,
		 0.54803315897767   ,
		 0.12279411220335   ,
		 0.21658961543220e1 ,
		 0.15841735109724e1 ,
		-0.23132705405503   ,
		 0.58116916431436e-1,
		-0.55369137205382   ,
		 0.48946615909422   ,
		-0.24275739843501e-1,
		 0.62494790501678e-1,
		-0.12175860225246   ,
		-0.37055685270086   ,
		-0.16775879700426e-1,
		-0.11960736637987   ,
		-0.45619362508778e-1,
		 0.35612789270346e-1,
		-0.74427727132052e-2,
		-0.17395704902432e-2,
		-0.21810121289527e-1,
		 0.24332166559236e-1,
		-0.37440133423463e-1,
		 0.14338715756878   ,
		-0.13491969083286   ,
		-0.23151225053480e-1,
		 0.12363125492901e-1,
		 0.21058321972940e-2,
		-0.33958519026368e-3,
		 0.55993651771592e-2,
		-0.30335118055646e-3,
		-0.21365488688320e3 ,
		 0.26641569149272e5 ,
		-0.24027212204557e5 ,
		-0.28341603423999e3 ,
		 0.21247284400179e3 ,
		-0.66642276540751   ,
		 0.72608632349897   ,
		 0.55068668612842e-1
	};
	const double d[40]={0,1,1,1,1,2,2,3,1,2,4,5,5,5,6,6,6,1,1,4,4,4,7,8,2,3,3,5,5,6,7,8,10,4,8,2,2,2,3,3};
	const double t[40]={0.0,0.00,0.75,1.00,2.00,0.75,2.00,0.75,1.50,1.50,2.50,0.00,1.50,2.00,0.00,1.00,2.00,3.00,6.00,
		3.00,6.00,8.00,6.00,0.00,7.00,12.00,16.00,22.00,24.00,16.00,24.00,8.00,2.00,28.00,14.00,1.00,0.00,1.00,3.00,3.00};
	const double c[35]={0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,4,4,4,4,4,4,5,6};

	const double alfa[40] ={0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 25,25,25,15,20};
	const double beta[40] ={0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 325,300,300,275,275};
	const double gamma[40]={0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 1.16,1.19,1.19,1.25,1.22};

	const double a[43]={0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 3.5,3.5,3.0};
	const double b[43]={0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0.875,0.925,0.875};
	const double B[43]={0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0.3,0.3,1.0};
	const double C[43]={0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 10.0,10.0,12.5};
	
	double DELTA=0,res=0;
	int i;
	for(i=1;i<=7;i++)
		res+=N[i]*pow(delta,d[i])*pow(tau,t[i]);
	for(i=8;i<=34;i++)
		res+=N[i]*pow(delta,d[i])*pow(tau,t[i])*exp(-pow(delta,c[i]));
	for(i=35;i<=39;i++)
		res+=N[i]*pow(delta,d[i])*pow(tau,t[i]) *exp(-alfa[i]*(delta-1.0)*(delta-1.0)-beta[i]*(tau-gamma[i])*(tau-gamma[i]));
	for(i=40;i<=42;i++){
		DELTA=pow((1.0-tau+0.7*pow((delta-1.0)*(delta-1.0),1/0.6)),2) + B[i]*pow(((delta-1)*(delta-1)),a[i]);
		res+=N[i]*pow(DELTA,b[i])*delta*exp(-C[i]*(delta-1)*(delta-1)-275*(tau-1)*(tau-1));
	}
	return res;
}


double NIST_CO2::dX_alpha_0(double tau, double delta){ return NR::dfridrX(alpha_0,tau,delta);}
double NIST_CO2::dY_alpha_0(double tau, double delta){ return NR::dfridrY(alpha_0,tau,delta);}
double NIST_CO2::dX_alpha_r(double tau, double delta){ return NR::dfridrX(alpha_r,tau,delta);}
double NIST_CO2::dY_alpha_r(double tau, double delta){ return NR::dfridrY(alpha_r,tau,delta);}


double NIST_CO2::pressure(double T, double D){
	NIST_CO2::ReferenceConstants();
	double tau,delta;
	delta=D/Dc;
	tau=Tc/T;
	return 1.0e-3*D*M*R*T*(1+delta*NR::dfridrY(alpha_r,tau,delta)); 
	//Eq. 56  MPa unit of pressure J.Phys.Chem.Ref.Data, vol.29, No.6,2000
}

double NIST_CO2::dpressure(double D){ 
	double res;
	res = NIST_CO2::PP - NIST_CO2::pressure(NIST_CO2::TT, D);
	if (NIST_CO2::pressure(NIST_CO2::TT, D) <= 0)
		res = NIST_CO2::PP;
	return res;
	//NIST_CO2::PP - NIST_CO2::pressure(NIST_CO2::TT, D);
}

double NIST_CO2::density(double T, double P){
	double Ttemp = NIST_CO2::TT, Ptemp = NIST_CO2::PP, res=0.;
	double Ps, DV, DL;
	NIST_CO2::ReferenceConstants();
	NIST_CO2::TT=T;
	NIST_CO2::PP=P;
	if (T>Tc+1)
		res = NR::zbrent(dpressure, 1.0e-6, 200, 1.0e-8);	
	else{
		Ps = NIST_CO2::vapor_pressure(T);
		DV = NIST_CO2::saturated_vapor_density(T);
		DL = NIST_CO2::saturated_liquid_density(T);
		if (P>=Ps)
			res = NR::zbrent(dpressure, DL*0.9, 200, 1.0e-8);
		else if (P<Ps)
			res = NR::zbrent(dpressure, 0.0001, DV*1.2, 1.0e-8);
	}
	NIST_CO2::TT = Ttemp;
	NIST_CO2::PP = Ptemp;
	return res;
}

double NIST_CO2::volume(double T, double P){
	return 1.0/density(T,P);
}

double NIST_CO2::compressibility_factor(double T, double P){
	double tau,delta;
	delta= NIST_CO2::density(T,P)/Dc;
	tau=Tc/T;
	return 1.0+delta*NR::dfridrY(alpha_r,tau,delta); //Eq.57
}

double NIST_CO2::interal_energy(double T, double P){
	double tau,delta;
	delta= NIST_CO2::density(T,P)/Dc;
	tau=Tc/T;
	return R*T*tau*(NR::dfridrX(alpha_0,tau,delta)+NR::dfridrX(alpha_r,tau,delta)) *M/1000 + 22.3036; //Eq.58
}

double NIST_CO2::enthalpy(double T, double P){
	double tau,delta,dr0,dr1,dr2;
	delta= NIST_CO2::density(T,P)/Dc;
	tau=Tc/T;	
	dr0=NR::dfridrX(alpha_0,tau,delta);
	dr1=NR::dfridrX(alpha_r,tau,delta);
	dr2=NR::dfridrY(alpha_r,tau,delta);
	return R*T*(tau*(dr0+dr1)+delta*dr2+1) *M/1000 + 22.3036;
}

double NIST_CO2::entropy(double T, double P){
	double tau,delta;
	delta= NIST_CO2::density(T,P)/Dc;
	tau=Tc/T;
	return R*(tau*(NR::dfridrX(alpha_0,tau,delta)+NR::dfridrX(alpha_r,tau,delta))-alpha_0(tau,delta)-alpha_r(tau,delta)) *M + 120.5441;//Eq.60
}

double NIST_CO2::gibbs_energy(double T, double P){
	double tau,delta;
	delta= NIST_CO2::density(T,P)/Dc;
	tau=Tc/T;
	return R*T*(1.0+alpha_0(tau,delta)+alpha_r(tau,delta)+delta*(NR::dfridrY(alpha_r,tau,delta))) *M/1000 -13.85963;//Eq.61
}

double NIST_CO2::helmholtz_energy(double T, double P){
	double tau,delta;
	delta= NIST_CO2::density(T,P)/Dc;
	tau=Tc/T;
	return R*T*(alpha_0(tau,delta)+alpha_r(tau,delta)) *M/1000 -13.85963;
}


double NIST_CO2::isochoric_heat_capacity(double T, double P){
	double tau,delta;
	delta= NIST_CO2::density(T,P)/Dc;
	tau=Tc/T;
	return -R*tau*tau*(NR::dfridrX(dX_alpha_0,tau,delta)+NR::dfridrX(dX_alpha_r,tau,delta)) *M;//Eq.62
}

double NIST_CO2::isobaric_heat_capacity(double T, double P){
	double tau,delta,Cv,term1,term2;
	delta= NIST_CO2::density(T,P)/Dc;
	tau=Tc/T;
	Cv=NIST_CO2::isochoric_heat_capacity(T,P);
	term1=1.0+delta*NR::dfridrY(alpha_r,tau,delta)-delta*tau*NR::dfridrY(dX_alpha_r,tau,delta);
	term1=term1*term1;
	term2=1.0+2.0*delta*NR::dfridrY(alpha_r,tau,delta)+delta*delta*NR::dfridrY(dY_alpha_r,tau,delta);
	return Cv+R*term1/term2*M;//Eq.63
}

double NIST_CO2::dQ(double dT){
	double H1,H2,T1,P1,T2,P2,dH,res;
	T1=NIST_CO2::TT;
	T2=NIST_CO2::TT+dT;
	P1=NIST_CO2::PP;
	P2=NIST_CO2::PP+NIST_CO2::dP;
	H1=NIST_CO2::enthalpy(T1,P1) ;
	H2=NIST_CO2::enthalpy(T2,P2);
	dH=H2-H1;
	res=dH-NIST_CO2::dP/NIST_CO2::DD;
	return res;
}

double NIST_CO2::adiabatic_temperature(double T, double P, double P1, int n){
	int i;
	double dT=0.,Tx,Px;
	NIST_CO2::TT=T;
	NIST_CO2::PP=P;
	NIST_CO2::dP=(P1-P)/n;
	NIST_CO2::DD=NIST_CO2::density(T,P);
	Tx=NIST_CO2::TT;
	Px=NIST_CO2::PP;
	if(fabs(P1-P)>1.0e-6)
		for(i=0;i<n;i++){
			if(P1>P) dT=NR::zbrent(dQ, 1.0e-6, 20.0, 1.0e-8);
			if(P1<P) dT=NR::zbrent(dQ, -1.0e-6, -20.0, 1.0e-8);
			Tx+=dT;
			Px+=NIST_CO2::dP;
			NIST_CO2::TT=Tx;
			NIST_CO2::PP=Px;
			NIST_CO2::DD=NIST_CO2::density(NIST_CO2::TT,NIST_CO2::PP);
		}
	return Tx;	
}

double NIST_CO2::adiabatic_temperature_id(double T, double P, double P1){
	double gamma,V,V1,Cv,Cp;
	Cp=37.11;//J K-1 mol-1, Table 2.6
	Cv=Cp-Rm;//Chapter.2 Atkins' Physical Chemistry, ninth edition
	gamma=Cp/Cv;
	V=Rm*T/P;//ideal gas
	V1=V*pow(P/P1,1/gamma); //P*V^gamma = P1*V1^gamma Eq(2.29) =>  V1/V=(P/P1)^(1/gamma)
	return P1*V1/Rm;
}

double NIST_CO2::vapor_pressure(double T){
	double res=0.0;
	const double a[4] = {-7.0602087, 1.9391218, -1.6463597, -3.2995634 };
	const double t[4] = {1.0, 1.5, 2.0, 4.0};
	int i;
	NIST_CO2::ReferenceConstants();
	if (T > Tc)
		T = Tc;
	for (i = 0; i < 4; i++){
		res = res + a[i]*pow((1.0 - T / NIST_CO2::Tc), t[i]);
	}
	res =  exp(NIST_CO2::Tc/T*res)*NIST_CO2::Pc;
	return res;
}

double NIST_CO2::saturated_liquid_density(double T){
	double res = 0.0;
	const double a[4] = { 1.9245108, -0.62385555, -0.32731127, 0.39245142 };
	double t[4];
	t[0] = 0.34;
	t[1] = 0.5;
	t[2] = 10 / 6;
	t[3] = 11 / 6;
	int i;
	NIST_CO2::ReferenceConstants();
	for (i = 0; i < 4; i++)
		res = res+ a[i] * pow((1.0 - T / NIST_CO2::Tc), t[i]);
	res = exp(res)*NIST_CO2::Dc;
	return res;
}

double NIST_CO2::saturated_vapor_density(double T){
	double res = 0.0;
	const double a[5] = { -1.7074879, -0.82274670, -4.6008549, -10.111178, -29.742252 };
	double t[5];
	t[0] = 0.34;
	t[1] = 0.5;
	t[2] = 1.0;
	t[3] = 7 / 3;
	t[4] = 14 / 3;
	int i;
	NIST_CO2::ReferenceConstants();
	for (i = 0; i < 5; i++)
		res = res+ a[i] * pow((1.0 - T / NIST_CO2::Tc), t[i]);
	res = exp(res)*NIST_CO2::Dc;
	return res;
}

double NIST_CO2::isenthalpic_temperature(double T, double P, double P1, int n){
	int i;
	double dT=0., Tx, Px;
	double uJT = NIST_CO2::uJT(T, P);
	NIST_CO2::TT = T;
	NIST_CO2::PP = P;
	NIST_CO2::dP = (P1 - P) / n;
	Tx = NIST_CO2::TT;
	Px = NIST_CO2::PP;
	if (fabs(P1 - P)>1.0e-6)
	for (i = 0; i<n; i++){
		if ((P1 - P)*uJT>0)
			dT = NR::zbrent(dH, 1.0e-6, 5.0, 1.0e-8);
		if ((P1 - P)*uJT<0)
			dT = NR::zbrent(dH, -1.0e-6, -5.0, 1.0e-8);
		Tx += dT;
		Px += NIST_CO2::dP;
		NIST_CO2::TT = Tx;
		NIST_CO2::PP = Px;
	}
	return Tx;
}


double NIST_CO2::dH(double dT){
	double H1, H2, T1, P1, T2, P2, dH;
	T1 = NIST_CO2::TT;
	T2 = NIST_CO2::TT + dT;
	P1 = NIST_CO2::PP;
	P2 = NIST_CO2::PP + NIST_CO2::dP;
	H1 = NIST_CO2::enthalpy(T1, P1);
	H2 = NIST_CO2::enthalpy(T2, P2);
	dH = H2 - H1;
	return dH;
}
double NIST_CO2::uJT(double T, double P){
	double Cp, uT, res;
	Cp = NIST_CO2::isobaric_heat_capacity(T, P);
	uT = NR::dfridrY(NIST_CO2::enthalpy, T, P);
	res = -1.0e3*uT / Cp; //unit convert to K MPa-1
	return res;
}





//R. Schmidt, W. Wagner, 1985 Fluid Phase Equilibria, 19:175-200
double NIST_O2::TT, NIST_O2::PP, NIST_O2::DD, NIST_O2::dP;
double NIST_O2::R, NIST_O2::Rm;
double NIST_O2::M, NIST_O2::Tc, NIST_O2::Pc, NIST_O2::Dc;

void NIST_O2::ReferenceConstants(void){//Table 2 and Table 5
	Rm = 8.31434 ;//kJ kmol-1 K-1 (molar gas constant)
	M  = 32.00   ;//kg kmol-1 (molecular weight)
	Tc = 154.581 ;//K
	Pc = 50.43   ;//bar
	Dc = 13.63   ;//mol dm-3
}

double NIST_O2::alpha_0(double tau, double delta){
	//Table 3
	const double k[10] = {0,-0.740775e-3,-0.664930e-4, 0.250042e1,-0.214487e2, 0.101258e1,-0.944365, 0.145066e2, 0.749148e2, 0.414817e1};
	double delta0,T0,p0;
	double res;
	T0=298.15;
	p0=1.10325;
	delta0=p0/(Rm*T0*Dc)*1.0e2;
	res=k[1]*pow(tau, 1.5)+k[2]*pow(tau,-2)+k[3]*log(tau)+k[4]*tau+k[5]*log(exp(k[7]*tau)-1)
		+k[6]*log(1+2/3*exp(-k[8]*tau))+k[9]+log(delta/delta0); //Eq(15)
	return res;
}


double NIST_O2::alpha_r(double tau, double delta){
	//Table 2
	const double N[33]={
		 0.0            ,
		 0.3983768749   ,
		-0.1846157454e1 ,
		 0.4183473197   ,
		 0.2370620711e-1,
		 0.9771730573e-1,
		 0.3017891294e-1,
		 0.2273353212e-1,
		 0.1357254086e-1,
		-0.4052698943e-1,
		 0.5454628515e-3,
		 0.5113182277e-3,
		 0.2953466883e-6,
		-0.8687645072e-4,
		-0.2127082589   ,
		 0.8735941958e-1,
		 0.1275509190   ,
		-0.9067701064e-1,
		-0.3540084206e-1,
		-0.3623278059e-1,
		 0.1327699290e-1,
		-0.3254111865e-3,
		-0.8313582932e-2,
		 0.2124570559e-2,
		-0.8325206232e-3,
		-0.2626173276e-4,
		 0.2599581482e-2,
		 0.9984649663e-2,
		 0.2199923153e-2,
		-0.2591350486e-1,
		-0.1259630848   ,
		 0.1478355637   ,
		-0.1011251078e-1};
	const double r[33]={0,1,1,1,2,2,2,3,3,3,6,7,7,8,1,1,2,2,3,3,5,6,7,8,10,2,3,3,4,4,5,5,5};
	const double s[33]={0,0,1.5,2.5,-0.5,1.5,2,0,1,2.5,0,2,5,2,5,6,3.5,5.5,3,7,6,8.5,4,6.5,5.5,22,11,18,11,23,17,18,23};
	double res=0;
	int i;
	for(i=1;i<=13;i++) //Eq(11)
		res+=N[i]*pow(delta,r[i])*pow(tau,s[i]);
	for(i=14;i<=24;i++)
		res+=N[i]*pow(delta,r[i])*pow(tau,s[i])*exp(-delta*delta);
	for(i=25;i<=32;i++)
		res+=N[i]*pow(delta,r[i])*pow(tau,s[i])*exp(-pow(delta,4));
	return res;
}



double NIST_O2::dX_alpha_0(double tau, double delta){ return NR::dfridrX(alpha_0,tau,delta);}
double NIST_O2::dY_alpha_0(double tau, double delta){ return NR::dfridrY(alpha_0,tau,delta);}
double NIST_O2::dX_alpha_r(double tau, double delta){ return NR::dfridrX(alpha_r,tau,delta);}
double NIST_O2::dY_alpha_r(double tau, double delta){ return NR::dfridrY(alpha_r,tau,delta);}


double NIST_O2::pressure(double T, double D){
	NIST_O2::ReferenceConstants();
	double tau,delta;
	delta=D/Dc;	
	tau=Tc/T;
	return 1.0e-3*D*Rm*T*(1+delta*NR::dfridrY(alpha_r,tau,delta)); 
}

double NIST_O2::dpressure(double D){ 
	return NIST_O2::PP-NIST_O2::pressure(NIST_O2::TT,D);
}

double NIST_O2::density(double T, double P){
	double Ttemp = NIST_O2::TT, Ptemp = NIST_O2::PP, res;
	NIST_O2::TT=T;
	NIST_O2::PP=P;
	if (P<200)	
		res= NR::zbrent(dpressure, 0.0001, 50.0, 1.0e-8);
	else 
		res = NR::zbrent(dpressure, 0.0005, 200.0, 1.0e-8);
	NIST_O2::TT = Ttemp;
	NIST_O2::PP = Ptemp;
	return res;
}

double NIST_O2::volume(double T, double P){
	return 1.0/density(T,P);
}


double NIST_O2::compressibility_factor(double T, double P){
	double tau,delta;
	delta= NIST_O2::density(T,P)/Dc;	
	tau=Tc/T;
	return 1.0+delta*NR::dfridrY(alpha_r,tau,delta); //Eq.57
}

double NIST_O2::interal_energy(double T, double P){
	double tau,delta;
	delta= NIST_O2::density(T,P)/Dc;	
	tau=Tc/T;
	return 1.0e-3*Rm*T*tau*(NR::dfridrX(alpha_0,tau,delta)+NR::dfridrX(alpha_r,tau,delta)) + 8.68; //Eq.58
}

double NIST_O2::enthalpy(double T, double P){
	double tau,delta,dr0,dr1,dr2;
	delta= NIST_O2::density(T,P)/Dc;	
	tau=Tc/T;
	dr0=NR::dfridrX(alpha_0,tau,delta);
	dr1=NR::dfridrX(alpha_r,tau,delta);
	dr2=NR::dfridrY(alpha_r,tau,delta);
	return 1.0e-3*Rm*T*(tau*(dr0+dr1)+delta*dr2+1)  + 8.68;
}

double NIST_O2::entropy(double T, double P){
	double tau,delta;
	delta= NIST_O2::density(T,P)/Dc;	
	tau=Tc/T;
	return Rm*(tau*(NR::dfridrX(alpha_0,tau,delta)+NR::dfridrX(alpha_r,tau,delta))-alpha_0(tau,delta)-alpha_r(tau,delta)) + 204.3322;//Eq.60
}

double NIST_O2::gibbs_energy(double T, double P){
	double tau,delta;
	delta= NIST_O2::density(T,P)/Dc;	
	tau=Tc/T;
	return 1.0e-3*Rm*T*(1.0+alpha_0(tau,delta)+alpha_r(tau,delta)+delta*(NR::dfridrY(alpha_r,tau,delta))) -52.61966;//Eq.61
}

double NIST_O2::helmholtz_energy(double T, double P){
	double tau,delta;
	delta= NIST_O2::density(T,P)/Dc;	
	tau=Tc/T;
	return 1.0e-3*Rm*T*(alpha_0(tau,delta)+alpha_r(tau,delta))  -52.61966;
}


double NIST_O2::isochoric_heat_capacity(double T, double P){
	double tau,delta;
	delta= NIST_O2::density(T,P)/Dc;	
	tau=Tc/T;
	return -Rm*tau*tau*(NR::dfridrX(dX_alpha_0,tau,delta)+NR::dfridrX(dX_alpha_r,tau,delta));//Eq.62
}

double NIST_O2::isobaric_heat_capacity(double T, double P){
	double tau,delta,Cv,term1,term2;
	delta= NIST_O2::density(T,P)/Dc;	
	tau=Tc/T;
	Cv=NIST_O2::isochoric_heat_capacity(T,P);
	term1=1.0+delta*NR::dfridrY(alpha_r,tau,delta)-delta*tau*NR::dfridrY(dX_alpha_r,tau,delta);
	term1=term1*term1;
	term2=1.0+2.0*delta*NR::dfridrY(alpha_r,tau,delta)+delta*delta*NR::dfridrY(dY_alpha_r,tau,delta);
	return Cv+Rm*term1/term2;//Eq.63
}

double NIST_O2::dQ(double dT){
	double H1,H2,T1,P1,T2,P2,dH,res;
	T1=NIST_O2::TT;
	T2=NIST_O2::TT+dT;
	P1=NIST_O2::PP;
	P2=NIST_O2::PP+NIST_O2::dP;
	H1=NIST_O2::enthalpy(T1,P1) ;
	H2=NIST_O2::enthalpy(T2,P2);
	dH=H2-H1;
	res=dH-NIST_O2::dP/NIST_O2::DD;
	return res;
}

double NIST_O2::adiabatic_temperature(double T, double P, double P1, int n){
	int i;
	double dT=0.,Tx,Px;
	NIST_O2::TT=T;
	NIST_O2::PP=P;
	NIST_O2::dP=(P1-P)/n;
	NIST_O2::DD=NIST_O2::density(T,P);
	Tx=NIST_O2::TT;
	Px=NIST_O2::PP;
	if(fabs(P1-P)>1.0e-6)
		for(i=0;i<n;i++){
			if(P1>P) dT=NR::zbrent(dQ, 1.0e-6, 20.0, 1.0e-8);
			if(P1<P) dT=NR::zbrent(dQ, -1.0e-6, -20.0, 1.0e-8);
			Tx+=dT;
			Px+=NIST_O2::dP;
			NIST_O2::TT=Tx;
			NIST_O2::PP=Px;
			NIST_O2::DD=NIST_O2::density(NIST_O2::TT,NIST_O2::PP);
		}
	return Tx;	
}

double NIST_O2::adiabatic_temperature_id(double T, double P, double P1){
	double gamma,V,V1,Cv,Cp;
    Cp=29.355;//J K-1 mol-1, Table 2.8
	Cv=Cp-Rm;//Chapter.2 Atkins' Physical Chemistry, ninth edition
	gamma=Cp/Cv;
	V=Rm*T/P;//ideal gas
	V1=V*pow(P/P1,1/gamma); //P*V^gamma = P1*V1^gamma Eq(2.29) =>  V1/V=(P/P1)^(1/gamma)
	return P1*V1/Rm;
}

double NIST_O2::isenthalpic_temperature(double T, double P, double P1, int n){
	int i;
	double dT=0., Tx, Px;
	double uJT = NIST_O2::uJT(T, P);
	NIST_O2::TT = T;
	NIST_O2::PP = P;
	NIST_O2::dP = (P1 - P) / n;
	Tx = NIST_O2::TT;
	Px = NIST_O2::PP;
	if (fabs(P1 - P)>1.0e-6)
	for (i = 0; i<n; i++){
		if ((P1 - P)*uJT>0)
			dT = NR::zbrent(dH, 1.0e-6, 5.0, 1.0e-8);
		if ((P1 - P)*uJT<0)
			dT = NR::zbrent(dH, -1.0e-6, -5.0, 1.0e-8);
		Tx += dT;
		Px += NIST_O2::dP;
		NIST_O2::TT = Tx;
		NIST_O2::PP = Px;
	}
	return Tx;
}


double NIST_O2::dH(double dT){
	double H1, H2, T1, P1, T2, P2, dH;
	T1 = NIST_O2::TT;
	T2 = NIST_O2::TT + dT;
	P1 = NIST_O2::PP;
	P2 = NIST_O2::PP + NIST_O2::dP;
	H1 = NIST_O2::enthalpy(T1, P1);
	H2 = NIST_O2::enthalpy(T2, P2);
	dH = H2 - H1;
	return dH;
}
double NIST_O2::uJT(double T, double P){
	double Cp, uT, res;
	Cp = NIST_O2::isobaric_heat_capacity(T, P);
	uT = NR::dfridrY(NIST_O2::enthalpy, T, P);
	res = -1.0e3*uT / Cp; //unit convert to K MPa-1
	return res;
}




//for ideal mixture system of the 5 gases componts
double NIST_M::TT, NIST_M::PP, NIST_M::VV, NIST_M::dP;
double NIST_M::mH2, NIST_M::mN2, NIST_M::mO2, NIST_M::mCH4, NIST_M::mCO2;

double NIST_M::enthalpy(double T, double P, double mH2, double mN2, double mO2, double mCH4, double mCO2){
	double res;
	res = 0;
	res += mH2*NIST_H2::enthalpy(T, P);
	res += mN2*NIST_N2::enthalpy(T, P);
	res += mO2*NIST_O2::enthalpy(T, P);
	res += mCH4*NIST_CH4::enthalpy(T, P);
	res += mCO2*NIST_CO2::enthalpy(T, P);
	return res;
}

double NIST_M::volume(double T, double P, double mH2, double mN2, double mO2, double mCH4, double mCO2){
	double res;
	res = 0;
	res += mH2*NIST_H2::volume(T, P);
	res += mN2*NIST_N2::volume(T, P);
	res += mO2*NIST_O2::volume(T, P);
	res += mCH4*NIST_CH4::volume(T, P);
	res += mCO2*NIST_CO2::volume(T, P);
	return res;
}

double NIST_M::interal_energy(double T, double P, double mH2, double mN2, double mO2, double mCH4, double mCO2){
	double res;
	res = 0;
	res += mH2*NIST_H2::interal_energy(T, P);
	res += mN2*NIST_N2::interal_energy(T, P);
	res += mO2*NIST_O2::interal_energy(T, P);
	res += mCH4*NIST_CH4::interal_energy(T, P);
	res += mCO2*NIST_CO2::interal_energy(T, P);
	return res;
}

double NIST_M::gibbs_energy(double T, double P, double mH2, double mN2, double mO2, double mCH4, double mCO2){
	double res;
	res = 0;
	res += mH2*NIST_H2::gibbs_energy(T, P);
	res += mN2*NIST_N2::gibbs_energy(T, P);
	res += mO2*NIST_O2::gibbs_energy(T, P);
	res += mCH4*NIST_CH4::gibbs_energy(T, P);
	res += mCO2*NIST_CO2::gibbs_energy(T, P);
	return res;
}

double NIST_M::helmholtz_energy(double T, double P, double mH2, double mN2, double mO2, double mCH4, double mCO2){
	double res;
	res = 0;
	res += mH2*NIST_H2::helmholtz_energy(T, P);
	res += mN2*NIST_N2::helmholtz_energy(T, P);
	res += mO2*NIST_O2::helmholtz_energy(T, P);
	res += mCH4*NIST_CH4::helmholtz_energy(T, P);
	res += mCO2*NIST_CO2::helmholtz_energy(T, P);
	return res;
}

double NIST_M::density(double T, double P, double mH2, double mN2, double mO2, double mCH4, double mCO2){
	double total_mass,total_volume,res;
    NIST_H2::ReferenceConstants();
    NIST_N2::ReferenceConstants();
    NIST_O2::ReferenceConstants();
    NIST_CH4::ReferenceConstants();
    NIST_CO2::ReferenceConstants();
	total_mass = 0;
	total_mass += mH2*NIST_H2::M;
	total_mass += mN2*NIST_N2::M;
	total_mass += mO2*NIST_O2::M;
	total_mass += mCH4*NIST_CH4::M;
	total_mass += mCO2*NIST_CO2::M;
	total_volume = NIST_M::volume(T, P, mH2, mN2, mO2, mCH4, mCO2);
	res = total_mass / total_volume;
	return res;
}

double NIST_M::entropy(double T, double P, double mH2, double mN2, double mO2, double mCH4, double mCO2){
	double res,mM,P1,P2,P3,P4,P5;
	mM = mH2 + mN2 + mO2 + mCH4 + mCO2;
	P1 = mH2 / mM;
	P2 = mN2 / mM;
	P3 = mO2 / mM;
	P4 = mCH4 / mM;
	P5 = mCO2 / mM;	
	res = 0;
	res += mH2*NIST_H2::entropy(T, P1);
	res += mN2*NIST_N2::entropy(T, P2);
	res += mO2*NIST_O2::entropy(T, P3);
	res += mCH4*NIST_CH4::entropy(T, P4);
	res += mCO2*NIST_CO2::entropy(T, P5);
	return res;
}

double NIST_M::isochoric_heat_capacity(double T, double P, double mH2, double mN2, double mO2, double mCH4, double mCO2){
	double res;
	res = 0;
	res += mH2*NIST_H2::isochoric_heat_capacity(T, P);
	res += mN2*NIST_N2::isochoric_heat_capacity(T, P);
	res += mO2*NIST_O2::isochoric_heat_capacity(T, P);
	res += mCH4*NIST_CH4::isochoric_heat_capacity(T, P);
	res += mCO2*NIST_CO2::isochoric_heat_capacity(T, P);
	return res;
}

double NIST_M::isobaric_heat_capacity(double T, double P, double mH2, double mN2, double mO2, double mCH4, double mCO2){
	double res;
	res = 0;
	res += mH2*NIST_H2::isobaric_heat_capacity(T, P);
	res += mN2*NIST_N2::isobaric_heat_capacity(T, P);
	res += mO2*NIST_O2::isobaric_heat_capacity(T, P);
	res += mCH4*NIST_CH4::isobaric_heat_capacity(T, P);
	res += mCO2*NIST_CO2::isobaric_heat_capacity(T, P);
	return res;
}

double NIST_M::dQ(double dT){
	double H1, H2, T1, P1, T2, P2, dH, res;
	T1 = NIST_M::TT;
	T2 = NIST_M::TT + dT;
	P1 = NIST_M::PP;
	P2 = NIST_M::PP + NIST_M::dP;
	H1 = NIST_M::enthalpy(T1, P1, NIST_M::mH2, NIST_M::mN2, NIST_M::mO2, NIST_M::mCH4, NIST_M::mCO2);
	H2 = NIST_M::enthalpy(T2, P2, NIST_M::mH2, NIST_M::mN2, NIST_M::mO2, NIST_M::mCH4, NIST_M::mCO2);
	dH = H2 - H1;
	res = dH - NIST_M::dP * NIST_M::VV;
	return res;
}

double NIST_M::dH(double dT){
	double H1, H2, T1, P1, T2, P2, dH;
	T1 = NIST_M::TT;
	T2 = NIST_M::TT + dT;
	P1 = NIST_M::PP;
	P2 = NIST_M::PP + NIST_M::dP;
	H1 = NIST_M::enthalpy(T1, P1, NIST_M::mH2, NIST_M::mN2, NIST_M::mO2, NIST_M::mCH4, NIST_M::mCO2);
	H2 = NIST_M::enthalpy(T2, P2, NIST_M::mH2, NIST_M::mN2, NIST_M::mO2, NIST_M::mCH4, NIST_M::mCO2);
	dH = H2 - H1;
	return dH;
}

double NIST_M::adiabatic_temperature(double T0, double P0, double P1, int n, double mH2, double mN2, double mO2, double mCH4, double mCO2){
	int i;
	double dT=0., Tx, Px;
	NIST_M::mH2 = mH2;
	NIST_M::mN2 = mN2;
	NIST_M::mO2 = mO2;
	NIST_M::mCH4 = mCH4;
	NIST_M::mCO2 = mCO2;
	NIST_M::TT = T0;
	NIST_M::PP = P0;
	NIST_M::dP = (P1 - P0) / n;
	NIST_M::VV = NIST_M::volume(T0, P0, mH2, mN2, mO2, mCH4, mCO2);
	Tx = NIST_M::TT;
	Px = NIST_M::PP;
	if (fabs(P1 - P0)>1.0e-6)
	for (i = 0; i<n; i++){
		if (P1>P0) dT = NR::zbrent(dQ, 1.0e-6, 20.0, 1.0e-8);
		if (P1<P0) dT = NR::zbrent(dQ, -1.0e-6, -20.0, 1.0e-8);
		Tx += dT;
		Px += NIST_M::dP;
		NIST_M::TT = Tx;
		NIST_M::PP = Px;
		NIST_M::VV = NIST_M::volume(NIST_M::TT, NIST_M::PP, mH2,mN2,mO2,mCH4,mCO2);
	}
	return Tx-T0;
}

double NIST_M::isenthalpic_temperature(double T0, double P0, double P1, int n, double mH2, double mN2, double mO2, double mCH4, double mCO2){
	int i;
	double dT=0., Tx, Px;
	NIST_M::mH2 = mH2;
	NIST_M::mN2 = mN2;
	NIST_M::mO2 = mO2;
	NIST_M::mCH4 = mCH4;
	NIST_M::mCO2 = mCO2;
	//double uJT = NIST_H2::uJT(T, P);
	NIST_M::TT = T0;
	NIST_M::PP = P0;
	NIST_M::dP = (P1 - P0) / n;
	Tx = NIST_M::TT;
	Px = NIST_M::PP;
	if (fabs(P1 - P0)>1.0e-6)
	for (i = 0; i<n; i++){
		if (dH(1.0e-6)*dH(5.0) < 0)
			dT = NR::zbrent(dH, 1.0e-6, 5.0, 1.0e-8);
		else
			dT = NR::zbrent(dH, -1.0e-6, -5.0, 1.0e-8);
		Tx += dT;
		Px += NIST_M::dP;
		NIST_M::TT = Tx;
		NIST_M::PP = Px;
	}
	return Tx-T0;
}


