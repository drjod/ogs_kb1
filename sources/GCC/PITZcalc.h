#include <string> 
#include <vector>
//using namespace std;
class PITZcalc
{
private:

public:
	PITZcalc(void);
	~PITZcalc(void);

//data

	static double  Z[40];
	static double  M[40], MN[20], LGAMMA[40], LGN[20], AW;	
	static    int  M1, M2, M3;
	static double  TT, PP;

	static double  BCX[4][20][40];
	//COMMON / MX3 / BCX(4,20,21:40)
	static double  THETA[40][40],  PSI[40][40][40];
	//COMMON / MX7 / THETA(40,40),PSI(40,40,40) (++- or --+)
	static double  LAM[40][20], LAMN[20][20], ZETA[20][40][20];
	//COMMON / MX6 / LAM(40,20),LAMN(20,20),ZETA(20,21:40,20)

	static double  I;
	static double  A0;
	static double  ALPHA[3];
	static double  BK[23], DK[23];

//method
	static int pitzer(void);
	static int activity_coefficient(double T, double P, std::vector<std::string> ListName, std::vector<double> ListMole, std::vector<double> &ListLGAMMA, double &ac);


	static double  D(double T, double P);


	static double  APHI(double T, double P);
	static double  PHIPHI(int J, int K);
	static double  PHI(int J, int K);
	static double  PHIP(int J, int K);
	static    int  BDK(double X);
	static double  JAY(double X);
	static double  JPRIME(double Y);
	static double  G(double Y);
	static double  GP(double Y);	
	static double  ETHETA(double ZJ, double ZK);
	static double  ETHETAP(double ZJ, double ZK);
	static   bool  TWOTWO(int J, int K);
	static double  BMXPHI(int J, int K);
	static double  BMX(int J, int K);
	static double  BMXP(int J, int K);
	static double  CMX(int J, int K);


	static void entrance(void);
};


