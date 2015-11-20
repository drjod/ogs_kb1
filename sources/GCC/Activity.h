#include <string>
#include <vector>
using namespace std;
class Activity
{
private:

public:
	Activity(void);
	~Activity(void);

//data
	static int activity_coefficient(int f, double T, double P, vector<string> ListName, vector<double> ListMole, vector<double> &ListLGAMMA, double &AW);
	//f - flag for activity coefficients methods
	//1 - D-H limiting law
	//2 - Revised D-H equation
	//3 - Guggenheim equation
	//4 - Davies equation
	//5 - Extend D-H equation
	//6 - Truesdell-Jones equation (T-J)
	//7 - Specific ion interaction theory (SIT)
	//8 - Bromley equation
	//9 - Helgeson equation
	//10- Pitzer equation
	//11- B-dot
	static double  Ar(double T, double P);
	static double  Br(double T, double P);
	static void entrance(void);
};
