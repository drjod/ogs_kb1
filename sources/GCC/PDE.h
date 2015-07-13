//using namespace std; 
class PDE
{
private:
public:
	PDE(void);
	~PDE(void);
/* Data */
	static double TT,PP,d0,d1,d2,d3,d4;

/* Methods */
	static void diffusion_CO2(void);
	static double diffusion_coefficient_CO2(double C);
	static void entrance(void);
};



