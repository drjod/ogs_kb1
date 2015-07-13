#include <string> 

//using namespace std;

class PITZdata
{
private:

public:
	PITZdata(void);
	~PITZdata(void);
//method
	static double charge(std::string N);
	static double pitzer_parameters(double T, double P, std::string param_switch);

};

