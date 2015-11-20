class NR 
{
private:

public:
	NR(void);
	~NR(void);
/* Methods */
	static double zbrent(double func(double), const double x1, const double x2, const double tol);
	static inline double SIGN(const double &a, const double &b);

	static double dfridr (double func(double),        const double x);
	static double dfridrX(double func(double,double), const double x, const double y);
	static double dfridrY(double func(double,double), const double x, const double y);

};

