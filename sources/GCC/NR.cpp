#include <iostream> 
#include <math.h>
#include <cmath>
#include <limits>
#include <stdio.h>
#include <string>
#include "NR.h"
using namespace std;


NR::NR(void){}
NR::~NR(void){}

//numerical recipe 2nd or 3rd edition
//numerical recipes 3rd edition section 5.7
double NR::zbrent(double func(double), const double x1, const double x2, const double tol)
{
	const int ITMAX=100;
	const double EPS=numeric_limits<double>::epsilon();
	double a=x1,b=x2,c=x2,d,e=0.0,fa=func(a),fb=func(b),fc,p,q,r,s,tol1,xm;
	//cout << "   x1  " << a  << "   x1  " << b << endl;
	//cout << " f(x1) " << fa << " f(x2) " << fb << endl;
	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
		throw("Root must be bracketed in zbrent");
	fc=fb;
	for (int iter=0;iter<ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (abs(fc) < abs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPS*abs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (abs(xm) <= tol1 || fb == 0.0) return b;
		if (abs(e) >= tol1 && abs(fa) > abs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=abs(p);
			double min1=3.0*xm*q-abs(tol1*q);
			double min2=abs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (abs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1,xm);
			fb=func(b);
	}
	throw("Maximum number of iterations exceeded in zbrent");
	return 0.0;
}
inline double NR::SIGN(const double &a, const double &b)
	{return (double)(b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));}




double dfridr(double func(double), const double x/*, const double h, double &err*/)
{
	double h=1.0e-3, err;

	const int ntab=10;
	const double con=1.4, con2=(con*con);
	const double big=numeric_limits<double>::max();
	const double safe=2.0;
	int i,j;
	double errt,fac,hh,ans=0.0;
	double a[ntab][ntab];
	if (h == 0.0) throw("h must be nonzero in dfridr.");
	hh=h;
	a[0][0]=(func(x+hh)-func(x-hh))/(2.0*hh);
	err=big;
	for (i=1;i<ntab;i++) {
		hh /= con;
		a[0][i]=(func(x+hh)-func(x-hh))/(2.0*hh);
		fac=con2;
		for (j=1;j<=i;j++) {
			a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
			fac=con2*fac;
			if(abs(a[j][i]-a[j-1][i])>abs(a[j][i]-a[j-1][i-1]))
				errt=abs(a[j][i]-a[j-1][i]);
			else
				errt=abs(a[j][i]-a[j-1][i-1]);
			if (errt <= err) {
				err=errt;
				ans=a[j][i];
			}
		}
		if (abs(a[i][i]-a[i-1][i-1]) >= safe*err) break;
	}
	return ans;
}

double NR::dfridrX(double func(double,double), const double x, const double y/*, const double h, double &err*/)
{
	double h=1.0e-3, err;

	const int ntab=10;
	const double con=1.4, con2=(con*con);
	const double big=numeric_limits<double>::max();
	const double safe=2.0;
	int i,j;
	double errt,fac,hh,ans=0.0;
	double a[ntab][ntab];
	if (h == 0.0) throw("h must be nonzero in dfridr.");
	hh=h;
	a[0][0]=(func(x+hh,y)-func(x-hh,y))/(2.0*hh);
	err=big;
	for (i=1;i<ntab;i++) {
		hh /= con;
		a[0][i]=(func(x+hh,y)-func(x-hh,y))/(2.0*hh);
		fac=con2;
		for (j=1;j<=i;j++) {
			a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
			fac=con2*fac;
			if(abs(a[j][i]-a[j-1][i])>abs(a[j][i]-a[j-1][i-1]))
				errt=abs(a[j][i]-a[j-1][i]);
			else
				errt=abs(a[j][i]-a[j-1][i-1]);
			if (errt <= err) {
				err=errt;
				ans=a[j][i];
			}
		}
		if (abs(a[i][i]-a[i-1][i-1]) >= safe*err) break;
	}
	return ans;
}

double NR::dfridrY(double func(double,double), const double x, const double y/*, const double h, double &err*/)
{
	double h=1.0e-3, err;

	const int ntab=10;
	const double con=1.4, con2=(con*con);
	const double big=numeric_limits<double>::max();
	const double safe=2.0;
	int i,j;
	double errt,fac,hh,ans=0.0;
	double a[ntab][ntab];
	if (h == 0.0) throw("h must be nonzero in dfridr.");
	hh=h;
	a[0][0]=(func(x,y+hh)-func(x,y-hh))/(2.0*hh);
	err=big;
	for (i=1;i<ntab;i++) {
		hh /= con;
		a[0][i]=(func(x,y+hh)-func(x,y-hh))/(2.0*hh);
		fac=con2;
		for (j=1;j<=i;j++) {
			a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
			fac=con2*fac;
			if(abs(a[j][i]-a[j-1][i])>abs(a[j][i]-a[j-1][i-1]))
				errt=abs(a[j][i]-a[j-1][i]);
			else
				errt=abs(a[j][i]-a[j-1][i-1]);
			if (errt <= err) {
				err=errt;
				ans=a[j][i];
			}
		}
		if (abs(a[i][i]-a[i-1][i-1]) >= safe*err) break;
	}
	return ans;
}
