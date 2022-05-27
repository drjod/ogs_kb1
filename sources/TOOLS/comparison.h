#ifndef WDCUTILITIES_H
#define WDCUTILITIES_H

#include <iostream>
#include <cmath>
#include <cassert>
#include <memory>

namespace wdc
{

struct Greater;
struct Smaller;

struct ComparisonMethod
{
        virtual bool execute(const double& x, const double& y) const = 0;
        virtual ~ComparisonMethod() {}
};

struct Greater : public ComparisonMethod
{
        double epsilon;
        Greater(const double& _epsilon) { epsilon = _epsilon; }
        bool execute(const double& x, const double& y) const 
        { return x > y + epsilon; }
};

struct Smaller : public ComparisonMethod
{
        double epsilon;
        Smaller(const double& _epsilon) { epsilon = _epsilon; }
        bool execute(const double& x, const double& y) const
        { return x < y - epsilon; }
};

struct Comparison
{
        std::unique_ptr<ComparisonMethod> method;

        Comparison() = default;
        Comparison(ComparisonMethod* _method) : method(std::unique_ptr<ComparisonMethod>(_method)) {}

        void configure(ComparisonMethod* _method)
	{
		method = std::unique_ptr<ComparisonMethod>(_method);
	}

        bool operator()(const double& x, const double& y) const
        {
                return method->execute(x, y);
        }

};

inline double make_confined(const double& value, const double& lower_limit, const double& upper_limit)
{
	return fmin(fmax(value, lower_limit), upper_limit);
}


inline int sign(const double& x)
{
	return (x>0.) ? 1 : ((x<0.) ? -1 : 0);
}


enum threshold_t{lower, upper};
// function that returns floating point number between [0.,1.] if a parameter (named value) is close to a threshold (upper or lower)
// returns:
//     0.0: if value is at threshold or beyond
//     1.0: if distance to threshold is larger than delta (but value has not reached threshold)
//     1.0 reduced by S-curve (for smoothing): if distance of value to the threshold is smaller than delta
inline double make_threshold_factor(const double& value, const double& threshold_value, const double& delta, const threshold_t& threshold)
{
	//assert(delta <= 0);
	if(delta <=0) return 1.;

	const int sigma = (threshold == lower)? 1: -1;
	const double U = sigma * (value - threshold_value) / delta;

	return (U <= 0.) ? 0. : (( U < 1. ) ? pow(U, 2*(1-U)) : 1. );
}


}  // end namespace wdc


#endif
