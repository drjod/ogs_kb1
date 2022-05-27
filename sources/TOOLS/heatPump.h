#ifndef HEAT_PUMP_H
#define HET_PUMP_H

#include "wdc_config.h"

namespace wdc
{

class HeatPump
{
protected:
	double COP;
	double heat_sink;
public:
	HeatPump() : COP(-1.) {}
	double get_COP() const { return COP; }
	double get_heat_sink() const { return heat_sink; }
	virtual double calculate_heat_source(const double& heat_sink,
			const double& T_source_in, const double& T_source_out) = 0;
	virtual double get_heat_sink(const double& heat_source) const = 0;
	virtual double get_parameter() const = 0;
	virtual ~HeatPump() = default;
};


class CarnotHeatPump : public HeatPump
{
	double T_sink;  // into sink
       	double	eta;  // Guetefaktor
public:
	CarnotHeatPump(const double& _T_sink, const double& _eta) : T_sink(_T_sink), eta(_eta) 
	{
		//WDC_LOG("Carnot - T_sink: " << T_sink << ", eta: " << eta);
	}
	double calculate_heat_source(const double& heat_sink, 
			const double& T_source_in, const double& T_source_out) override;
	double get_heat_sink(const double& heat_source) const override { return heat_source * COP / (COP-1); }
	double get_parameter() const override { return eta; }
};


class NoHeatPump : public HeatPump
{
public:
	double calculate_heat_source(const double& heat_sink, 
			const double& T_source_in, const double& T_source_out) override;
	double get_heat_sink(const double& heat_source) const override { return heat_source; }
	double get_parameter() const override { return -1.; }
};

}

#endif

