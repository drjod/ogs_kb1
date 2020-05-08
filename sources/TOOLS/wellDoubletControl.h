#ifndef WELL_DOUBLET_CONTROL_H
#define WELL_DOUBLET_CONTROL_H

#include <string>

#include "wdc_config.h"
#include "comparison.h"

namespace wdc
{


// const double c_well_shutdown_temperature_range = 10.;

// for schemes 1, 2
// well doublet is shutdown if storage becomes full (when storing) or empty (when extracting) 
// i.e.  flow and powerrate gruadually become zero if
// 	storing:  T_UA (cold well) exceeds value_threshold - c_well_shutdown_temperature_range
//	extracting: T_HE (warm well) falls below value_threshlod + c_well_shutdown_temperature_range 
const double c_powerrate_adaption_factor = 1.;
// take 1. as value
const double c_flowrate_adaption_factor = .5;
// used in schemes 1 2 - it is modified during iteration with a mutable variable
// Q_W = Q_W (1 +/- a (T_1 - value_target) / (T_HE - T_UA)) for scheme 1
// it is multiplied with itself, if a threshold is hit
// therefore: DO NOT USE 1 as value


class WellDoubletControl
{
public:
	enum storage_state_t { powerrate_to_adapt, on_demand, target_not_achievable };
	struct result_t
	{
		double Q_H, Q_W;  // power rate Q_H (given input, potentially adapted),
                                // flow rate Q_W either calculated
                                // (schemes A and C) or set as input (Scheme B)
        	double T_HE, T_UA;  // temperature at heat exchanger and in upwind aquifer
        	storage_state_t storage_state;  // says if storage meets demand or not
	};

	struct balancing_properties_t  // helper struct to pass properties around
	{
		double T_HE, T_UA, volumetricHeatCapacity_HE, volumetricHeatCapacity_UA;
	};
	struct accuracies_t
	{
		double temperature;  // 1.e-1  // for thresholds
		double powerrate;  // 10.  // for simulator
		double flowrate;    // 1.e-5  // used for comparison with threshold and zero, it is also minimum absolute flowrate
	};
private:
	result_t result;  // for the client
	int _scheme_ID;
protected:
	accuracies_t accuracies; // const
	double well_shutdown_temperature_range;  // 10. - to shut down if storage is full or empty 

	double Q_H_old;
	double Q_W_old;  // looks redundant

	WellDoubletControl(int __scheme_ID, double _well_shutdown_temperature_range, accuracies_t _accuracies) : 
		_scheme_ID(__scheme_ID), well_shutdown_temperature_range(_well_shutdown_temperature_range), 
				accuracies(_accuracies), value_target(0.){} 

	void set_flowrate(const double& _Q_W)
	{ 
		result.Q_W = _Q_W; 
		LOG("\t\t\tset flow rate\t" << _Q_W);
	}

	void set_powerrate(const double& _Q_H) 
	{ 
		result.Q_H = _Q_H;
		LOG("\t\t\tadapt power rate\t" << _Q_H);
		result.storage_state = powerrate_to_adapt;
	}

	void set_storage_state(storage_state_t _storage_state)
	{
		result.storage_state = _storage_state;
		if(result.storage_state == powerrate_to_adapt)
			LOG("\t\tpowerrate to adapt");
		if(result.storage_state == target_not_achievable)
			LOG("\t\ttarget not achievable");
	}

	double volumetricHeatCapacity_HE, volumetricHeatCapacity_UA;  // parameter
	double value_target, value_threshold;  // constraints
	// 0: T_HE_target, Q_W_max 
	// 1: Q_W_target, T_HE_max 
	// 2: DT_target, Q_W_max

	enum {storing, extracting} operationType;

	wdc::Comparison beyond, notReached;

	void set_balancing_properties(const balancing_properties_t& balancing_properites);
					// called in evaluate_simulation_result
	virtual void estimate_flowrate() = 0;
	void write_outputFile() const;
public:
	int scheme_ID() const { return _scheme_ID; }

	virtual ~WellDoubletControl() = default;

	result_t get_result() const { return result; }
	virtual void configure_scheme() = 0;
	virtual void evaluate_simulation_result(const balancing_properties_t& balancing_properites) = 0;
	
	void configure(const double& _Q_H,  
		const double& _value_target, const double& _value_threshold, 
		const balancing_properties_t& balancing_properites);
			// constraints are set at beginning of time step

	static WellDoubletControl* create_wellDoubletControl(const int& selection, 
		const double& well_shutdown_temperature_range, const accuracies_t& _accuracies);
			// instance is created before time-stepping
	
	void print_temperatures() const;
	bool powerrate_converged() const { return fabs(result.Q_H - Q_H_old) < accuracies.powerrate; }
	//virtual bool converged(double _T_HE, double accuracy) const = 0;
	virtual bool flowrate_converged() const = 0;
	bool converged() const { return flowrate_converged() && powerrate_converged(); }
	accuracies_t get_accuracies() const { return accuracies; } 
};


class WellScheme_0 : public WellDoubletControl
{
        void estimate_flowrate() override;
        void adapt_powerrate();

public:
	void configure_scheme() override;
	WellScheme_0(const double& _well_shutdown_temperature_range, const accuracies_t& _accuracies) : 
		WellDoubletControl(0, _well_shutdown_temperature_range, _accuracies) {}

	void evaluate_simulation_result(const balancing_properties_t& balancing_properites) override;

	//wdc_state get_state(double, double) const override { return wdc_state; }
		// convergence exclusively decided by simulator in scheme 0 - no iterations in wdc
	bool flowrate_converged() const override { return true; }
};

class WellScheme_1 : public WellDoubletControl
{
	double flowrate_adaption_factor, deltaTsign_stored;
	// to store values from the last interation when adapting flowrate
        void estimate_flowrate() override;
       	void adapt_flowrate();
        void adapt_powerrate();
public:
	void configure_scheme() override;
	WellScheme_1(const double& _well_shutdown_temperature_range, const accuracies_t& _accuracies) : 
		WellDoubletControl(1, _well_shutdown_temperature_range, _accuracies) {}

	void evaluate_simulation_result(const balancing_properties_t& balancing_properites) override;
	//bool converged(double _T_HE, double _accuracy) const override { return (fabs(_T_HE - value_target) < _accuracy || flag_converged); }
	bool flowrate_converged() const override
	{ 
		if(notReached(get_result().T_HE, value_target) && get_result().storage_state == on_demand)
			return false;
		return fabs(get_result().Q_W - Q_W_old) < accuracies.flowrate; 
	}
};


class WellScheme_2 : public WellDoubletControl
{
	double flowrate_adaption_factor, deltaTsign_stored;
	// to store values from the last interation when adapting flowrate
        void estimate_flowrate() override;
       	void adapt_flowrate();
        void adapt_powerrate();
public:
	void configure_scheme() override;
	WellScheme_2(const double& _well_shutdown_temperature_range, const accuracies_t& _accuracies) : 
		WellDoubletControl(2, _well_shutdown_temperature_range, _accuracies) {}

	void evaluate_simulation_result(const balancing_properties_t& balancing_properites) override;
	//bool converged(double _T_HE, double _accuracy) const override { return (fabs(_T_HE - value_target) < _accuracy || flag_converged); }
	bool flowrate_converged() const override 
	{ 
		if(beyond(get_result().T_HE, value_target) && get_result().storage_state == on_demand)
			return false;
		return fabs(get_result().Q_W - Q_W_old) < accuracies.flowrate; 
	}
};


} // end namespace wdc

#endif
