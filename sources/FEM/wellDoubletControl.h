#ifndef WELLDOUBLETCONTROL_H
#define WELLDOUBLETCONTROL_H

#define ACCURACY_FLOWRATE_TARGET 1.e-4
#define ACCURACY_TEMPERATURE_THRESHOLD 1.e-2
#define FLOWRATE_ADAPTION_FACTOR 0.7
// Q_w = Q_w (1 +/- a (T_1 - value_target)) 

#include "comparison.h"
#include <string>

class WellDoubletControl
{
	std::string name;  // for output file
	int timeStep;
	double simulationTime;
public:
	enum iterationState_t {searchingFlowrate, 
				searchingPowerrate, converged};
	// schemes A/C: first flow rate is adapted, then powerrate
	// scheme B: only powerrate is adapted
	struct result_t
	{
		double Q_H, Q_w;  // power rate Q_H (given input, potentially adapted),
                                // flow rate Q_w either calculated
                                // (schemes A and C) or set as input (Scheme B)
        	double T1, T2;  // temperature at well 1 ("warm well") 
                        // and 2 ("cold well") obtained from simulation results
        	bool flag_powerrateAdapted;  // says if storage meets requirement or not
	};
protected:
	result_t result;
	double heatCapacity1, heatCapacity2;  // Schemes A/B: 1: at warm well, 2: ar cold well
				// Scheme C: 1 for both wells

	char scheme_identifier; // 'A', 'B', or 'C'	
	double value_target, value_threshold;
	// A: T1_target, Q_w_max 
	// B: Q_w_target, T1_max 
	// C: DT_target, Q_w__max

	enum {storing, extracting} operationType;
	iterationState_t iterationState;

	wdc::Comparison beyond, notReached;

	void set_heatFluxes(const double& _T1, const double& _T2,
				const double& _heatCapacity1, const double& _heatCapacity2);  
					// called in evaluate_simulation_result
	virtual void set_flowrate() = 0;
	void write_outputFile() const;
public:
	WellDoubletControl() = default;
	virtual ~WellDoubletControl() = default;

	const WellDoubletControl::result_t& get_result() const { return result; }

	virtual void configureScheme() = 0;
	virtual const iterationState_t& evaluate_simulation_result(
				const double& _T1, const double& _T2,
				const double& _heatCapacity1, const double& _heatCapacity2) = 0;
	
	void configure(const char* _name, int _timeStep, double _simulationTime,
		const double& _Q_H,  
		const double& _value_target, const double& _value_threshold,
		const double& _T1, const double& _T2,
		const double& _heatCapacity1, const double& _heatCapacity2);
			// constraints are set at beginning of time step

	static WellDoubletControl* create_wellDoubletControl(
				const char& selection);
			// instance is created before time-stepping
	
	void print_temperatures() const;
};

class WellSchemeAC : public WellDoubletControl
{
        double temperature_well1() const { return result.T1; }  // to evaluate target
        double temperature_difference_well1well2() const { return result.T1 - result.T2; }
                                                        // to evaluate target
        double (WellSchemeAC::*simulation_result_aiming_at_target) () const;
	// to differentiate between scheme A and C
        // scheme A: pointing at temperature_Well1(),
        // scheme C: pointing at temperature_difference(),
	double factor, deltaTsign_stored;
	// to store values from the last interation when adapting flowrate
        void set_flowrate();
       	void adapt_flowrate();
        void adapt_powerrate();
	void configureScheme();
public:
	WellSchemeAC(const char& _scheme_identifier)
	{ scheme_identifier = _scheme_identifier; }

	const iterationState_t& evaluate_simulation_result(
				const double& _T1, const double& _T2,
				const double& _heatCapacity1, const double& _heatCapacity2);
};


class WellSchemeB : public WellDoubletControl
{
	void configureScheme();
public:
	WellSchemeB(const char& _scheme_identifier)
	{ scheme_identifier = _scheme_identifier; }

	const iterationState_t& evaluate_simulation_result(
				const double& _T1, const double& _T2,
				const double& _heatCapacity1, const double& _heatCapacity2);

        void set_flowrate();
        void adapt_powerrate();
};

#endif
