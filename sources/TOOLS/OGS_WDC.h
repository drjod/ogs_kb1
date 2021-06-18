// interface class and wrapper for WDC
// including:
//   parameter lists:  time, indicator, powerrate, targetValue, thresholdValue
//   doublet geometry mesh nodes: warm well1, cold well2, heat exchanger
//
// provides mesh nodes where temperature (and volumetric heat capacity) are measured to balance heat input and output
//
// CSourceTerm::apply_wellDoubletControl

#ifndef OGS_WDC_H
#define OGS_WDC_H


class OGS_WDC;
#include "rf_pcs.h"
#include "rf_node.h"

#include <list>
#include <vector>
#include <utility>
#include <iostream>
#include <stdexcept>
#include <memory>
#include "wellDoubletControl.h"


class OGS_WDC
{
public:
	struct parameter_group_t
	{
		double time;
		int indicator;
		double powerrate;
		double target_value;
		double threshold_value;
		parameter_group_t(const double& _time, const int& _indicator,
				const double& _powerrate, const double& _target_value, const double& _threshold_value) :
			time(_time), indicator(_indicator), powerrate(_powerrate),
			target_value(_target_value), threshold_value(_threshold_value) {}
		// constructor used with emplace_back below - do not reorder without considering this
	};

	struct doublet_mesh_nodes_t
	{	// geometry
		std::vector<size_t> well1_aquifer;
		std::vector<size_t> well2_aquifer;
		std::vector<size_t> heatExchanger;
		std::vector<double> well1_aquifer_area_fraction;  // sums up to 1
		std::vector<double> well2_aquifer_area_fraction;  // sums up to 1
		std::vector<double> heatExchanger_area_fraction;  // sums up to 1
	};

	struct heat_pump_parameter_t
	{
		int _type;
		double T_sink;
		double eta;
	};

	struct result_t  // for scheme 3
	{
		int scheme_ID;
		double flow_rate;
		double power_rate;
		double power_rate_target;
		double T_HE;
		double T_UA;
	};

	void write_logfileheader(const std::size_t&);
	void write_logfile(const double&, const std::size_t&, const CRFProcess*);
private:
	std::shared_ptr<wdc::WellDoubletControl> wellDoubletControl;
	std::list<parameter_group_t> parameter_list;
	heat_pump_parameter_t heat_pump_parameter;

	bool is_initialized; // new WDC for each new time step - WDC creation is controlled by this flag
	bool is_evaluated;  // to do wdc evaluation only once each iteration beginning with second iteration
	bool logging;
	doublet_mesh_nodes_t doublet_mesh_nodes;
	size_t nodes_counter;
	//double heatExchangerArea;
	double well_shutdown_temperature_range;
	double accuracy_temperature, accuracy_powerrate, accuracy_flowrate;

	void create_new_WDC(const wdc::WellDoubletControl::balancing_properties_t& balancing_properties);
	void evaluate_simulation_result(const wdc::WellDoubletControl::balancing_properties_t& balancing_properties);

	result_t wdc_result;  // for scheme 3
	bool shut_in;
public:
	double get_well_shutdown_temperature_range() { return well_shutdown_temperature_range; }
	//void set_result(const double& _result) { result = _result; }   // for scheme 3
	result_t get_result() { return wdc_result; }
	std::list<parameter_group_t> get_parameter_list() { return parameter_list; }

	OGS_WDC(const double& _well_shutdown_temperature_range, const double& _accuracy_temperature,
			const double& _accuracy_powerrate, const double& _accuracy_flowrate,
			std::size_t ndx, const bool& logging) :
		is_initialized(false), is_evaluated(true), logging(logging),
		doublet_mesh_nodes({std::vector<size_t>(), std::vector<size_t>(),
					std::vector<size_t>(), std::vector<double>(), std::vector<double>(), std::vector<double>()}),
		nodes_counter(0), // heatExchangerArea(1.),
		well_shutdown_temperature_range(_well_shutdown_temperature_range),
		accuracy_temperature(_accuracy_temperature), accuracy_powerrate(_accuracy_powerrate), accuracy_flowrate(_accuracy_flowrate)
	{ write_logfileheader(ndx); }

	std::shared_ptr<wdc::WellDoubletControl> get_WellDoubletControl() const { return wellDoubletControl; }

	void discard(const double& time, const std::size_t ndx, const CRFProcess* m_pcs)  // must be called, when iteration loop between LIQUID and HEAT has converged(than new WDC will be created in next time step)
		{ is_initialized = false; write_logfile(time, ndx, m_pcs); }
	void set_unevaluated() { is_evaluated = false; }  // call this at end of each iteration - than wdc will be evaluated once in next time step
	template<typename... Args>
	void add_parameterGroup(Args&&...args) { parameter_list.emplace_back(std::forward<Args>(args)...); }
	doublet_mesh_nodes_t get_doublet_mesh_nodes() const { return doublet_mesh_nodes; }
	void set_doublet_mesh_nodes(doublet_mesh_nodes_t _doublet_mesh_nodes) { doublet_mesh_nodes = _doublet_mesh_nodes; }
																		// called in set functions for source terms
	void set_heat_pump_parameter(const int& _type, const double& T_sink, const double& eta)
		{ heat_pump_parameter._type = _type ; heat_pump_parameter.T_sink = T_sink; heat_pump_parameter.eta = eta ;}
	int get_aquifer_mesh_nodes(const double& current_time, // for time step
			const bool & wdc_flag_extract_and_reinject,
			std::vector<size_t>& heatExchanger_aquifer_mesh_nodes,
			std::vector<double>& heatExchanger_aquifer_mesh_nodes_area_fraction,
			std::vector<size_t>& upwind_aquifer_mesh_nodes,
			std::vector<double>& upwind_aquifer_mesh_nodes_area_fraction);
	void set_heat_exchanger_mesh_nodes(std::vector<size_t> mesh_nodes, std::vector<double> mesh_nodes_area_fraction)
			// just for output with WDC at polylines (heat exchanger switches)
			{ doublet_mesh_nodes.heatExchanger = mesh_nodes;
				doublet_mesh_nodes.heatExchanger_area_fraction = mesh_nodes_area_fraction; }
	double call_WDC(CRFProcess* m_pcs,
			const wdc::WellDoubletControl::balancing_properties_t& balancing_properties,
			std::vector<size_t> heatExchanger_aquifer_mesh_nodes);

	double get_extremum(const CRFProcess* m_pcs, const int& ndx, const std::vector<size_t> nodes) const;

	//void set_heatExchangerArea(double _heatExchangerArea) { heatExchangerArea = _heatExchangerArea; }
};


#endif
