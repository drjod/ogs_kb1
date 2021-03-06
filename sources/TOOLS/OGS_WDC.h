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
	{	// geometric locations
		long well1_aquifer;
		long well2_aquifer;
		std::vector<size_t> heatExchanger;
	};
private:
	std::shared_ptr<wdc::WellDoubletControl> wellDoubletControl;
	std::list<parameter_group_t> parameter_list;

	bool is_initialized; // new WDC for each new time step - WDC creation is controlled by this flag
	bool is_evaluated;  // to do wdc evaluation only once each iteration beginning with second iteration
	doublet_mesh_nodes_t doublet_mesh_nodes;
	size_t nodes_counter;
	double heatExchangerArea;
	double well_shutdown_temperature_range;
	double accuracy_temperature, accuracy_powerrate, accuracy_flowrate;

	void create_new_WDC(const wdc::WellDoubletControl::balancing_properties_t& balancing_properties);
	void evaluate_simulation_result(const wdc::WellDoubletControl::balancing_properties_t& balancing_properties);
public:
	OGS_WDC(const double& _well_shutdown_temperature_range, const double& _accuracy_temperature,
			const double& _accuracy_powerrate, const double& _accuracy_flowrate) :
		is_initialized(false), is_evaluated(true),
		doublet_mesh_nodes({-1, -1, std::vector<size_t>()}),
		nodes_counter(0), heatExchangerArea(1.),
		well_shutdown_temperature_range(_well_shutdown_temperature_range),
		accuracy_temperature(_accuracy_temperature), accuracy_powerrate(_accuracy_powerrate), accuracy_flowrate(_accuracy_flowrate) {}
	std::shared_ptr<wdc::WellDoubletControl> get_WellDoubletControl() const { return wellDoubletControl; }
	void discard() { is_initialized = false; }  // must be called, when iteration loop between LIQUID and HEAT has converged(than new WDC will be created in next time step)
	void set_unevaluated() { is_evaluated = false; }  // call this at end of each iteration - than wdc will be evaluated once in next time step
	template<typename... Args>
	void add_parameterGroup(Args&&...args) { parameter_list.emplace_back(std::forward<Args>(args)...); }
	doublet_mesh_nodes_t get_doublet_mesh_nodes() const { return doublet_mesh_nodes; }
	void set_doublet_mesh_nodes(doublet_mesh_nodes_t _doublet_mesh_nodes) { doublet_mesh_nodes = _doublet_mesh_nodes; }
																		// called in set functions for source terms
	long get_upwind_aquifer_mesh_node(const double& current_time); // for time step
	double call_WDC(const CRFProcess* m_pcs, const wdc::WellDoubletControl::balancing_properties_t& balancing_properties);

	double get_extremum(const CRFProcess* m_pcs, const int& ndx, const std::vector<size_t> nodes) const;

	void set_heatExchangerArea(double _heatExchangerArea) { heatExchangerArea = _heatExchangerArea; }
};


#endif
