#include "OGS_WDC.h"


// configuration for time step (parameter list, measurement mesh nodes)
// must be called (at least) at beginning of each time step before WDC is created (in CSourceTerm::apply_wellDoubletControl())
OGS_WDC::measurement_mesh_nodes_t OGS_WDC::update_measurement_mesh_nodes(const double& current_time)
{
	// delete parameter list entry if old such that current entry is at begin of list (list entries have to be ordered in time)
	while(current_time > parameter_list.begin()->time)
		parameter_list.erase(parameter_list.begin());

	return OGS_WDC::measurement_mesh_nodes_t {
			// warm well1
			(parameter_list.begin()->powerrate > 0) ?
						get_doublet_mesh_nodes().heatExchanger : // storing
						get_doublet_mesh_nodes().well1_aquifer, // extracting
			// cold well 2
			(parameter_list.begin()->powerrate > 0) ?
						get_doublet_mesh_nodes().well2_aquifer : // storing
						get_doublet_mesh_nodes().heatExchanger // extracting
			};
}

// creates and configures new wellDoubletControl
// must be called in first LIQUID_FLOW call in each time step
// initial flow rate estimation is available after that
// Requirements:
// 1. measurement mesh nodes must be up to date by calling OGS_WDC::update_measurement_mesh_nodes()
// 2. when flow and transport have converged OGS_WDC::discard() must be called
void OGS_WDC::update_WellDoubletControl(const WellDoubletControl::balancing_properties_t& balancing_properties)
{
	// !!!!! parameterList must have been updated before (by calling update_measurement_mesh_nodes())
	std::cout << "\t\t\tWDC\n";
	wellDoubletControl.reset(WellDoubletControl::create_wellDoubletControl(parameter_list.begin()->indicator));
	wellDoubletControl->configure(parameter_list.begin()->powerrate,
						parameter_list.begin()->target_value, parameter_list.begin()->threshold_value,
						balancing_properties);
	is_initialized = true;
}

// called with heat transport (but not in FCT correction)
// heat power rate is available after that
// also flow rates are updated if required
void OGS_WDC::evaluate_simulation_result(const WellDoubletControl::balancing_properties_t& balancing_properties)
{
	std::cout << "\t\t\tWDC\n";
	wellDoubletControl->evaluate_simulation_result(balancing_properties);
}

// interface function between WDC and CSourceTerm::apply_wellDoubletControl()
// processes WDC calls
double OGS_WDC::get_result(const CRFProcess* m_pcs, const WellDoubletControl::balancing_properties_t& balancing_properties)
{
	switch(m_pcs->getProcessType())
	{
		case FiniteElement::LIQUID_FLOW:
			// result is flow rate
			if(!is_initialized)
			{
				update_WellDoubletControl(balancing_properties);
			}
			return wellDoubletControl->get_result().Q_w;

		case FiniteElement::HEAT_TRANSPORT:
			// result is power rate (flow rate might be updated)
			if(m_pcs->iter_outer_cpl > 0 && !m_pcs->inFemFCTmode())
			{
				evaluate_simulation_result(balancing_properties);
			}
			return wellDoubletControl->get_result().Q_H;

		default:
			throw std::runtime_error("WellDoubletControl - PCS not supported");
	}
	return 0.;
}

