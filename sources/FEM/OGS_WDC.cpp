#include "OGS_WDC.h"


// configuration for time step (parameter list, measurement mesh nodes)
// must be called (at least) at beginning of each time step before WDC is created (in CSourceTerm::apply_wellDoubletControl())
OGS_WDC::measurement_mesh_nodes_t OGS_WDC::update_measurement_mesh_nodes(const double& current_time)
{
	// delete parameter list entry if old such that current entry is at begin of list (list entries have to be ordered in time)
	while(current_time > wellDoubletData.parameter_list.begin()->time)
		wellDoubletData.parameter_list.erase(
				wellDoubletData.parameter_list.begin());

	return OGS_WDC::measurement_mesh_nodes_t {
			// warm well1
			(wellDoubletData.parameter_list.begin()->powerrate > 0) ?
						get_doublet_mesh_nodes().heatExchanger : // storing
						get_doublet_mesh_nodes().well1_aquifer, // extracting
			// cold well 2
			(wellDoubletData.parameter_list.begin()->powerrate > 0) ?
						get_doublet_mesh_nodes().well2_aquifer : // storing
						get_doublet_mesh_nodes().heatExchanger // extracting
			};
}

// creates and configures new wellDoubletControl
// must be called in first LIQUID_FLOW call in each time step
// initial flow rate estimation is available after that
// Requirements:
// 1. measurement mesh nodes must be up to date by calling OGS_WDC::update_measurement_mesh_nodes(const double& current_time)
// 2. when flow and transport have converged OGS_WDC::discard() must be called
void OGS_WDC::update_WellDoubletControl(const double& temperature_well1, const double& temperature_well2,
									const double& volumetricHeatCapacity_well1, const double& volumetricHeatCapacity_well2)
{
	// !!!!! parameterList must have been updated before (by calling update_measurement_mesh_nodes())
	std::cout << "\t\t\tWDC\n";
	wellDoubletControl.reset(WellDoubletControl::create_wellDoubletControl(wellDoubletData.parameter_list.begin()->indicator));
	wellDoubletControl->configure(wellDoubletData.parameter_list.begin()->powerrate,
									wellDoubletData.parameter_list.begin()->target_value,
									wellDoubletData.parameter_list.begin()->threshold_value,
									temperature_well1, temperature_well2, volumetricHeatCapacity_well1, volumetricHeatCapacity_well2);
	is_initialized = true;
}


// called with heat transport (but not in FCT correction)
// heat power rate is available afer that
// also flow rates are updated if required
void OGS_WDC::evaluate_simulation_result(const double& temperature_well1, const double& temperature_well2,
									const double& volumetricHeatCapacity_well1, const double& volumetricHeatCapacity_well2)
{
	std::cout << "\t\t\tWDC\n";
	wellDoubletControl->evaluate_simulation_result(temperature_well1, temperature_well2, volumetricHeatCapacity_well1, volumetricHeatCapacity_well2);
}



// interface function to CSourceTerm::apply_wellDoubletControl()
// processes WDC calls
double OGS_WDC::get_result(const CRFProcess* m_pcs, const double& temperature_well1, const double& temperature_well2,
		const double& volumetricHeatCapacity_well1, const double& volumetricHeatCapacity_well2)
{
	switch(m_pcs->getProcessType())
	{
		case FiniteElement::LIQUID_FLOW:
			// result is flow rate
			if(!is_initialized)
			{
				std::cout << "\t\t\tWDC\n";
				update_WellDoubletControl(temperature_well1, temperature_well2, volumetricHeatCapacity_well1, volumetricHeatCapacity_well2);
			}
			return wellDoubletControl->get_result().Q_w;

		case FiniteElement::HEAT_TRANSPORT:
			// result is power rate (flow rate updated if required)
			if(m_pcs->iter_outer_cpl > 0 && !m_pcs->inFemFCTmode())
			{
				std::cout << "\t\t\tWDC\n";
				evaluate_simulation_result(temperature_well1, temperature_well2, volumetricHeatCapacity_well1, volumetricHeatCapacity_well2);
			}
			return wellDoubletControl->get_result().Q_H;

		default:
			throw std::runtime_error("WellDoubletControl - PCS not supported");
	}
	return 0.;
}

