#include "OGS_WDC.h"


// configuration for time step (parameter list, measurement mesh nodes)
// must be called (at least) at beginning of each time step before WDC is created (in CSourceTerm::apply_wellDoubletControl())
long OGS_WDC::get_upwind_aquifer_mesh_node(const double& current_time)
{
	// delete parameter list entry if old such that current entry is at begin of list (list entries have to be ordered in time)
	while(current_time > parameter_list.begin()->time)
		parameter_list.erase(parameter_list.begin());

	if(parameter_list.empty())
		throw std::runtime_error("ERROR in WDC - No Parameter");

	if(parameter_list.begin()->powerrate > 0)
		return get_doublet_mesh_nodes().well2_aquifer;  // storing
	else
		return get_doublet_mesh_nodes().well1_aquifer; // extracting
}

// creates and configures new wellDoubletControl
// must be called in first LIQUID_FLOW call in each time step
// initial flow rate estimation is available after that
// Requirements:
// 1. measurement mesh nodes must be up to date by calling OGS_WDC::update_measurement_mesh_nodes()
// 2. when flow and transport have converged OGS_WDC::discard() must be called
void OGS_WDC::create_new_WDC(const wdc::WellDoubletControl::balancing_properties_t& balancing_properties)
{
	// !!!!! parameterList must have been updated before (by calling update_measurement_mesh_nodes())
	std::cout << "\t\t\tWDC\n";
	wellDoubletControl.reset(wdc::WellDoubletControl::create_wellDoubletControl(parameter_list.begin()->indicator,
			{ accuracy_temperature, accuracy_powerrate, accuracy_flowrate}));
	wellDoubletControl->configure(parameter_list.begin()->powerrate,
						parameter_list.begin()->target_value, parameter_list.begin()->threshold_value,
						balancing_properties);
	is_initialized = true;
}

// called with heat transport (but not in FCT correction)
// heat power rate is available after that
// also flow rates are updated if required
void OGS_WDC::evaluate_simulation_result(const wdc::WellDoubletControl::balancing_properties_t& balancing_properties)
{
	std::cout << "\t\t\tWDC\n";
	wellDoubletControl->evaluate_simulation_result(balancing_properties);
}

// interface function between WDC and CSourceTerm::apply_wellDoubletControl()
// processes WDC calls and returns source term value (flow rate in case of flow and power rate in case of heat)
double OGS_WDC::call_WDC(const CRFProcess* m_pcs, const wdc::WellDoubletControl::balancing_properties_t& balancing_properties)
{
	double result;
	switch(m_pcs->getProcessType())
	{
		case FiniteElement::LIQUID_FLOW:
			// result is flow rate
			if(!is_initialized)// && nodes_counter == 0)
			{
				create_new_WDC(balancing_properties);
				is_evaluated = true;
			}
			else if(!is_evaluated)
			{
				evaluate_simulation_result(balancing_properties);
				is_evaluated = true;
			}

			result = wellDoubletControl->get_result().Q_W;

			// std::cout << "T_HE:\t" << balancing_properties.T_HE << std::endl;
			// std::cout << "T_UA:\t" << balancing_properties.T_UA << std::endl;
			// std::cout << "Q_w: " << wellDoubletControl->get_result().Q_W << std::endl;
			break;
		case FiniteElement::HEAT_TRANSPORT:
			// result is power rate (flow rate might be updated)
			/*if(m_pcs->iter_outer_cpl > 0 && !m_pcs->inFemFCTmode() && nodes_counter == 0)
			{
				evaluate_simulation_result(balancing_properties);
			}*/
			result = wellDoubletControl->get_result().Q_H / heatExchangerArea;
			// std::cout << "T_HE:\t" << balancing_properties.T_HE << std::endl;
			// std::cout << "T_UA:\t" << balancing_properties.T_UA << std::endl;
			// std::cout << "Q_H / A:\t" << wellDoubletControl->get_result().Q_H / heatExchangerArea << std::endl;
			break;
		default:
			throw std::runtime_error("WellDoubletControl - PCS not supported");
	}
	nodes_counter++;
	if(nodes_counter == doublet_mesh_nodes.heatExchanger.size())
		nodes_counter = 0;
	return result;
}


/*
double OGS_WDC::call_WDC(const CRFProcess* m_pcs, const wdc::WellDoubletControl::balancing_properties_t& balancing_properties)
{
	double result;
	switch(m_pcs->getProcessType())
	{
		case FiniteElement::LIQUID_FLOW:
			// result is flow rate
			if(!is_initialized && nodes_counter == 0)
			{
				create_new_WDC(balancing_properties);
			}
			result = wellDoubletControl->get_result().Q_W;

			std::cout << "T_HE:\t" << balancing_properties.T_HE << std::endl;
			std::cout << "T_UA:\t" << balancing_properties.T_UA << std::endl;
			std::cout << "Q_w: " << wellDoubletControl->get_result().Q_W << std::endl;
			break;
		case FiniteElement::HEAT_TRANSPORT:
			// result is power rate (flow rate might be updated)
			if(m_pcs->iter_outer_cpl > 0 && !m_pcs->inFemFCTmode() && nodes_counter == 0)
			{
				evaluate_simulation_result(balancing_properties);
			}
			result = wellDoubletControl->get_result().Q_H / heatExchangerArea;
			std::cout << "T_HE:\t" << balancing_properties.T_HE << std::endl;
			std::cout << "T_UA:\t" << balancing_properties.T_UA << std::endl;
			std::cout << "Q_H / A:\t" << wellDoubletControl->get_result().Q_H / heatExchangerArea << std::endl;
			break;
		default:
			throw std::runtime_error("WellDoubletControl - PCS not supported");
	}
	nodes_counter++;
	if(nodes_counter == doublet_mesh_nodes.heatExchanger.size())
		nodes_counter = 0;
	return result;
}
*/
// for ST on polyline
double OGS_WDC::get_extremum(const CRFProcess* m_pcs, const int& ndx, const std::vector<size_t> nodes) const
{
	/*double result{0.};
	for(auto it = nodes.begin(); it != nodes.end(); ++it)
		result += m_pcs->GetNodeValue(*it, ndx);
	return result / nodes.size();
	*/
	if(parameter_list.begin()->powerrate > 0)
		return m_pcs->GetMaxNodeValue(nodes, ndx);
	else
		return m_pcs->GetMinNodeValue(nodes, ndx);

}

