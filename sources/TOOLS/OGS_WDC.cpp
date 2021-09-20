#include "OGS_WDC.h"
#include "rf_pcs.h"

#include <fstream>
#include <sstream>

#if defined(USE_MPI)
#include <mpi.h>
#endif


double interpolate( vector<double> &xData, vector<double> &yData, double x, bool extrapolate )
{
   int size = xData.size();

   int i = 0;                                                                  // find left end of interval for interpolation
   if ( x >= xData[size - 2] )                                                 // special case: beyond right end
   {
      i = size - 2;
   }
   else
   {
      while ( x > xData[i+1] ) i++;
   }
   double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1];      // points on either side (unless beyond ends)
   if ( !extrapolate )                                                         // if beyond ends of array and not extrapolating
   {
      if ( x < xL ) yR = yL;
      if ( x > xR ) yL = yR;
   }

   double dydx = ( yR - yL ) / ( xR - xL );                                    // gradient

   return yL + dydx * ( x - xL );                                              // linear interpolation
}


// configuration for time step (parameter list, measurement mesh nodes)
// must be called (at least) at beginning of each time step before WDC is created (in CSourceTerm::apply_wellDoubletControl())
// return value 1: stroing, -1: retrieving
int OGS_WDC::get_aquifer_mesh_nodes(const double& current_time,
		const bool & wdc_flag_extract_and_reinject,
		std::vector<size_t>& heatExchanger_aquifer_mesh_nodes,
		std::vector<double>& heatExchanger_aquifer_mesh_nodes_area_fraction,
		std::vector<size_t>& upwind_aquifer_mesh_nodes,
		std::vector<double>& upwind_aquifer_mesh_nodes_area_fraction)
{
	// delete parameter list entry if old such that current entry is at begin of list (list entries have to be ordered in time)
	while(parameter_list.size()>0 && current_time > parameter_list.begin()->time)
		parameter_list.erase(parameter_list.begin());

	if(parameter_list.empty())
		throw std::runtime_error("ERROR in WDC - Parameter list empty");

	if(wdc_flag_extract_and_reinject)
	{	// heat exchanger at warm well when storing and at cold well when retrieving

		if(parameter_list.begin()->powerrate > 0)
		{	// storing
			heatExchanger_aquifer_mesh_nodes = get_doublet_mesh_nodes().well1_aquifer;
			heatExchanger_aquifer_mesh_nodes_area_fraction =  get_doublet_mesh_nodes().well1_aquifer_area_fraction;

			upwind_aquifer_mesh_nodes = get_doublet_mesh_nodes().well2_aquifer;
			upwind_aquifer_mesh_nodes_area_fraction =  get_doublet_mesh_nodes().well2_aquifer_area_fraction;
			return 1;
		}
		else
		{	// retrieving
			heatExchanger_aquifer_mesh_nodes = get_doublet_mesh_nodes().well2_aquifer;
			heatExchanger_aquifer_mesh_nodes_area_fraction =  get_doublet_mesh_nodes().well2_aquifer_area_fraction;

			upwind_aquifer_mesh_nodes = get_doublet_mesh_nodes().well1_aquifer;
			upwind_aquifer_mesh_nodes_area_fraction = get_doublet_mesh_nodes().well1_aquifer_area_fraction;
			return -1;
		}
	}
	else
	{	// heat exchanger at fixed position given by $GEO_TYPE
		heatExchanger_aquifer_mesh_nodes = get_doublet_mesh_nodes().heatExchanger;
		heatExchanger_aquifer_mesh_nodes_area_fraction =  get_doublet_mesh_nodes().heatExchanger_area_fraction;
		if(parameter_list.begin()->powerrate > 0)
		{  // storing
			upwind_aquifer_mesh_nodes = get_doublet_mesh_nodes().well2_aquifer;
			upwind_aquifer_mesh_nodes_area_fraction =  get_doublet_mesh_nodes().well2_aquifer_area_fraction;
			return 1;
		}
		else
		{  // retrieving
			upwind_aquifer_mesh_nodes = get_doublet_mesh_nodes().well1_aquifer;
			upwind_aquifer_mesh_nodes_area_fraction =  get_doublet_mesh_nodes().well1_aquifer_area_fraction;
			return -1;
		}
	}

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
	std::cout << "\tWDC create"
#if defined(USE_MPI)
				<< " - pcs: " << myrank
#endif
	<< '\n';

	wellDoubletControl.reset(wdc::WellDoubletControl::create_wellDoubletControl(parameter_list.begin()->indicator,
			well_shutdown_temperature_range,
			{ accuracy_temperature, accuracy_powerrate, accuracy_flowrate}));

	double target_value;
	if(!heat_pump_flag || parameter_list.begin()->powerrate > 0.)
		target_value = parameter_list.begin()->target_value;
	else
	{
		const double T_UA = (balancing_properties.T_UA > 200.)?  balancing_properties.T_UA - 273.15 : balancing_properties.T_UA;

		target_value = (parameter_list.begin()->indicator == 2)?
				-interpolate(heat_pump_parameter.T_source, heat_pump_parameter.Delta_T, T_UA, true) :
				parameter_list.begin()->target_value;

		const double eta = interpolate(heat_pump_parameter.T_source, heat_pump_parameter.eta, T_UA, true);
		wellDoubletControl->set_heatPump(heat_pump_flag, parameter_list.begin()->temperature_sink, eta);
	}

	wellDoubletControl->configure(parameter_list.begin()->powerrate,
						target_value, parameter_list.begin()->threshold_value,
						balancing_properties);
	temperature_sink = parameter_list.begin()->temperature_sink;
	is_initialized = true;
}

// called with heat transport (but not in FCT correction)
// heat power rate is available after that
// also flow rates are updated if required
void OGS_WDC::evaluate_simulation_result(const wdc::WellDoubletControl::balancing_properties_t& balancing_properties)
{
	std::cout << "\tWDC - evaluate"
#if defined(USE_MPI)
				<< " - pcs: " << myrank
#endif
	<< '\n';

	double target_value = parameter_list.begin()->target_value;

	if(heat_pump_flag & parameter_list.begin()->powerrate < 0.)
	{
		const double T_UA = (balancing_properties.T_UA > 200.)?  balancing_properties.T_UA - 273.15 : balancing_properties.T_UA;
		const double eta = interpolate(heat_pump_parameter.T_source, heat_pump_parameter.eta, T_UA, true);

		wellDoubletControl->set_heatPump(heat_pump_flag, parameter_list.begin()->temperature_sink, eta);

		if(parameter_list.begin()->indicator == 2)
			target_value = -interpolate(heat_pump_parameter.T_source, heat_pump_parameter.Delta_T, T_UA, true);

	}

	wellDoubletControl->evaluate_simulation_result(balancing_properties);
}

// interface function between WDC and CSourceTerm::apply_wellDoubletControl()
// processes WDC calls and returns source term value (flow rate in case of flow and power rate in case of heat)
double OGS_WDC::call_WDC(CRFProcess* m_pcs,
		const wdc::WellDoubletControl::balancing_properties_t& balancing_properties,
		std::vector<size_t> heatExchanger_aquifer_mesh_nodes)
{
	double result;
	switch(m_pcs->getProcessType())
	{
		case FiniteElement::LIQUID_FLOW:
			// result is flow rate
			if(parameter_list.begin()->indicator == 0 ||
					parameter_list.begin()->indicator == 1 ||
					parameter_list.begin()->indicator == 2) // via STs
			{
				wdc_result.scheme_ID = parameter_list.begin()->indicator;
				{
					if(!is_initialized)// && nodes_counter == 0)
					{
						create_new_WDC(balancing_properties);
						is_evaluated = true;
					}

					if(!is_evaluated)
					{
						evaluate_simulation_result(balancing_properties);
						is_evaluated = true;
					}

					result = wellDoubletControl->get_result().Q_W;
				}
			}
			else if(parameter_list.begin()->indicator == 3) // via ST and BC
			{	// set ST
				// result = flow rate

				if(fabs(parameter_list.begin()->powerrate) > 1.e-10)
				{
					result = ( fabs(balancing_properties.T_HE - balancing_properties.T_UA) > 1.e-10) ? 
						fabs(parameter_list.begin()->powerrate) / ( (balancing_properties.T_HE - balancing_properties.T_UA) * 
							balancing_properties.volumetricHeatCapacity_HE ) : 
						parameter_list.begin()->threshold_value;

					if (fabs(result) > fabs(parameter_list.begin()->threshold_value) - 1.e-10)
					{
						result = parameter_list.begin()->threshold_value;
						wdc_result.power_rate  = fabs(result) * 
							(balancing_properties.T_HE - balancing_properties.T_UA) *  balancing_properties.volumetricHeatCapacity_HE;
					}
					else
						wdc_result.power_rate = parameter_list.begin()->powerrate;
				}	

				shut_in = false;
				if(fabs(parameter_list.begin()->powerrate) <= 1.e-10 
						|| parameter_list.begin()->powerrate * (parameter_list.begin()->target_value - balancing_properties.T_UA) < 0 
						|| fabs(parameter_list.begin()->target_value - balancing_properties.T_UA) < well_shutdown_temperature_range)
				{
					result = 0.;
					wdc_result.power_rate = 0.;
					shut_in = true;
				}
	

				wdc_result.scheme_ID = 3;
				wdc_result.flow_rate = result;
				wdc_result.power_rate_target = parameter_list.begin()->powerrate;
				wdc_result.T_HE = balancing_properties.T_HE;
				wdc_result.T_UA = balancing_properties.T_UA;

				//if(!is_evaluated)
				//std::cout << "\t\tWDC - Powerrate " << wdc_result.power_rate <<  " - Set flow rate: " << result << '\n';
				is_evaluated = true;
			}
			else
				std::runtime_error("WDC indicator unknown");

			break;
		case FiniteElement::HEAT_TRANSPORT:
			is_evaluated = false;
			if(parameter_list.size() > 0)
			{
				if(parameter_list.begin()->indicator == 0 ||
						parameter_list.begin()->indicator == 1 ||
						parameter_list.begin()->indicator == 2) // via STs
				{
					result = wellDoubletControl->get_result().Q_H; // / heatExchangerArea;
				}
				else if(parameter_list.begin()->indicator == 3 && !shut_in) // via ST and BC
				{ 	// set BC
					for(int i=0; i < heatExchanger_aquifer_mesh_nodes.size(); i++)
						m_pcs->set_BCNode(heatExchanger_aquifer_mesh_nodes[i], parameter_list.begin()->target_value);
				}
				else
					std::runtime_error("WDC indicator unknown");
			}

			break;
		default:
			throw std::runtime_error("WellDoubletControl - PCS not supported");
	}	
	nodes_counter++;
	if(nodes_counter == doublet_mesh_nodes.heatExchanger.size())
		nodes_counter = 0;
			
	return result;
}


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

// for output if parallel
void OGS_WDC::write_logfileheader(const std::size_t& ndx)
{
	if(logging)
	{
		std::stringstream ss;
		ss << "WDC_logfile_" << ndx <<
#if defined(USE_MPI)
				"_" << myrank <<
#endif
				".txt";

		std::fstream fout(ss.str(), std::fstream::out);
 		fout << "scheme_ID = 0, 1, 2:\n\ttime\tscheme_ID\tpower_adaption_flag\tpower_rate\tsystem_power_rate\tsystem_target_power_rate\tflow_rate\ttemperature_warm_well\ttemperature_cold_well\ttemperature_heat_exchanger\ttemperature_sink\tCOP\teta\n";
		fout << "scheme_ID = 3:\n\ttime\tscheme_ID\tpower_rate\tpower_rate_target\tflow_rate\tT_HE (heat exchanger)\tT_UA (upwind aquifer)\n";

	}
}

// for output if parallel
void OGS_WDC::write_logfile(const double& time, const std::size_t& ndx, const CRFProcess* m_pcs)
{
	if(logging)
	{
		std::stringstream ss;
		ss << "WDC_logfile_" << ndx <<
#if defined(USE_MPI)
				"_" << myrank <<
#endif
				".txt";
		std::fstream fout(ss.str(), std::fstream::out | std::fstream:: app);
		if(wdc_result.scheme_ID == 0 || wdc_result.scheme_ID == 1 || wdc_result.scheme_ID == 2) // ST/ST
		{
		    if(wellDoubletControl)
		    {
                       const wdc::WellDoubletControl::result_t& result = wellDoubletControl->get_result();

                       const double system_powerrate = (result.Q_H>0.)?
                               result.Q_H: wellDoubletControl->get_system_powerrate();
                       const double COP = (result.Q_H>0)? -1: wellDoubletControl->get_COP();
                       const double eta = (result.Q_H>0)? -1: wellDoubletControl->get_heatPumpParameter();

                       fout << time
                               << '\t' << wellDoubletControl->get_scheme_ID()
                               << '\t' << result.storage_state   // 0: powerrate_to_adapt, 1: on_demand
                               << '\t' << result.Q_H
                               << '\t' << system_powerrate
                               << '\t' << wellDoubletControl->get_system_target_powerrate()
                               << '\t' << result.Q_W
                               << '\t' << m_pcs->GetWeightedAverageNodeValue(doublet_mesh_nodes.well1_aquifer,
                                               doublet_mesh_nodes.well1_aquifer_area_fraction, 1)
			       << '\t' << m_pcs->GetWeightedAverageNodeValue(doublet_mesh_nodes.well2_aquifer,
                                               doublet_mesh_nodes.well2_aquifer_area_fraction, 1)
                               << '\t' << m_pcs->GetWeightedAverageNodeValue(doublet_mesh_nodes.heatExchanger,
                                               doublet_mesh_nodes.heatExchanger_area_fraction, 1)
                               //m_pcs->ogs_WDC_vector[i]->get_extremum(m_pcs, 1, doublet_mesh_nodes.heatExchanger)
			       << '\t' << get_temperature_sink()
                               << '\t' << COP
							   << '\t' << eta
                               << '\n';
		  }		
		}
		else if(wdc_result.scheme_ID == 3) // BC/ST
		{
			fout << time << '\t'
				<< wdc_result.scheme_ID
				<< '\t' << wdc_result.power_rate
				<< '\t' << wdc_result.power_rate_target
				<< '\t' << wdc_result.flow_rate
				<< '\t' << wdc_result.T_HE
				<< '\t' << wdc_result.T_UA
				<< '\n';
		}
//		else
//			throw std::runtime_error("WDC scheme unknown");

	}
}


void OGS_WDC::set_heat_pump_parameter(const bool & _heat_pump_flag, const std::string& heat_pump_file_name)
{
	heat_pump_flag = _heat_pump_flag;

	if(heat_pump_flag)
	{
		std::string line;
		std::fstream infile(heat_pump_file_name);
		if(infile.is_open())
		{
			std::getline(infile, line); // header

			while (std::getline(infile, line))
			{
				std::istringstream iss(line);
				double T_source, eta, Delta_T;
				iss >> T_source >> eta >> Delta_T;

				heat_pump_parameter.T_source.push_back(T_source);
				heat_pump_parameter.eta.push_back(eta);
				heat_pump_parameter.Delta_T.push_back(Delta_T);
			}
		}
		else
			std::runtime_error("Failed to open heat pump file");
	}
}
