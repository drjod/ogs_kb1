#include "wellDoubletControl.h"
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <utility>


#define LOG(x) std::cout << x << std::endl


void WellDoubletControl::print_temperatures() const
{
	std::cout << "\t\t\tT1: " << result.T1 << " - T2: " <<
		result.T2 << std::endl;
}


void WellDoubletControl::write_outputFile() const
{
	std::fstream stream;
	std::string smiley;

	if(timeStep == 0)
	{
		stream.open(name + std::string(".txt"), std::ios::out);
		stream << "Time step; Simulation time t; scheme; " <<
		"Smiley as power rate adaption identifier; Power rate Q_H; " << 
		"Flow rate Q_w; T_1 at warm well 1; T_2 at cold well" << std::endl;
	}
	else
	{
		stream.open(name + std::string(".txt"), std::ios::out | std::ios::app);
	}

	if(result.flag_powerrateAdapted)
		smiley = ":-(";
	else
		smiley = ":-)";

	stream << timeStep << "\t" << simulationTime << "\t" << scheme_identifier <<
		"\t" << smiley << "\t" << result.Q_H << "\t" << result.Q_w <<
		"\t" << result.T1 << "\t" << result.T2 << std::endl; 

	stream.close();	
}

void WellDoubletControl::configure(
	const char* _name, int _timeStep, double _simulationTime,
	const double& _Q_H,
	const double& _value_target, const double& _value_threshold,
	const double& _T1, const double& _T2, 
	const double& _heatCapacity1, const double& _heatCapacity2) 
{
	name = _name;
	timeStep = _timeStep;
	simulationTime = _simulationTime;

	set_heatFluxes(_T1, _T2, _heatCapacity1, _heatCapacity2);
	// set input values for well doublet control, e.g. from file
	result.Q_H = _Q_H;  // stored (Q_H>0) or extracted (Q_H<0) heat
	result.flag_powerrateAdapted = false;

	value_target = _value_target;
	value_threshold = _value_threshold;
	
	if(_Q_H > 0.)
	{
		LOG("\t\t\tset power rate\t\t" <<_Q_H << " - storing");
		operationType = storing;
	}
	else
	{
		LOG("\t\t\tset power rate\t\t" << _Q_H << " - extracting");
		operationType = extracting;
	}

	// the scheme-dependent stuff
	configureScheme();  // iterationState & comparison functions 
			//for temperature target (A, C), temperature constraint (B)
	set_flowrate();  // an estimation for scheme A and a target for scheme B
}

void WellDoubletControl::set_heatFluxes(const double& _T1, const double& _T2, 
			const double& _heatCapacity1, const double& _heatCapacity2) 
{
	result.T1 = _T1;  // warm well
	result.T2 = _T2;  // cold well
	heatCapacity1 = _heatCapacity1;  // warm well
	heatCapacity2 = _heatCapacity2;  // cold well
	LOG("\t\t\tset temperatures\twarm: " << _T1 << "\t\tcold: " << _T2);
	LOG("\t\t\tset heat capacities\twarm: " << _heatCapacity1
					<< "\tcold: " << _heatCapacity2);
}



WellDoubletControl* WellDoubletControl::create_wellDoubletControl(
				const char& selection)
{
	switch(selection)
	{
		case  'A':
			return new WellSchemeAC('A');
			break;
		case  'B':
			return new WellSchemeB('B');
			break;
		case  'C':
			return new WellSchemeAC('C');
			break;
		default:
			throw std::runtime_error("Well scheme not set");
	}
	return 0;
}


//////////////////////////////////////////////////////////


void WellSchemeAC::configureScheme()
{
	iterationState = searchingFlowrate;
	deltaTsign_stored = 0, factor = 0.1;  
		// initialization - for adapting flowrate

	if(scheme_identifier == 'A')  // T1 at warm well
	{
		LOG("\t\t\tconfigure scheme A");
		simulation_result_aiming_at_target =
			&WellSchemeAC::temperature_well1;
	} 
	else if(scheme_identifier == 'C')  // T1 - T2
	{
		LOG("\t\t\tconfigure scheme C");
		simulation_result_aiming_at_target =
			&WellSchemeAC::temperature_difference_well1well2;
	}

	if(operationType == storing)
	{
		LOG("\t\t\t\tfor storing");
		beyond = std::move(wdc::Comparison(
			new wdc::Greater(ACCURACY_TEMPERATURE_THRESHOLD)));
		notReached =  std::move(wdc::Comparison(
			new wdc::Smaller(ACCURACY_TEMPERATURE_THRESHOLD)));
	}
	else
	{
		LOG("\t\t\t\tfor extracting");
		beyond = std::move(wdc::Comparison(
			new wdc::Smaller(ACCURACY_TEMPERATURE_THRESHOLD)));
		notReached = std::move(wdc::Comparison(
			new wdc::Greater(ACCURACY_TEMPERATURE_THRESHOLD)));
	}
}

const WellDoubletControl::iterationState_t& 
	WellSchemeAC::evaluate_simulation_result(
		const double& _T1, const double& _T2,
		const double& _heatCapacity1, const double& _heatCapacity2)
{
	set_heatFluxes(_T1, _T2, _heatCapacity1, _heatCapacity2);
 
	double simulation_result_aiming_at_target = 
		(this->*(this->simulation_result_aiming_at_target))();
	// first adapt flow rate if temperature 1 at warm well is not
	// at target value
	if(iterationState == searchingFlowrate)
	{
		if(beyond(simulation_result_aiming_at_target, value_target))
		{
			if(fabs(result.Q_w - value_threshold) < ACCURACY_FLOWRATE_TARGET)
        		{  // cannot store / extract the heat
                		iterationState = searchingPowerrate;  
                		LOG("\t\t\tstop adapting flow rate");
			}
			else
			{
				adapt_flowrate();
			}
		}
		else if(notReached(
			simulation_result_aiming_at_target, value_target))
		{
			if(fabs(result.Q_w) < ACCURACY_FLOWRATE_TARGET)
          		{  // limited by flowrate
                		LOG("\t\t\tstop adapting flow rate");
                		iterationState = converged;  
				write_outputFile();
				LOG("\tconverged");
			}
			else
			{
				adapt_flowrate();
			}
		}
		else  // T1 at threshold
		{
			iterationState = converged;
			write_outputFile();
		}
	}

	if(iterationState == searchingPowerrate)
	{	
		// then adapt power rate if temperature 1 at warm well is beyond target
		// (and flow rate adaption as not succedded before)
		if(beyond(simulation_result_aiming_at_target, value_target))
		{
			adapt_powerrate();
		}
		else
		{
			iterationState = converged;
			write_outputFile();
			LOG("\tconverged");
		}
	}

	return iterationState;
}

void WellSchemeAC::set_flowrate()
{
	double temp, denominator = 
		heatCapacity1 * result.T1 - heatCapacity2 * result.T2;

	if(fabs(result.T1 - result.T2) < 1.e-10)
	{
		std::cout << "WARNING - well 1 is not warmer than well 2\n";
		if(operationType == storing)
		{
			LOG("\tset flowrate zero");
			result.Q_w = 0.;
		}
		else
		{
			LOG("\tset flowrate to threshold value");
			result.Q_w = value_threshold;
		}
	}
	else
	{
		temp = result.Q_H / denominator;

        	if(operationType == WellDoubletControl::storing)
			result.Q_w = wdc::confined(temp , 0., value_threshold);
		else
			result.Q_w = wdc::confined(temp, value_threshold, 0.);
      		
		LOG("\t\t\testimate flow rate\t" << result.Q_w);
	}
}

void WellSchemeAC::adapt_flowrate()
{
        double deltaT =  
		((this->*(this->simulation_result_aiming_at_target))() -
				value_target);
  
	// avoid that T1 jumps around threshold (deltaT flips sign)
	if(deltaTsign_stored != 0  // == 0: take initial value for factor
			&& deltaTsign_stored != wdc::sign(deltaT))
		factor *= FLOWRATE_ADAPTION_FACTOR;

	deltaTsign_stored = wdc::sign(deltaT);

        if(operationType == WellDoubletControl::storing)
        {
                result.Q_w = wdc::confined(result.Q_w * (1 + factor * deltaT),
							0., value_threshold);
        }
        else
        {	// decrease absolute value of flow rate
                result.Q_w = wdc::confined(result.Q_w * (1 - factor * deltaT),
							value_threshold, 0.);
        }
        LOG("\t\t\tadapt flow rate to\t" << result.Q_w);
}

void WellSchemeAC::adapt_powerrate()
{
        result.Q_H -= fabs(result.Q_w) * heatCapacity1 * (
                (this->*(this->simulation_result_aiming_at_target))() -
                        // Scheme A: T1, Scheme C: T1 - T2 
                value_target); 
        LOG("\t\t\tadapt power rate to\t" << result.Q_H);
        result.flag_powerrateAdapted = true;
}

void WellSchemeB::configureScheme()
{
	iterationState = searchingPowerrate;  // not used right now

	if(operationType == storing)
	{
		LOG("\t\t\tconfigure scheme B for storing");
		beyond = wdc::Comparison(new wdc::Greater(0.));
	}
	else
	{
		LOG("\t\t\tconfigure scheme B for extracting");
		beyond = wdc::Comparison(new wdc::Smaller(0.));
	}
}

const WellDoubletControl::iterationState_t& WellSchemeB::evaluate_simulation_result(
		const double& _T1, const double& _T2,
		const double& _heatCapacity1, const double& _heatCapacity2)
{
	set_heatFluxes(_T1, _T2, _heatCapacity1, _heatCapacity2);

	if(beyond(result.T1, value_threshold))
	{
		adapt_powerrate();
	}
	else
	{
		write_outputFile();
		iterationState = converged;
		LOG("\tconverged");
	}

	return iterationState;
}


void WellSchemeB::set_flowrate()
{
        result.Q_w = value_target;
        LOG("\t\t\tset flow rate\t" << result.Q_w);
}


void WellSchemeB::adapt_powerrate()
{
        result.Q_H -= fabs(result.Q_w) * heatCapacity1 * (
                                        result.T1 - value_threshold);
        LOG("\t\t\tadapt power rate\t" << result.Q_H);
        result.flag_powerrateAdapted = true;
}


