#include "OGS_contraflow.h"

void OGS_contraflow::initialize()
{
	contraflow = new contra::Contraflow(indicator, segment_data_vec, piping_data, fluid_data);
}

contra::Result OGS_contraflow::call_contraflow(const double& current_time, stru3::DVec T_s)
{

	// delete parameter list entry if old such that current entry is at begin of list (list entries have to be ordered in time)
	while(input_list.size()>0 && current_time > input_list.begin()->time)
		input_list.erase(input_list.begin());

	if(input_list.empty())
	{
		throw std::runtime_error("ERROR in contraflow - Input list empty");
	}

	if(input_list.front().Q > 10e-10)
		contraflow->calculate(input_list.front().Q, input_list.front().mode, input_list.front().var, T_s);
	return contraflow->get_result();

}
