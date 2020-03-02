#ifndef OGS_CONTRAFLOW_H
#define OGS_CONTRAFLOW_H

#include <vector>
#include <list>
#include "contraflow.h"
#include <sstream>
#include <string>



class OGS_contraflow
{
public:
	OGS_contraflow(int _indicator, contra::PipingData _piping_data, contra::FluidData _fluid_data) :
		indicator(_indicator), piping_data(_piping_data), fluid_data(_fluid_data)
	{}
	~OGS_contraflow() { delete contraflow; }

	struct input_group_t
	{
		double time;
		int mode;
		double Q;
		double var;
		input_group_t(const double& _time, const int& _mode, const double& _Q, const double& _var) :
					time(_time), mode(_mode), Q(_Q), var(_var) {}
	};

	contra::Contraflow* get_Contraflow() { return contraflow; }
	std::list<input_group_t> get_input_list() { return input_list; }
	std::vector<contra::SegmentData> get_segment_data_vec() { return segment_data_vec; }
	std::vector<long> get_nodes_vec() { return nodes_vec; }
	void add_node(int node_nr) { nodes_vec.push_back(node_nr); }

	void add_segment_data_group(contra::SegmentData segment_data)
			{ segment_data_vec.push_back(segment_data); }
	void add_input_group(double _time, int _mode, double _Q, double _var)
			{ input_list.push_back(input_group_t(_time, _mode, _Q, _var)); }

	void initialize();
	contra::Result call_contraflow(const double& current_time, stru3::DVec T_s);

private:
	std::string _itos(int i) // convert int to string
	{
	    std::stringstream s;
	    s << i;
	    return s.str();
	}

	contra::Contraflow* contraflow;  // change to shared ptr
	int indicator;  // 0: U, 1: 2U, 2: coax
	std::vector<contra::SegmentData> segment_data_vec;
	contra::PipingData piping_data;
	contra::FluidData fluid_data;
	std::list<input_group_t> input_list;

	std::vector<long> nodes_vec;
};

#endif
