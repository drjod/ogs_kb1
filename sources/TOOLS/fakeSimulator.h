#ifndef FAKESIMULATOR_H
#define FAKESIMULATOR_H

#include <string>
#include "wellDoubletControl.h"
#include "parameter.h"
#include "timer.h"



struct Simulator
{
	virtual ~Simulator() {}
	
	virtual const bool& get_flag_iterate() const = 0;
	virtual const wdc::WellDoubletControl* get_wellDoubletControl() const = 0;
	virtual void create_wellDoubletControl(const int& selection) = 0;

        virtual void initialize_temperatures() = 0;
	virtual void calculate_temperatures(
				const double& Q_H, const double& Q_W) = 0;
        virtual void update_temperatures() = 0;
	virtual double calculate_error() = 0;

	virtual void execute_timeStep(const double& well2_temperature) = 0;
	virtual void simulate(
		const int& wellDoubletControlScheme, const double& Q_H, 
		const double& value_target, const double& value_threshold) = 0;
};


class FakeSimulator : public Simulator
{
        double temperatures[c_gridSize]; 
        double temperatures_previousIteration[c_gridSize];  // to calculate error 
        double temperatures_previousTimestep[c_gridSize];

        wdc::WellDoubletControl* wellDoubletControl;
        bool flag_iterate;  // to convert threshold value into a target value

public:
	FakeSimulator() : wellDoubletControl(nullptr) {}
	~FakeSimulator() 
	{ if(wellDoubletControl != nullptr) delete wellDoubletControl; }
			// a wellDoubletControl instance is constructed 
			// (and destructed) each time step

	const bool& get_flag_iterate() const override { return flag_iterate; }
	const wdc::WellDoubletControl* get_wellDoubletControl() const override
	{ return wellDoubletControl; }
	void create_wellDoubletControl(const int& selection) override;
				// is done at the begiining of each time step

        void initialize_temperatures() override;
	void calculate_temperatures(const double& Q_H, const double& Q_W);
						// solves advection equation
        void update_temperatures() override;
	double calculate_error() override;

	void execute_timeStep(const double& well2_temperature) override;
	void simulate(const int& wellDoubletControlScheme, const double& Q_H, 
		const double& value_target, const double& value_threshold) override;
		// values are passed to execute_timeStep(), they are constant
		// now but will be timestep-dependent in a real application
	template <typename T> void log_file(T toLog);

	friend std::ostream& operator<<(std::ostream& stream,
					const FakeSimulator& simulator);
};


#endif
