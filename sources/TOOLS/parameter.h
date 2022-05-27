#ifndef PARAMETER_H
#define PARAMETER_H

// ICs and BCs
const double c_temperature_storage_initial = 50;
const double c_temperature_upwindAquifer_storing = 10;
const double c_temperature_upwindAquifer_extracting = 50;

// aquifer parameter
const double c_heatCapacity = 5.e6;
const double c_porosity = 0.5;

// mesh and geometry
const int c_gridSize = 11;
const int c_heatExchanger_nodeNumber = 5;

// numerics
const int c_minNumberOfIterations = 3;
const int c_maxNumberOfIterations = 200;
const int c_numberOfTimeSteps = 10;
const double c_timeStepSize = 1.e2;
const double c_accuracy_temperature = 0.01;
const double c_accuracy_flowrate = 1.e-6;
const double c_accuracy_powerrate = 10;

#endif
