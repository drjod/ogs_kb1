/**************************************************************************
 FEMLib - Object: Source Terms ST
 Task:
 Programing:
 01/2004 OK Implementation
 last modified
 **************************************************************************/

#include "rf_fct.h"

#include "fem_ele_std.h"
#include "makros.h"
// C++ STL
//#include <fstream>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <numeric>
#include <limits>
#include <set>
#include <stdexcept>

#include "display.h"
#include "files0.h"
#include "mathlib.h"

#include "MshEditor.h" //NB
#include "PointWithID.h" // NB

// GeoSys-GeoLib
#include "GEOObjects.h"
#include "SensorData.h"
//#include "GeoType.h"

// GeoSys-MshLib
#include "fem_ele.h"

#include "tools.h"                                //GetLineFromFile
/* Tools */
#if !defined(USE_PETSC) && !defined(NEW_EQS) // && defined(other parallel libs)//03~04.3012. WW
//#ifndef NEW_EQS                                   //WW. 06.11.2008
#include "matrix_routines.h"
#endif

#ifdef NEW_EQS
#include "equation_class.h"
using Math_Group::CSparseMatrix;
#endif

// GeoSys-FEMLib
//OK_IC #include "rfsousin.h"
#include "rf_st_new.h"
#include "rf_tim_new.h"

// Math
#include "matrix_class.h"

// BaseLib
#include "FileTools.h"

//#include "pcs_dm.h"

// FEM
//#include "problem.h"
// For analytical source terms
#include "rf_mfp_new.h"
#include "rf_msp_new.h"
#include "rf_mmp_new.h"
#include "rf_node.h"
#include "rfmat_cp.h"


// Base
#include "quicksort.h"

// MathLib
#include "InterpolationAlgorithms/InverseDistanceInterpolation.h"
#include "InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

// FileIO
#include "FEMIO/GeoIO.h"
#include "FEMIO/ProcessIO.h"
#include "readNonBlankLineFromInputStream.h"
#include "XmlIO/RapidXMLInterface.h"

#include "SourceTerm.h"


using FiniteElement::CElement;
using MeshLib::CElem;
using MeshLib::CEdge;
using MeshLib::CNode;
using Math_Group::vec;

#ifndef GRAVITY_CONSTANT
#define GRAVITY_CONSTANT 9.81
#endif

#if defined(USE_MPI)
bool global_flag_keep_values = false;  // to keep ST values and make them consistent over all processors
#endif

int scaling_node_group_running = 0;

std::vector<CSourceTerm*> st_vector;
std::list<CSourceTermGroup*> st_group_list;
std::vector<std::string> analytical_processes;
std::vector<std::string> analytical_processes_polylines;
std::vector<NODE_HISTORY*> node_history_vector;   //CMCD
/**************************************************************************
 FEMLib-Method:
 Task: ST constructor
 Programing:
 01/2004 OK Implementation
 **************************************************************************/
CSourceTerm::CSourceTerm() :
	ProcessInfo(), GeoInfo(), _coupled (false), _sub_dom_idx(-1), dis_linear_f(NULL), GIS_shape_head(NULL), _distances(NULL)
                                                  // 07.06.2010, 03.2010. WW
{
	geoInfo_connected = new GeoInfo();
	geoInfo_threshold = new GeoInfo();
	geoInfo_storageRateInlet = new GeoInfo();
	geoInfo_storageRateOutlet = new GeoInfo();

	//geoInfo_wellDoublet_HE = new GeoInfo();
	geoInfo_wellDoublet_well1_aquifer = new GeoInfo();
	geoInfo_wellDoublet_well2_aquifer = new GeoInfo();
	geoInfo_wellDoublet_well1_liquidBC = new GeoInfo();
	geoInfo_wellDoublet_well2_liquidBC = new GeoInfo();
	ogs_WDC = nullptr; // JOD 2018-08-09
	ogs_contraflow = nullptr;

   CurveIndex = -1;
   //KR critical_depth = false;
   //	COUPLING_SWITCH = false;
   geo_node_value = 0.0;
   nodes = NULL;                                  //OK
   analytical = false;                            //CMCD
   pressureBoundaryCondition = false;
   //  display_mode = false; //OK
   TimeInterpolation = 0;                   //BG
   time_contr_curve = -1;                //SB
   time_contr_function = "";
   _isConstrainedST = false;

   everyoneWithEveryone = false;  //  JOD 2015-11-18
   threshold.type = Threshold::no;
   storageRate.apply = false;
   assign_to_element_edge = false;

   connected_geometry = false;
   connected_geometry_verbose_level = 0;
   connected_geometry_exchange_term = 1.;  // value 1 used as default to use it with BOREHOLE
   connected_geometry_offset = 0.;
   connected_geometry_mode = -1;
   connected_geometry_minimum_velocity_abs = -1;               // JOD 2015-11-18
   connected_geometry_ref_element_number = -1;
   connected_geometry_reference_direction[0] = connected_geometry_reference_direction[1] = connected_geometry_reference_direction[2] = 0;
   connected_geometry_couplingType = 1; // matrix couplig as default

   threshold_geometry = false;
   storageRate_geometry = false;
   ignore_axisymmetry = false;

   wdc_connector_materialGroup = -1;
   wdc_flag_extract_and_reinject = false;

   variable_storage = false;

   scaling_mode = 0;
   scaling_verbosity = 0;
   keep_values = false;
   borehole_mode = -1;
   borehole_modified_aquifer_parameter = -1.;
   average_mode = -1;
   verbosity = 0;
   pcs_type_name_cond2 = "";
}

// KR: Conversion from GUI-ST-object to CSourceTerm
CSourceTerm::CSourceTerm(const SourceTerm* st)
	: ProcessInfo(st->getProcessType(),st->getProcessPrimaryVariable(),NULL),
	  GeoInfo(st->getGeoType(),st->getGeoObj()),
	  DistributionInfo(st->getProcessDistributionType()),
	  _distances(NULL)
{
	geoInfo_connected = new GeoInfo();
	geoInfo_threshold = new GeoInfo();

	setProcess( PCSGet( this->getProcessType() ) );
	this->geo_name = st->getGeoName();
	const std::vector<size_t> dis_nodes = st->getDisNodes();
	const std::vector<double> dis_values = st->getDisValues();

	if (this->getProcessDistributionType() == FiniteElement::CONSTANT_NEUMANN)
	{
		this->geo_node_value = dis_values[0];
	}
	else if (this->getProcessDistributionType() == FiniteElement::LINEAR_NEUMANN)
	{
		for (size_t i=0; i<dis_values.size(); i++)
		{
			this->PointsHaveDistribedBC.push_back(static_cast<int>(dis_nodes[i]));
			this->DistribedBC.push_back(dis_values[i]);
		}
	}
	else if (this->getProcessDistributionType() == FiniteElement::DIRECT
			|| this->getProcessDistributionType() == FiniteElement::RECHARGE_DIRECT)
	{
		// variable "fname" needs to be set, this must be done from outside!
	}
	else
		std::cout << "Error in CBoundaryCondition() - DistributionType \""
		          << FiniteElement::convertDisTypeToString(this->getProcessDistributionType())
				  << "\" currently not supported." << "\n";
}

/**************************************************************************
 FEMLib-Method:
 Task: BC deconstructor
 Programing:
 04/2004 OK Implementation
 **************************************************************************/
CSourceTerm::~CSourceTerm()
{
	delete geoInfo_connected;
	delete geoInfo_threshold;
	delete geoInfo_storageRateInlet;
	delete geoInfo_storageRateOutlet;
	//delete geoInfo_wellDoublet_HE;
	delete geoInfo_wellDoublet_well1_aquifer;
	delete geoInfo_wellDoublet_well2_aquifer;
	delete geoInfo_wellDoublet_well1_liquidBC;
	delete geoInfo_wellDoublet_well2_liquidBC;

	delete _distances;
	for (size_t i=0; i<this->_weather_stations.size(); i++) // KR / NB clear climate data information
	   delete this->_weather_stations[i];

   DeleteHistoryNodeMemory();
   //    dis_file_name.clear();
   node_number_vector.clear();
   node_value_vector.clear();
   node_renumber_vector.clear();
   PointsHaveDistribedBC.clear();
   DistribedBC.clear();
   element_st_vector.clear();
   //WW----------22.02.2007-------------------
   // TF 06/2010
   size_t size(normal2surface.size());
   for (size_t i = 0; i < size; i++)
      delete normal2surface[i];
   size = pnt_parameter_vector.size();
   for (size_t i = 0; i < size; i++)
      delete pnt_parameter_vector[i];
   if(GIS_shape_head)                             // 07.06.2010. WW
   {
      delete [] GIS_shape_head;
      GIS_shape_head = NULL;
   }
   //WW
   if(dis_linear_f) delete dis_linear_f;
   dis_linear_f = NULL;

   //WW---------------------------------------
}


const std::string& CSourceTerm::getGeoName() const
{
   return geo_name;
}


double CSourceTerm::getCoupLeakance () const
{
   return _coup_leakance;
}


/**************************************************************************
 FEMLib-Method:
 Task: ST read function
 Programing:
 01/2004 OK Implementation
 11/2004 MB neues read Konzept
 02/2005 MB River condition
 03/2005 WW Node force released by excavation
 11/2005 CMCD Analytical source for matrix
 04/2006 OK CPL
 04/2006 OK MSH_TYPE
06/2010 TF modification of the signature, added geo_obj and unique_name
**************************************************************************/
std::ios::pos_type CSourceTerm::Read(std::ifstream *st_file,
		const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name)
{
   char line[MAX_ZEILE];
   std::string line_string, sub_string;
   bool new_keyword = false;

   std::stringstream in;
                                                  // JOD 
   channel = 0, node_averaging = 0, air_breaking = false;
   no_surface_water_pressure = 0, explicit_surface_water_pressure = false;
   distribute_volume_flux = false;
   std::ios::pos_type position;

   // read loop
   while (!new_keyword)
   {
      position = st_file->tellg();
      if (!GetLineFromFile(line, st_file))
         break;
      line_string = line;
      if (line_string.find("#") != std::string::npos)
      {
         new_keyword = true;
         break;
      }
      remove_white_space(&line_string);           //OK

      /* search for keywords */
                                                  // subkeyword found
      if (line_string.find("$PCS_TYPE") != std::string::npos)
      {
    	  FileIO::ProcessIO::readProcessInfo (*st_file, _pcs_type);
//         in.str(GetLineFromFile1(st_file));
//         std::string tmp;
//         in >> tmp;
//         setProcessType (convertProcessType (tmp));
//         in.clear();
         continue;
      }

                                                  // subkeyword found
      if (line_string.find("$PRIMARY_VARIABLE") != std::string::npos)
      {
    	  in.str (readNonBlankLineFromInputStream (*st_file));
//         in.str(GetLineFromFile1(st_file));
         std::string tmp;
         in >> tmp;

         if ( this->getProcessType() == FiniteElement::MASS_TRANSPORT )
         {
             // HS set the pointer to MCP based on component name.
             // first do a check whether this name is existing and unique.
             if ( cp_name_2_idx.count( tmp ) == 1 )
             {
                 setProcess(cp_vec[cp_name_2_idx[tmp]]->getProcess() );
                 setProcessPrimaryVariable( FiniteElement::CONCENTRATION );
             }
             else
             {
                 DisplayErrorMsg("Error: In reading ST file, the input component names are not found in MCP file!!!");
                 exit(1);
             }
         }
         else
         {
             setProcess( PCSGet( this->getProcessType() ) );
             setProcessPrimaryVariable (FiniteElement::convertPrimaryVariable (tmp));
         }
         in.clear();
         continue;
      }

      if (line_string.find("$COMP_NAME") != std::string::npos)
      {
    	  in.str(readNonBlankLineFromInputStream (*st_file));
//         in.str(GetLineFromFile1(st_file));
         std::string tmp;
         in >> tmp;
         // HS set the pointer to MCP based on component name.
         // first do a check whether this name is existing and unique.
         if ( cp_name_2_idx.count( tmp ) == 1 )
         {
             setProcess(cp_vec[cp_name_2_idx[tmp]]->getProcess() );
             setProcessPrimaryVariable( FiniteElement::CONCENTRATION );
         }
         else
         {
             DisplayErrorMsg("Error: In reading ST file, the input component names are not found in MCP file!!!");
             exit(1);
         }
         in.clear();
         continue;
      }

      if (line_string.find("$GEO_TYPE") != std::string::npos)
      {
         ReadGeoType(st_file, geo_obj, unique_name);
         continue;
      }

                                                  //05.09.2008 WW
      if (line_string.find("$DIS_TYPE") != std::string::npos)
      {
         //10.04.2008. WW  if(line_string.compare("$DIS_TYPE")==0) {
         if (line_string.find("CONDITION") != std::string::npos)
         {
            _coupled = true;
            ReadDistributionType(st_file);
            in.str(readNonBlankLineFromInputStream(*st_file));
            in >> line_string >> pcs_type_name_cond >> pcs_type_name_cond2;
            in.clear();
            in.str(readNonBlankLineFromInputStream(*st_file));    //
            in >> pcs_pv_name_cond;
            in.clear();
//            in.str(GetLineFromFile1(st_file));
            in.str(readNonBlankLineFromInputStream(*st_file));
            in >> _coup_leakance >> st_rill_height >> coup_given_value >> coup_residualPerm;
            in.clear();
         }                                        //05.09.2008 WW
         else
         {
            ReadDistributionType(st_file);
            continue;
         }
      }

      if (line_string.find("$NODE_AVERAGING") != std::string::npos)
      {
         in.clear();
         node_averaging = true;
         continue;
      }

	  if (line_string.find("$DISTRIBUTE_VOLUME_FLUX") != std::string::npos) //  JOD 5.3.07
      {
         in.clear();
         distribute_volume_flux = true;
         continue;
      }
                                               
      if (line_string.find("$NEGLECT_SURFACE_WATER_PRESSURE") != std::string::npos)
      {       // JOD 4.10.01
         in.clear();
         no_surface_water_pressure = true;
         continue;
      }

	  if (line_string.find("$EXPLICIT_SURFACE_WATER_PRESSURE") != std::string::npos) 
      {       // JOD 5.3.07
         in.clear();
         explicit_surface_water_pressure = true;
         continue;
      }

      if (line_string.find("$CHANNEL") != std::string::npos)
      {
         in.clear();
         in.str(readNonBlankLineFromInputStream(*st_file));
         in >> channel_width;
         channel = 1;
         continue;
      }

	  if (line_string.find("$AIR_BREAKING") != std::string::npos) // JOD 5.3.07
      {
         in.clear();
         in.str(readNonBlankLineFromInputStream(*st_file));
         in >> air_breaking_factor >> air_breaking_capillaryPressure >> air_closing_capillaryPressure;
         continue;
      }

      if (line_string.find("$TIM_TYPE") != std::string::npos)
      {
         in.str(readNonBlankLineFromInputStream(*st_file));
         in >> tim_type_name;
         if (tim_type_name.find("CURVE") != std::string::npos)
         {
        	 //				dis_type = 0;
            in >> CurveIndex;
         }
         in.clear();
         continue;
      }

  	  //defines if time dependent source terms are use as piecewise constant or linear interpolated; BG 05/2011
      if (line_string.find("$TIME_INTERPOLATION") != std::string::npos)
      {
         in.str(GetLineFromFile1(st_file));
         in >> interpolation_method;
         if (interpolation_method.find("LINEAR") != std::string::npos)
         {
            this->TimeInterpolation = 0;
         }
         if (interpolation_method.find("PIECEWISE_CONSTANT") != std::string::npos)
         {
            this->TimeInterpolation = 1;
         }
         in.clear();
         continue;
      }

      if (line_string.find("$FCT_TYPE") != std::string::npos)
      {
         in.str(readNonBlankLineFromInputStream(*st_file));
         in >> fct_name;                          //sub_line
                                                  //WW
         if (fct_name.find("METHOD") != std::string::npos)
         in >> fct_method;
         in.clear();
      }
  		//....................................................................
		//active state of the bc is time controlled  SB
		if (line_string.find("$TIME_CONTROLLED_ACTIVE") != std::string::npos)
		{
			in.str(readNonBlankLineFromInputStream(*st_file));
			std::string test;
			in >> test;
			if ( test.compare("CURVE")== 0){ // Curve number
				in >> time_contr_curve;
			}
			else if(test.compare("FUNCTION")==0){
				// std::cout << " CSourceTerm::Read():Time_CONTROLLED_ACTIVE: Reading FUNCTION name " << "\n";
				in >> time_contr_function;
			}
			else
				std::cout << " CSourceTerm::Read():Time_CONTROLLED_ACTIVE: Read Error; Give CURVE and number or FUNCTION and name " << "\n";
			in.clear();
		}

      if (line_string.find("$MSH_TYPE") != std::string::npos)
      {
         in.str(readNonBlankLineFromInputStream(*st_file));
         std::string sub_string;
         in >> sub_string;                        //sub_line
         msh_type_name = "NODE";
         if (sub_string.find("NODE") != std::string::npos)
         {
            in >> msh_node_number;
            in.clear();
         }
         continue;
      }

	  if (line_string.find("$CONSTRAINED") != std::string::npos)
	  {
		  Constrained temp;

		  _isConstrainedST = true;
		  in.str(readNonBlankLineFromInputStream(*st_file));
		  std::string tempst;

		  in >> tempst;	//PROCESS_TYPE associated with PRIMARY_VARIABLE
		  temp.constrainedProcessType = FiniteElement::convertProcessType(tempst);
		  if (!(temp.constrainedProcessType == FiniteElement::MASS_TRANSPORT ||
			  temp.constrainedProcessType == FiniteElement::HEAT_TRANSPORT ||
			  temp.constrainedProcessType == FiniteElement::LIQUID_FLOW ||
			  temp.constrainedProcessType == FiniteElement::RICHARDS_FLOW)) {
			  _isConstrainedST = false;
			  break;
		  }

		  in >> tempst;	//PRIMARY_VARIABLE to be constrained
		  temp.constrainedPrimVar = FiniteElement::convertPrimaryVariable(tempst);

		  in >> temp.constrainedValue;	//Constrained Value

		  in >> tempst;	//Constrain direction (greater/smaller than value)
		  temp.constrainedDirection = convertConstrainedType(tempst);
		  temp.constrainedVariable = ConstrainedVariable::INVALID_CONSTRAINED_VARIABLE;
		  if (!(temp.constrainedDirection == ConstrainedType::SMALLER || temp.constrainedDirection == ConstrainedType::GREATER))
		  {
			  std::cout << "No valid constrainedDirection for " << FiniteElement::convertProcessTypeToString(temp.constrainedProcessType)
				  << " (" << tempst << ")" << "\n";
			  _isConstrainedST = false;
		  }

		  in >> tempst;	//full constrain option 
		  if (tempst == "COMPLETE_CONSTRAIN")
		  {
			  temp._isCompleteConstrained = true;
		  }

		  if (_isConstrainedST)
			  this->_constrainedST.push_back(temp);
		  in.clear();
	  }
// CB JOD MERGE //
	  //....................................................................
	  if (line_string.find("$CONNECTED_GEOMETRY") != std::string::npos) // SB 02/2015    JOD 2015-11-18
	  {
		  //ReadGeoType(st_file, geo_obj, unique_name);

	      FileIO::GeoIO::readGeoInfo(geoInfo_connected, *st_file, connected_geometry_name, geo_obj, unique_name);
	      //in.str(readNonBlankLineFromInputStream(*st_file));
		  //in >> connected_geometry_type >> connected_geometry_name;
		  this->connected_geometry = true;
		  in.clear();
		  continue;
	  }
	  //....................................................................
	  if (line_string.find("$THRESHOLD_GEOMETRY") != std::string::npos) // JOD 2018-02-20
	  {
	      FileIO::GeoIO::readGeoInfo(geoInfo_threshold, *st_file, threshold_geometry_name, geo_obj, unique_name);
	      //in.str(readNonBlankLineFromInputStream(*st_file));
		  //in >> connected_geometry_type >> connected_geometry_name;
		  this->threshold_geometry = true;
		  in.clear();
		  continue;
	  }
	  //....................................................................
	  if (line_string.find("$STORAGE_RATE_GEOMETRY") != std::string::npos) // JOD 2018-02-22
	  {
	      FileIO::GeoIO::readGeoInfo(geoInfo_storageRateInlet, *st_file,
	    		  	  storageRate.inlet_geometry_name, geo_obj, unique_name);
	      FileIO::GeoIO::readGeoInfo(geoInfo_storageRateOutlet, *st_file,
	    		  	  storageRate.outlet_geometry_name, geo_obj, unique_name);

		  this->storageRate_geometry = true;
		  in.clear();
		  continue;
	  }
	  //....................................................................
	  if (line_string.find("$CONNECT_PARAMETER") != std::string::npos) //  JOD 2015-11-18
	  {
		  in.str(readNonBlankLineFromInputStream(*st_file));
		  in >> connected_geometry_exchange_term >> connected_geometry_verbose_level; // >> connected_geometry_offset;
		  this->connected_geometry = true;
		  in.clear();
		  continue;
	  }
	  //....................................................................
	  if (line_string.find("$CONNECT_MODE") != std::string::npos)  //  JOD 2015-11-18
	  {
		  in.str(readNonBlankLineFromInputStream(*st_file));
		  in >> connected_geometry_mode;  // 0: NNNC symmetric, 1: NNNC non-symmetric (downwind fixed), 2 NNNC non-symmetric (downwind)
		  in >> connected_geometry_couplingType; // 0: RHS, 1: matrix entry
		  if ((connected_geometry_mode == 2))
			  in >> connected_geometry_ref_element_number >> connected_geometry_reference_direction[0] >>
			  connected_geometry_reference_direction[1] >> connected_geometry_reference_direction[2] >>
			  connected_geometry_minimum_velocity_abs;
		  in.clear();

		  if(connected_geometry_mode == 0)
			  std::cout << "Symmetric NNNC -";
		  else if(connected_geometry_mode == 1)
			  std::cout << "Non-symmetric NNNC -";
		  else if(connected_geometry_mode == 2)
		  {
			  std::cout << "Non-symmetric NNNC with element " << connected_geometry_ref_element_number << " and direction "
			  	  << connected_geometry_reference_direction[0] << " "
				  << connected_geometry_reference_direction[1] << " "
				  << connected_geometry_reference_direction[2] << " - ";
		  }
		  else
			  throw std::runtime_error("Error - Connected geometry mode not supported");

		  if(connected_geometry_couplingType == 0)
			  std::cout << " via RHS\n";
		  if(connected_geometry_couplingType == 1)
			  std::cout << " via matrix\n";  // default
		  else
			  throw std::runtime_error("Error - Connected geometry coupling type not supported");

		  continue;
	  }
	  //....................................................................
	  if (line_string.find("$EVERYONE_WITH_EVERYONE") != std::string::npos)
	  {       //  JOD 2015-11-18
		  in.clear();
		  everyoneWithEveryone = true;
		  continue;
	  }
	  //....................................................................
	  if (line_string.find("$THRESHOLD") != std::string::npos)
	  {       //  JOD 2018-1-31
		  int threshold_type, threshold_scheme;

		  in.str(readNonBlankLineFromInputStream(*st_file));
		  in >> threshold_type >> threshold.process >> threshold.value >> threshold_scheme;

		  threshold.type = static_cast<Threshold::Type>(threshold_type);  //  enum Type { no, lower, upper};
		  threshold.scheme = static_cast<Threshold::Scheme>(threshold_scheme);  //  enum Scheme {_explicit, _implicit};

		  if (threshold.scheme == Threshold::_implicit)
			  in >> threshold.delta;  // smoothed threshold
		  in >> threshold.verbosity;
		  in.clear();
		  continue;
	  }
	  //....................................................................
	  if (line_string.find("$STORAGE_RATE") != std::string::npos)
	  {       //  JOD 2018-1-31
		  in.str(readNonBlankLineFromInputStream(*st_file));
	  	  in >> storageRate.process >> storageRate.inputValue >> storageRate.absMaximum >> storageRate.verbosity;
		  in.clear();
		  storageRate.apply = true;
	  	  continue;
	  }
	  //....................................................................
	  if (line_string.find("$WELL_DOUBLET_GEOMETRY") != std::string::npos) // JOD 2018-06-13
	  {
	      //FileIO::GeoIO::readGeoInfo(geoInfo_wellDoublet_HE, *st_file,
	    	//	  	  well1_geometry_name_HE, geo_obj, unique_name);
		  if(getGeoType() == GEOLIB::POINT)
		  {
			  FileIO::GeoIO::readGeoInfo(geoInfo_wellDoublet_well1_aquifer, *st_file,
						  well1_geometry_name_aquifer, geo_obj, unique_name);
			  FileIO::GeoIO::readGeoInfo(geoInfo_wellDoublet_well1_liquidBC, *st_file,
						  well1_geometry_name_liquidBC, geo_obj, unique_name);
			  FileIO::GeoIO::readGeoInfo(geoInfo_wellDoublet_well2_aquifer, *st_file,
						  well2_geometry_name_aquifer, geo_obj, unique_name);
			  FileIO::GeoIO::readGeoInfo(geoInfo_wellDoublet_well2_liquidBC, *st_file,
						  well2_geometry_name_liquidBC, geo_obj, unique_name);
		  }
		  else if(getGeoType() == GEOLIB::POLYLINE)
		  {
			  FileIO::GeoIO::readGeoInfo(geoInfo_wellDoublet_well2_aquifer, *st_file,
			  						  well2_geometry_name_aquifer, geo_obj, unique_name);
			  wdc_flag_extract_and_reinject = true;
		  }
		  else
			  throw std::runtime_error("Error in WDC - GeoType unknown or not supported");
		  in.clear();
		  continue;
	  }
	  //....................................................................
	  if (line_string.find("$WELL_DOUBLET_PARAMETER") != std::string::npos) // JOD 2018-06-13
	  {
		  std::cout << "Read WDC parameter\n";
		  int numberOfParameterSets;
		  double accuracy_temperature = 0.01, accuracy_powerrate = 1e3, accuracy_flowrate = 1.e-5;  // default values
		  double well_shutdown_temperature_range = 5.;
		  connected_geometry_mode = 1;  // NNNC downwind fixed
		  connected_geometry_couplingType = 1; // NNNC via matrix entry as default

		  std::string tmp;
		  int heat_pump_flag;
		  std::string heat_pump_file_name;  // for carnot

		  in.str(readNonBlankLineFromInputStream(*st_file));
		  in >> tmp >> heat_pump_flag >> heat_pump_file_name;
		  if(heat_pump_flag == 0)
		  {
			  std::cout << "\tNo heat pump";
		  }
		  else if(heat_pump_flag == 1)
		  {
			  in >> heat_pump_file_name;
			  std::cout << "\tCarnot heat pump - File name: " << heat_pump_file_name <<  "\n";
		  }
		  else
			  throw std::runtime_error("ERROR in heat pump specification");
		  in.clear();

		  bool logging = false;

		  in.str(readNonBlankLineFromInputStream(*st_file));
		  in >> numberOfParameterSets >> well_shutdown_temperature_range
		  	  >> accuracy_temperature >> accuracy_powerrate >> accuracy_flowrate
			  >> connected_geometry_couplingType  // as in $CONNECT_MODE
			  >> logging;
		  in.clear();

		  if(CRFProcess* m_pcs = PCSGet(convertProcessTypeToString(getProcessType())))
		  {
			  ogs_WDC = new OGS_WDC(well_shutdown_temperature_range,
					  accuracy_temperature, accuracy_powerrate, accuracy_flowrate,
					  m_pcs->ogs_WDC_vector.size(), logging);

			  ogs_WDC->set_heat_pump_parameter(heat_pump_flag, heat_pump_file_name);

			  while(numberOfParameterSets--)  // no check if number is right
			  {
				  double tmp0, tmp2, tmp3, tmp4, tmp5=-1.;
				  int tmp1;

				  in.str(readNonBlankLineFromInputStream(*st_file));
				  in >> tmp0 >> tmp1 >> tmp2 >> tmp3 >> tmp4 >> tmp5;
				  in.clear();
				  ogs_WDC->add_parameterGroup
				  (
					tmp0,  // time
					tmp1,  // scheme indicator [0, 1, 2]
					tmp2,  // powerrate
					tmp3,  // target value
					tmp4,  // threshold value
					tmp5   // temperature sink for heat pump
				  );
			  }

			  m_pcs->ogs_WDC_vector.push_back(ogs_WDC);
		  }
		  else
			  throw std::runtime_error("No PCS for WellDoubletControl");

		  continue;
	  }
	  //....................................................................
	  if (line_string.find("$WELL_DOUBLET_CONNECTOR") != std::string::npos) // JOD 2018-10-25
	  {
		  in.str(readNonBlankLineFromInputStream(*st_file));
		  in >> wdc_connector_materialGroup;
		  in >> wdc_connector_normaldirectionVector[0] >> wdc_connector_normaldirectionVector[1]
					>> wdc_connector_normaldirectionVector[2];
		  in.clear();
		  continue;
	  }
	  //....................................................................
	  if (line_string.find("$SCALING") != std::string::npos) // JODSH 2021-11-04
	  {
		  in.str(readNonBlankLineFromInputStream(*st_file));
		  in >> scaling_mode >> scaling_verbosity;  // 1: with permeability, 2: with permeability and viscosity
		  std::cout << "\tScaling mode " << scaling_mode;
		  if(scaling_mode == 0)
			  std::cout << " - No scaling\n";
		  else if(scaling_mode == 1)
			  std::cout << " - With horizontal permeability\n";
		  else if(scaling_mode == 2)
			  std::cout << " - With horizontal permeability and viscosity\n";
		  else
			  throw std::runtime_error("Scaling mode not supported");

		  scaling_node_group = scaling_node_group_running;
		  scaling_node_group_running++;
		  in.clear();

		  continue;
	  }
	  //....................................................................
	  if (line_string.find("$CONTRAFLOW_PIPES") != std::string::npos) // JOD 2019-07-30
	  {
		  std::cout << "CONTRAFLOW_PIPES\n";
		 int tmp0;
		 double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11;

		 in.str(readNonBlankLineFromInputStream(*st_file));  // pipe
		 in >> tmp0 >> tmp1 >> tmp2 >> tmp3 >> tmp4 >> tmp5 >> tmp6 >> tmp7;
		 in.clear();

		 in.str(readNonBlankLineFromInputStream(*st_file));  // fluid
		 in >> tmp8 >> tmp9 >> tmp10 >> tmp11;
		 in.clear();

		 ogs_contraflow = new OGS_contraflow(tmp0, // indicator - 0: U, 1: 2U, 2: coax
				 {	// pipe
						tmp1,  // d_0_i
						tmp2,  // d_0_o
						tmp3,  // d_1_i
						tmp4,  // d_1_o
						tmp5,  // w
						tmp6,  // lambda_0
						tmp7  // lambda_1
				 },
				 {	// fluid
						 tmp8,  // lambda
						 tmp9,	// mu
						 tmp10,	// c
						 tmp11	// rho
				 }

		 );

		 int numberOfSegments;
		 in.str(readNonBlankLineFromInputStream(*st_file));
		 in >> numberOfSegments;
		 in.clear();

		 while(numberOfSegments--)  // no check if number is right
		 {
			  in.str(readNonBlankLineFromInputStream(*st_file));
			  in >> tmp0 >> tmp1 >> tmp2 >> tmp3;
			  in.clear();
			  //new_ogs_WellDoubletControl.wellDoubletData.parameter_list.emplace_back(
			  ogs_contraflow->add_segment_data_group
			  ({
				tmp0,  // N
				tmp1,  // L
				tmp2,  // D
				tmp3  // lambda_g
			  });
		 }
		ogs_contraflow->initialize();


		 if(CRFProcess* m_pcs = PCSGet(convertProcessTypeToString(getProcessType())))
		 {
			 m_pcs->ogs_contraflow_vector.push_back(ogs_contraflow);
		 }
		 else
			 throw std::runtime_error("No PCS for WellDoubletControl");

		 std::cout << "Set contraflow source term\n";

		 continue;
	  }
	  //....................................................................
	  if (line_string.find("$CONTRAFLOW_INPUT") != std::string::npos) // JOD 2019-07-30
	  {
		 std::cout << "CONTRAFLOW_INPUT\n";
		  int numberOfInputSets;

		  in.str(readNonBlankLineFromInputStream(*st_file));
		  in >> numberOfInputSets;

		  in.clear();
		  while(numberOfInputSets--)  // no check if number is right
		  {
			  int tmp1;
			  double tmp0, tmp2, tmp3;

			  in.str(readNonBlankLineFromInputStream(*st_file));
			  in >> tmp0 >> tmp1 >> tmp2 >> tmp3;
			  in.clear();

			  ogs_contraflow->add_input_group
			  (
				tmp0,  // time
				tmp1,	// mode 0: provide feed in temperature T_in, 1: provide temperature difference dT
				tmp2,  // Q
				tmp3  // T_in / dT
			  );
		  }
		  continue;
	  }
	  //....................................................................
	  if (line_string.find("$ASSIGN_TO_ELEMENT_EDGE") != std::string::npos)
	  {
		 in.clear();
		 assign_to_element_edge = true;
		 continue;
	  }
	  //....................................................................
	  if (line_string.find("$VARIABLE_STORAGE") != std::string::npos)
	  {
		 in.clear();
		 variable_storage = true;
		 continue;
		 std::cout << "VARIABLE STORAGE\n";
		 std::cout << "\tWARNING: Time step size form LIQUID_FLOW used\n";
	  }
	 //....................................................................
      if (line_string.find("$IGNORE_AXISYMMETRY") != std::string::npos)
      {
         in.clear();
         ignore_axisymmetry = true;
         continue;
      }
      //....................................................................
      if (line_string.find("$KEEP_VALUES") != std::string::npos)
      {
         in.clear();
         keep_values = true;
#if defined(USE_MPI)
	  global_flag_keep_values = true;
#endif
	  std::cout << "\tKeep values\n";
         continue;
      }
      //....................................................................
	  if (line_string.find("$BOREHOLE") != std::string::npos && // JOD 2021-12-06
			  line_string.find("GIVEN_VALUE") == std::string::npos &&
			  line_string.find("CONDUCTIVITY") == std::string::npos)
	  {
		in.str(readNonBlankLineFromInputStream(*st_file));
		in >> borehole_mode;

		switch(borehole_mode)
		{
			case -1:
				std::cout << " - Switched off";
				break;
			case 0: // conductive - peaceman
				// use the following :
				// $CONNECT_PARAMETERS
 				//   1. 0  ; exchange coefficient (is 1. for borehole in mode 0),   verbosity level = 0, 1, 1
				// $CONNECT_MODE
				//  0   ; symmetric
				in >> borehole_data.radius >> connected_geometry_verbose_level >> connected_geometry_couplingType;
				std::cout << "\tBorehole - Peaceman connection - Radius " << borehole_data.radius;
				connected_geometry_mode = 0; // NNNC symmetric
				break;
			case 1:  // advective (currently only heat) transport - with LIQUID_FLOW source / sink term
				in >> connected_geometry_verbose_level >> connected_geometry_couplingType;
				std::cout << "\tBorehole - advective connection";
				connected_geometry_mode = 1; // NNNC non-symmetric
				break;
			default:
				throw std::runtime_error("Borehole mode not supported");
		}

		  if(connected_geometry_couplingType == 0)
			  std::cout << " - coupling via RHS\n";
		  else if(connected_geometry_couplingType == 1)  // default
			  std::cout << " - coupling via matrix\n";
		  else
			  throw std::runtime_error("Error - Connected geometry coupling type not supported");
		
		connected_geometry = true;  // connected by NNNC (also if not used with $CONNECTED_GEOMETRY)

		in.clear();
		continue;
	  }
      	  //....................................................................
	  if (line_string.find("$BOREHOLE_GIVEN_VALUE") != std::string::npos) // JOD 2022-02-10
	  {
		in.str(readNonBlankLineFromInputStream(*st_file));
		in >> borehole_data.value;
		connected_geometry_couplingType = 2; // for given source / sink term via RHS
		in.clear();

		std::cout << "\tBorehole - No coupling, but given primary value: " << borehole_data.value << '\n';
		continue;
	  }
	  //....................................................................
	  if (line_string.find("$BOREHOLE_MODIFIED_AQUIFER_PARAMETER") != std::string::npos) // JOD 2022-06-17
	  {
		in.str(readNonBlankLineFromInputStream(*st_file));
		in >> borehole_modified_aquifer_parameter;
		in.clear();

		std::cout << "\tBorehole - Modified aquifer parameter: " << borehole_modified_aquifer_parameter << '\n';
		continue;
	  }
	  //....................................................................
	  if (line_string.find("$AVERAGE_MODE") != std::string::npos) //JOD-2021-11-12
	  {
		  in.str(readNonBlankLineFromInputStream(*st_file));
		  in >> average_mode >> average_verbosity;  //  0: equal, 1: node area, 2: LIQUID_FLOW ST
		  std::cout << "\tAverage mode " << average_mode;
		  if(average_mode == 0)
			  std::cout << " - equal\n";
		  else if(average_mode == 1)
			  std::cout << " - With node area\n";
		  else if(average_mode == 2)
			  std::cout << " - With LIQUID_FLOW sink term (take care that values are kept)";
		  else
			  throw std::runtime_error("Average mode not supported");
		  in.clear();
		  continue;
	  }
	  //....................................................................
	  if (line_string.find("$VERBOSITY") != std::string::npos) //JOD-2021-11-12
	  {
		  in.str(readNonBlankLineFromInputStream(*st_file));
		  in >> verbosity;
		  in.clear();
		  continue;
	  }
	  //....................................................................

   } // end!new_keyword
   return position;
}


/**************************************************************************
 FEMLib-Method:
 Task: for CSourceTerm::Read
 Programing:
 11/2007 JOD Implementation
 02/2009 WW  Add a functionality to directly assign source terms to element nodes.
 **************************************************************************/
void CSourceTerm::ReadDistributionType(std::ifstream *st_file)
{
   std::stringstream in;
   // 03.2010 WW
   std::string aline;
   std::stringstream ss;
   int abuff, nLBC = 0;
   double bbuff;

   std::string dis_type_name;
   in.str(GetLineFromFile1(st_file));
   in >> dis_type_name;

   this->setProcessDistributionType (FiniteElement::convertDisType(dis_type_name));

   if (dis_type_name.compare (convertDisTypeToString (this->getProcessDistributionType())) != 0)
   {
      std::cerr << "Error in CSourceTerm::ReadDistributionType (): dist_type_name #" << dis_type_name << "#, new: " << convertDisTypeToString (this->getProcessDistributionType()) << "\n";
      exit (1);
   }

   if (   this->getProcessDistributionType() == FiniteElement::CONSTANT
       || this->getProcessDistributionType() == FiniteElement::CONSTANT_NEUMANN
       || this->getProcessDistributionType() == FiniteElement::CONSTANT_GEO      )
   {
      in >> geo_node_value;
      in.clear();
      if (geo_node_value==.0)
          std::cout << "Warning in CSourceTerm::ReadDistributionType (): Zero is set to the distribution type " << dis_type_name << ", which actually does nothing. There might be a mistake in yoru ST file. \n";
   }

   //outflux to surrounding depending on current value of solution at boundary (e.g. heat transfer to surrounding dependent on current wall temperature)
   if (this->getProcessDistributionType() == FiniteElement::TRANSFER_SURROUNDING){
	   in >> transfer_coefficient;
	   in >> value_surrounding;
	   in.clear();
   }

   //	if (dis_type_name.find("ANALYTICAL") != std::string::npos) {
   if (this->getProcessDistributionType() == FiniteElement::ANALYTICAL)
   {
      in >> analytical_material_group;            //Which material group is it being applied to
      in >> analytical_diffusion;                 //D value
      in >> analytical_porosity;                  //n value of matrix
      in >> analytical_tortousity;                //t value of matrix
      in >> analytical_linear_sorption_Kd;        //Linear sorption coefficient
      in >> analytical_matrix_density;            //Density of solid
      in >> number_of_terms;                      //no timesteps to consider in solution
      in >> resolution;                           //every nth term will be considered
      in >> factor;                               //to convert temperature to energy
      analytical = true;
      analytical_processes.push_back(convertPrimaryVariableToString (getProcessPrimaryVariable()));
      //		if (geo_type_name.compare("POLYLINE") == 0)
      if (this->getGeoType() == GEOLIB::POLYLINE)
         analytical_processes_polylines.push_back(geo_name);
      in.clear();
   }

	// If a linear function is given. 25.08.2011. WW  
	if (getProcessDistributionType() == FiniteElement::FUNCTION) 
	{ 
	  in.clear(); 
	  dis_linear_f = new LinearFunctionData(*st_file); 
	} 

   if (this->getProcessDistributionType() == FiniteElement::LINEAR || this->getProcessDistributionType() == FiniteElement::LINEAR_NEUMANN)
   {
      in >> nLBC;
      in.clear();
      for (int i = 0; i < nLBC; i++)
      {
         in.str(GetLineFromFile1(st_file));
         in >> abuff >> bbuff;
         in.clear();
         PointsHaveDistribedBC.push_back(abuff);
         DistribedBC.push_back(bbuff);
      }

      //      Read LINENODES AND VALUES......
      in.clear();
   }

   if (this->getProcessDistributionType() == FiniteElement::CRITICALDEPTH)
   {
      //KR critical_depth = true;
      in >> geo_node_value;
      in.clear();
      in.str(GetLineFromFile1(st_file));
      in >> st_rill_height;
      in.clear();
      //		dis_type = 6;
   }

   if (this->getProcessDistributionType() == FiniteElement::NORMALDEPTH)
   {
      dis_type_name = "NORMALDEPTH";
      in >> geo_node_value;
      in.clear();
      in.str(GetLineFromFile1(st_file));
      in >> normaldepth_slope >> st_rill_height;
      in.clear();
   }


   if (this->getProcessDistributionType() == FiniteElement::GREEN_AMPT)
   {
      dis_type_name = "GREEN_AMPT";
      in >> geo_node_value;
      in.clear();
      in.str(GetLineFromFile1(st_file));
      in >> sorptivity >> constant >> rainfall >> moistureDeficit;
      in.clear();
   }
   // Soure terms are assign to element nodes directly. 23.02.2009. WW
   if (this->getProcessDistributionType() == FiniteElement::DIRECT)
   {
      dis_type_name = "DIRECT";
      in >> fname;
      fname = FilePath+fname;
      in.clear();
   }
   if (this->getProcessDistributionType() == FiniteElement::RECHARGE_DIRECT)
   {
      dis_type_name = "RECHARGE_DIRECT";
      in >> fname;
      fname = FilePath+fname;
      in.clear();
   }

   // Soure terms from precipitation are assign to element nodes directly.03.2010. WW
   if(dis_type_name.find("PRECIPITATION")!=std::string::npos)
   {
      dis_type_name = "PRECIPITATION";
      in >> fname;
      fname = FilePath+fname;
      std::ifstream ins(fname.c_str());
      if(!ins.good())
      {
         std::cout<<"Could not find file "<<fname<<"\n";
         exit(0);
      }
      double timess;
      GIS_shape_head = new double[6];             // 07.06.2010. WW
      for(int i=0; i<6; i++)
      {
         getline(ins, aline);
         ss.str(aline);
         ss>> aline >> GIS_shape_head[i];
         ss.clear();

      }
      while(!ins.eof())
      {
         getline(ins, aline);
         if(aline.find("#STOP")!=std::string::npos)
            break;
         ss.str(aline);
         ss>> timess >> aline;
         precip_times.push_back(timess);
         precip_files.push_back(aline);

         ss.clear();
      }
      in.clear();
   }
if (this->getProcessDistributionType() == FiniteElement::CLIMATE)
   {
      dis_type_name = "CLIMATE";
      in >> fname; // base filename for climate input
      in.clear();

	  std::vector<GEOLIB::Point*> *stations (FileIO::RapidXMLInterface::readStationFile(FilePath + fname));

	  const size_t nStations(stations->size());
	  for (size_t i=0; i<nStations; i++)
			_weather_stations.push_back(static_cast<GEOLIB::Station*>((*stations)[i]));

	  delete stations;
   }


   if (this->getProcessDistributionType() == FiniteElement::RECHARGE)
   {
	   dis_type_name = "RECHARGE";
	   in >> geo_node_value;
	   in.clear();
   }


}


/**************************************************************************
 FEMLib-Method:
 Task: for CSourceTerm::Read
 Programing:
 11/2007 JOD Implementation
 06/2010 TF modification of the signature, added geo_obj and unique_name
 **************************************************************************/
void CSourceTerm::ReadGeoType(std::ifstream *st_file,
const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name)
{
   FileIO::GeoIO::readGeoInfo(this, *st_file, geo_name, geo_obj, unique_name);

   if (getProcessPrimaryVariable() == FiniteElement::EXCAVATION) //WW
   {
      std::stringstream strstr;
      strstr.str(GetLineFromFile1(st_file));
      //size_t tmp_geo_type;
      std::string sub_string;
      strstr >> sub_string >> _sub_dom_idx;
      strstr.clear();
   }
}


/**************************************************************************
 FEMLib-Method:
 Task: ST read function
 Programing:
 01/2004 OK Implementation
 06/2010 TF modification of the signature, added geo_obj and unique_name
 **************************************************************************/
bool STRead(const std::string &file_base_name,
const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name)
{
   char line[MAX_ZEILE];
   std::string line_string, st_file_name;
   std::ios::pos_type position;

   // File handling
   st_file_name = file_base_name + ST_FILE_EXTENSION;
   std::ifstream st_file(st_file_name.data(), std::ios::in);

   if (!st_file.good())
   {
      std::cout << "! Warning in STRead: No source terms !" << "\n";
      return false;
   }

   // Keyword loop
   std::cout << "STRead ... " << std::flush;
   while (!st_file.eof())
   {
      st_file.getline(line, MAX_ZEILE);
      line_string = line;
                                                  //Code included to make dynamic memory for analytical solution
      if (line_string.find("#STOP") != std::string::npos)
      {

         size_t no_source_terms(st_vector.size());
         size_t no_an_sol = 0, number_of_terms = 0;
                                                  //Need to find the memory size limits for the anal. solution.
         for (size_t i = 0; i < no_source_terms; i++)
         {
            if (st_vector[i]->isAnalytical())
            {
               no_an_sol++;
               number_of_terms = std::max(st_vector[i]->getNumberOfTerms(), number_of_terms);
            }
         }
         if (no_an_sol > 0)
         {
            for (size_t i = 0; i < no_source_terms; i++)
            {
               st_vector[i]->setNumberOfAnalyticalSolutions (no_an_sol);
               st_vector[i]->setMaxNumberOfTerms (number_of_terms);
            }
         }

         std::cout << "done, read " << st_vector.size() << " source terms" << "\n";

         return true;
      }
      //----------------------------------------------------------------------
                                                  // keyword found
      if (line_string.find("#SOURCE_TERM") != std::string::npos)
      {
         CSourceTerm *st(new CSourceTerm());
         std::ios::pos_type pos (st_file.tellg());
         position = st->Read(&st_file, geo_obj, unique_name);
         if (pos != position)
         {
        	if (st->getProcessPrimaryVariable()==FiniteElement::DISPLACEMENT_N)
        	{
        		st->SetPressureBoundaryConditionFlag(true);
        		CSourceTerm *st2(new CSourceTerm (*st)); // copies st's properties onto st2
        		st2->setProcessPrimaryVariable(FiniteElement::DISPLACEMENT_X);
        		st_vector.push_back(st2);
        		st->setProcessPrimaryVariable(FiniteElement::DISPLACEMENT_Y);
        	}
            st_vector.push_back(st);

         }
         /*else   removed by  JOD 2015-11-18
         {
            std::cerr << "WARNING: in STRead: could not read source term" << "\n";
            delete st;
         }*/
         st_file.seekg(position, std::ios::beg);
      }                                           // keyword found
   }                                              // eof

   std::cout << "done, read " << st_vector.size() << " source terms" << "\n";

   return true;
}


/**************************************************************************
 FEMLib-Method:
 Task: ST to right-hand-side vector
 Programing:
 01/2004 OK Implementation
 **************************************************************************/
void ST2RHS(std::string pcs_function, double* rhs_vector)
{
   int j;
   long i;
   long node_number_vector_length;
   CSourceTerm *m_st = NULL;
   long no_st = (long) st_vector.size();
   for (j = 0; j < no_st; j++)
   {
      m_st = st_vector[j];
      if (pcs_function.compare(convertPrimaryVariableToString (m_st->getProcessPrimaryVariable())) == 0)
      {
         node_number_vector_length = (long) m_st->node_number_vector.size();
         for (i = 0; i < node_number_vector_length; i++)
         {
            rhs_vector[m_st->node_number_vector[i]]
               = m_st->node_value_vector[i];
         }
      }
   }
}


/**************************************************************************
 FEMLib-Method:
 Task: ST to mesh nodes
 Programing:
 01/2004 OK Implementation
 **************************************************************************/
//void ST2NOD()
//{
//	CGLVolume *m_volume = NULL;
//	double st_value;
//
//	CGLPolyline *ply (NULL);
//	size_t points_along_polyline;
//
//	// Nodes
//	size_t no_st (st_vector.size());
//	for (size_t j = 0; j < no_st; j++) {
//		CSourceTerm *m_st = st_vector[j];
//		switch (m_st->getGeoType()) {
//		case GS_POINT:
//			m_st->node_number_vector.push_back(m_st->getGeoObjIdx());
//			break;
//		case GS_POLYLINE:
//			CGLPolyline *ply = GEOGetPLYByName(m_st->geo_prop_name);
//			points_along_polyline = ply->point_vector.size();
//			for (size_t i = 0; i < points_along_polyline; i++) {
//				m_st->node_number_vector.push_back(ply->point_vector[i]->id);
//			}
//			break;
//		case GS_SURFACE:
//			Surface *m_surface (GEOGetSFCByName(m_st->geo_prop_name));
//			long points_in_surface;
//			long* nodes (NULL);
//			nodes = GetPointsIn(m_surface, &points_in_surface);
//			//MB patch areas
//			for (size_t i = 0; i < static_cast<size_t>(points_in_surface); i++) {
//				m_st->node_number_vector.push_back(nodes[i]);
//			}
//			delete [] nodes;
//			break;
//		case GS_VOLUME:
//			m_volume = GEOGetVOL(m_st->geo_prop_name);
//			break;
//		default:
//			break;
//		} // switch
//	} // while
//	//========================================================================
//	size_t st_point_number;
//	size_t geo_point_number;
//	// Values
//	for (size_t j = 0; j < no_st; j++) {
//		CSourceTerm *m_st = st_vector[j];
//		switch (m_st->dis_type) {
//		case CONSTANT:
//			st_point_number = m_st->node_number_vector.size();
//			for (size_t i = 0; i < st_point_number; i++) {
//				m_st->node_value_vector.push_back(m_st->dis_prop[0]);
//			}
//			break;
//		case LINEAR: // for polylines
//			CGLPolyline *ply = GEOGetPLYByName(m_st->geo_prop_name);//CC
//
//			geo_point_number = ply->point_vector.size();
//			//if(!(st_point_number==geo_point_number)) Warning !
//			for (size_t i = 0; i < geo_point_number; i++) {
//				st_value = ply->point_vector[i]->propert; // i.e. node property
//				m_st->node_value_vector.push_back(st_value);
//			}
//			break;
//		} // switch
//	} // while
//}

/**************************************************************************
 FEMLib-Method: STWrite
 Task: master write function
 Programing:
 04/2004 OK Implementation
 last modification:
 05/2010 TF
 **************************************************************************/
void STWrite(std::string base_file_name)
{
   // File handling
   std::string st_file_name = base_file_name + ST_FILE_EXTENSION;
   std::fstream st_file(st_file_name.data(), std::ios::trunc | std::ios::out);
   st_file.setf(std::ios::scientific, std::ios::floatfield);
   st_file.precision(12);
   if (!st_file.good())
      return;

   st_file
      << "GeoSys-ST: Source Terms ------------------------------------------------\n";

   size_t no_st(st_vector.size());
   for (size_t j = 0; j < no_st; j++)
   {
      st_vector[j]->Write(&st_file);
   }
   st_file << "#STOP";
   st_file.close();
}


/**************************************************************************
 FEMLib-Method:
 Task: write function
 Programing:
 02/2004 OK Implementation
 04/2005 OK PRIMARY_VARIABLE
 06/2005 OK RIVER
 last modification:
 **************************************************************************/
void CSourceTerm::Write(std::fstream* st_file)
{
   //KEYWORD
   *st_file << "#SOURCE_TERM" << "\n";
   //--------------------------------------------------------------------
   //NAME+NUMBER
   *st_file << " $PCS_TYPE" << "\n";
   *st_file << "  ";
   //	*st_file << pcs_type_name << "\n";
   *st_file << convertProcessTypeToString (getProcessType()) << "\n";
   *st_file << " $PRIMARY_VARIABLE" << "\n";
   *st_file << "  ";
   *st_file << convertPrimaryVariableToString (getProcessPrimaryVariable()) << "\n";
   //--------------------------------------------------------------------
   //GEO_TYPE
   if (this->getProcessDistributionType() != FiniteElement::DIRECT
		   || this->getProcessDistributionType() != FiniteElement::RECHARGE_DIRECT)
   {
	   *st_file << " $GEO_TYPE" << "\n";
	   *st_file << "  ";
	   *st_file << getGeoTypeAsString() << " " << geo_name << "\n";
   }
   //--------------------------------------------------------------------
   // TIM_TYPE
   if (tim_type_name.size() > 0)                  //OK
   {
      *st_file << " $TIM_TYPE" << "\n";
      *st_file << "  ";
      *st_file << tim_type_name << "\n";
   }
   //--------------------------------------------------------------------
   //DIS_TYPE
   *st_file << " $DIS_TYPE" << "\n";
   *st_file << "  ";
   *st_file << convertDisTypeToString(this->getProcessDistributionType());
   switch (this->getProcessDistributionType())
   {
      case FiniteElement::CONSTANT:
         *st_file << " " << geo_node_value;
         *st_file << "\n";
         break;
      case  FiniteElement::CONSTANT_NEUMANN:
         *st_file << " " << geo_node_value;
         *st_file << "\n";
         break;
      case FiniteElement::LINEAR:
         *st_file << " " << (int) PointsHaveDistribedBC.size() << "\n";
         for (long i = 0; i < (long) PointsHaveDistribedBC.size(); i++)
         {
            *st_file << "  " << PointsHaveDistribedBC[i] << " ";
            *st_file << "  " << DistribedBC[i] << "\n";
         }
         break;
      case  FiniteElement::LINEAR_NEUMANN:
         *st_file << " " << PointsHaveDistribedBC.size() << "\n";
         for (size_t i = 0; i < PointsHaveDistribedBC.size(); i++)
         {
            *st_file << "  " << PointsHaveDistribedBC[i] << " ";
            *st_file << "  " << DistribedBC[i] << "\n";
         }
         break;
      case  FiniteElement::DIRECT:
      case  FiniteElement::RECHARGE_DIRECT:
    	  *st_file << " " << this->fname << "\n";
    	  break;
      default:
         std::cerr << "this distributition type is not handled in CSourceTerm::Write" << "\n";
   }
   //--------------------------------------------------------------------
}

/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 04/2004 OK Implementation
 11/2007 JOD Reaktivation
 last modification:
 **************************************************************************/
//void CSourceTerm::SetDISType()
//{
////	if (this->getProcessDistributionType() == CONSTANT)
////		dis_type = 1;
////////	if (dis_type_name.compare("CONSTANT_GEO") == 0)
////////		dis_type = 12; //SB flux is distributed along polyline. To do. 4.10.06
////	if (this->getProcessDistributionType() == LINEAR)
////		dis_type = 2;
////	if (this->getProcessDistributionType() == CONSTANT_NEUMANN)
////		dis_type = 3;
////	if (this->getProcessDistributionType() == LINEAR_NEUMANN)
////		dis_type = 4;
////	if (this->getProcessDistributionType() == RIVER) {
////		dis_type = 5;
////	}
//
////	if (this->getProcessDistributionType() == CRITICALDEPTH)
////		dis_type = 6;
////	if (this->getProcessDistributionType() == SYSTEM_DEPENDENT)
////		dis_type = 7; //YD
////	if (this->getProcessDistributionType() == NORMALDEPTH)
////		dis_type = 8; //JOD MB
////	if (this->getProcessDistributionType() == ANALYTICAL)
////		dis_type = 9;//CMCD 02 2006
//////	if (dis_type_name.compare("PHILIP") == 0)
////////	if (this->getProcessDistributionType() == PHILIP)
//////		dis_type = 10; // JOD
////	if (this->getProcessDistributionType() == GREEN_AMPT)
////		dis_type = 11; // JOD
//	// if(dis_type_name.compare("CONSTANT")==0) dis_type = 0;
//	// if(dis_type_name.compare("LINEAR")  ==0) dis_type = 1;
//}

/**************************************************************************
 FEMLib-Method:
 Task: set ST group member
 Programing:
 02/2004 OK Implementation
 09/2004 WW Face integration of Neumann boundary condition for all element type
 09/2004 WW Interpolation for piece-wise linear distributed source term or BC
 03/2005 OK LINE sources
 02/2005 MB River condition, CriticalDepth
 08/2006 WW Re-implementing edge,face and domain integration versatile for all element types
 04/2006 OK MSH types
02/2009 WW Direct assign node source terms
**************************************************************************/
void CSourceTermGroup::Set(CRFProcess* m_pcs, const int ShiftInNodeVector,
		std::string this_pv_name)
{
  bool set = false; // CB

   if (this_pv_name.size() != 0)                  //WW
      pcs_pv_name = this_pv_name;
   m_msh = FEMGet(pcs_type_name);                 //m_pcs->m_msh; //

   if (m_msh)                                     //WW
   {
      /// In case of P_U coupling monolithic scheme
      if (m_pcs->type == 41)                      //WW Mono
      {
                                                  //Deform
         if (pcs_pv_name.find("DISPLACEMENT") != std::string::npos)
            m_pcs->m_msh->SwitchOnQuadraticNodes(true);
         else
            m_pcs->m_msh->SwitchOnQuadraticNodes(false);
      } else if (m_pcs->type == 4)
      m_pcs->m_msh->SwitchOnQuadraticNodes(true);
      else
         m_pcs->m_msh->SwitchOnQuadraticNodes(false);
      //====================================================================
      long no_st = (long) st_vector.size();
      for (long i = 0; i < no_st; i++)
      {
         CSourceTerm *source_term (st_vector[i]);
         if(m_pcs->getProcessType() != source_term->getProcessType() )
            continue;

         source_term->setSTVectorGroup(i);

         // 07.01.2011. WW
         if(source_term->getProcessDistributionType()==FiniteElement::PRECIPITATION)
            continue;

         if (source_term->isCoupled())
            m_msh_cond = FEMGet(source_term->pcs_type_name_cond);

         if (source_term->getProcessType() == FiniteElement::MASS_TRANSPORT)
             //if ( cp_vec[cp_name_2_idx[convertPrimaryVariableToString(source_term->getProcessPrimaryVariable())]]->getProcess() != m_pcs )
             if ( cp_vec[source_term->getProcessCompVecIndex()]->getProcess() != m_pcs ) //CB cannot match CONCENTRATION with comp name!!
                 continue;
          //-- 23.02.3009. WW
         if (source_term->getProcessDistributionType()==FiniteElement::DIRECT
        		 || source_term->getProcessDistributionType()==FiniteElement::CLIMATE)
		 { //NB For climate ST, the source terms (recharge in this case) will also be assigned directly to surface nodes
           source_term->DirectAssign(ShiftInNodeVector);
           continue;
         }

         if (source_term->getProcessDistributionType()==FiniteElement::RECHARGE_DIRECT)
         {
        	 source_term->DirectAssign(ShiftInNodeVector);
        	 set = true;
         }
		 
         //CB cannot match CONCENRATION with comp name!!
         // modified to fix bug with MASS_TRANSPORT
         if (convertProcessTypeToString (source_term->getProcessType ()).compare(pcs_type_name) == 0)
         {
           if(convertPrimaryVariableToString(source_term->getProcessPrimaryVariable()).compare(pcs_pv_name) == 0)
             set = true;
           else if ((source_term->getProcessType() == FiniteElement::MASS_TRANSPORT) 
                && (cp_vec[source_term->getProcessCompVecIndex()]->compname.compare(pcs_pv_name) == 0))
             set = true;
         }		 
		 

		// 05/2012 BG, it does not work for mass transport if the second part with "CONCENTRATION" is not added !! problem identified by GK und HS
         //if ((convertProcessTypeToString (source_term->getProcessType ()).compare(pcs_type_name) == 0)
         //   && ((convertPrimaryVariableToString(source_term->getProcessPrimaryVariable()).compare(pcs_pv_name) == 0) || (convertPrimaryVariableToString(source_term->getProcessPrimaryVariable()).compare("CONCENTRATION1") == 0)))
         // if ( source_term->getProcess() == m_pcs )
         //CB cannot match CONCENRATION with comp name!!
         // modified to fix bug with MASS_TRANSPORT
         if(set)
         {
             set=false;
             source_term->setProcess (m_pcs);      // HS: 01.09.2009
             if (source_term->getGeoType() == GEOLIB::POINT)
                 SetPNT(m_pcs, source_term, ShiftInNodeVector);
             if (source_term->getGeoType () == GEOLIB::POLYLINE) {
                 SetPLY(source_term, ShiftInNodeVector);
             }
             if (source_term->getGeoType () == GEOLIB::SURFACE)
                 SetSFC(source_term, ShiftInNodeVector);
             if (source_term->getGeoType () == GEOLIB::GEODOMAIN)
                 SetDMN(source_term, ShiftInNodeVector);
             if (source_term->fct_name.size() > 0)
                 fct_name = source_term->fct_name;
			 // Recovery this functionality. 12.08.2011 WW
			// MSH types //OK4310
			if(source_term->msh_type_name.compare("NODE")==0)
				source_term->SetNOD();


			if (source_term->getProcessDistributionType()==FiniteElement::RECHARGE
					|| source_term->getProcessDistributionType()==FiniteElement::RECHARGE_DIRECT)	//MW
				MshEditor::sortNodesLexicographically(m_pcs->m_msh);

         }                                        // end pcs name & pv
      }                                           // end st loop
   }                                              // end msh
   else
      std::cout << "Warning in CSourceTermGroup::Set - no MSH data" << "\n";

   // CB MERGE // 
   /*  for (long i = 0; i < st_vector.size(); i++){
     if (/---/0==0){ // Adapt to connected geometries
       WriteNodeConnections();
       break;  // only once, loop over st inside
     }
   }*/

}

// CB
void CSourceTermGroup::WriteNodeConnections()
{
  return;
  
  std::string m_file_name = FileName + "_" + "_ST_Node_Connections.asc";
  std::ofstream os;
  CSourceTerm *source_term;

  std::vector < string > surfaces;

  for (long i = 0; i < st_vector.size(); i++)
  {
    bool found = false;
    source_term=st_vector[i];
      
    // check for this st
    if (/*source_term->connected_gemoetry*/0 == 0){
      int j;
      for (j = 0; j < surfaces.size(); j++){
        if ((surfaces[j].compare(source_term->geo_name) != 0))
        {
          // not found yet
        }
        else
          found = true; // data already written, skip writung
      }
      if (!found){
        // write data
        if (j == 0)
          os.open(m_file_name.c_str(), ios::trunc | ios::out);
        else
          os.open(m_file_name.c_str(), ios::app | ios::out);
        if (!os.good())
        {
          cout << "Failure to open file: " << m_file_name << "\n";
          abort();
        }
        else{  // adapt to conneced geometries
          for (int k = 0; k < source_term->st_node_ids.size(); k++){
            os << source_term->st_node_ids[k] << " " << source_term->st_node_ids[k] << "\n";
          }
          os.close();
        }
        // store geo name to avoid duplicate output
        surfaces.push_back(source_term->geo_name);
      }
    }
  }
  //clean up
  surfaces.clear();
}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 01/2005 OK Implementation
 last modified:
 **************************************************************************/
void STDelete()
{
   long i;
   int no_st = (int) st_vector.size();
   for (i = 0; i < no_st; i++)
   {
      delete st_vector[i];
   }
   st_vector.clear();
}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 01/2005 OK Implementation
 last modified:
 **************************************************************************/
void STGroupsDelete()
{
   CSourceTermGroup* m_st_group = NULL;
   std::list<CSourceTermGroup*>::const_iterator p = st_group_list.begin();
   while (p != st_group_list.end())
   {
      m_st_group = *p;
      delete m_st_group;
      //st_group_list.remove(*p);
      ++p;
   }
   st_group_list.clear();
}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 03/2005 OK Implementation
 last modification:
 **************************************************************************/
// 05/2010 TF not needed for compilation
//void STCreateFromLIN(vector<CGLLine*>lin_properties_vector)
//{
//  long i;
//  CSourceTerm *m_st = NULL;
//  CGLLine *m_lin = NULL;
//  long lin_properties_vector_size = (long)lin_properties_vector.size();
//  for(i=0;i<lin_properties_vector_size;i++){
//    m_st = new CSourceTerm();
//    m_lin = lin_properties_vector[i];
//    m_st->pcs_pv_name = "PRESSURE1"; // ToDo
//    m_st->geo_type_name = "LINE";
//    m_st->setGeoName (m_lin->name);
//    m_st->geo_id = m_lin->gli_line_id;
//    m_st->dis_type_name = "CONSTANT_NEUMANN";
//    m_st->geo_node_value = m_lin->value;
//    m_st->tim_type_name = m_lin->name;
//    st_vector.push_back(m_st);
//  }
//}

/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 03/2005 OK Implementation
 last modification:
 05/2010 TF restructured a little bit
 **************************************************************************/
CSourceTerm* STGet(std::string geo_name)
{
   size_t st_vector_size(st_vector.size());
   for (size_t i = 0; i < st_vector_size; i++)
   {
      if (st_vector[i]->getGeoName().compare(geo_name) == 0)
         return st_vector[i];
   }
   return NULL;
}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 04/2005 OK Implementation
 last modification:
 **************************************************************************/
CSourceTermGroup* STGetGroup(std::string pcs_type_name, std::string pcs_pv_name)
{
   CSourceTermGroup *m_st_group = NULL;
   std::list<CSourceTermGroup*>::const_iterator p_st_group = st_group_list.begin();
   while (p_st_group != st_group_list.end())
   {
      m_st_group = *p_st_group;
      if (m_st_group->pcs_type_name.compare(pcs_type_name) == 0
         && (m_st_group->pcs_pv_name.compare(pcs_pv_name) == 0))
         return m_st_group;
      ++p_st_group;
   }
   return NULL;
}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 02/2006 MB Implementation
 **************************************************************************/
                                                  //WW
double GetConditionalNODValue(CSourceTerm* st, CNodeValue* cnodev)
{
   int nidx;
   double value_cond = 0.0;
   double NodeReachLength;
   CRFProcess* pcs_cond (PCSGet(st->pcs_type_name_cond));
   long node_cond;

   //WW  node_cond = group_vector[i]->msh_node_number_conditional;
                                                  //WW
   node_cond = cnodev->msh_node_number_conditional;
   nidx = pcs_cond->GetNodeValueIndex(st->pcs_pv_name_cond) + 1;
   value_cond = pcs_cond->GetNodeValue(node_cond, nidx);

   if (st->pcs_pv_name_cond.find("FLUX") != std::string::npos)
   {
      //WW    NodeReachLength = group_vector[i]->node_value;
      NodeReachLength = cnodev->node_value;       //WW
      value_cond = value_cond * NodeReachLength;
   }

   return value_cond;
}


/**************************************************************************
 FEMLib-Method:
 Task: Calculate Philips (1957) two-term infiltration flux term
 calculated separately for each node
 Programing:
 05/2007 JOD Implementation
 09/2010 KR cleaned up code
 **************************************************************************/
void GetPhilipNODValue(double &value, const CSourceTerm* m_st)
{
   double infiltration = m_st->constant + m_st->sorptivity / sqrt(aktuelle_zeit);
   infiltration = std::min(m_st->rainfall, infiltration);
   value = infiltration * (-value);
}


/**************************************************************************
 FEMLib-Method:
 Task: Calculate Green-Ampt infiltration flux term
 for homogeneous soil, includes water depth
 writes cumulative infiltration in COUPLINGFLUX
 solution is sensitive to time step
 infiltration is estimated with first order derivative
 of cumulative infiltration
 Programing:
 05/2007 JOD Implementation

**************************************************************************/
void GetGreenAmptNODValue(double &value, CSourceTerm* m_st, long msh_node)
{

   double F, Fiter, Fold, infiltration;
   double conductivity, suction, Theta, wdepth;
   double a, b;
   CFEMesh* m_msh = NULL;
   CRFProcess* m_pcs_this = NULL;
   m_pcs_this = PCSGet(convertProcessTypeToString (m_st->getProcessType()));
   m_msh = m_pcs_this->m_msh;

   double area = value;

   wdepth = std::max(0., m_pcs_this->GetNodeValue(msh_node,
      m_pcs_this->GetNodeValueIndex("HEAD") + 1)
      - m_msh->nod_vector[msh_node]->getData()[2]);
   conductivity = m_st->constant;
   suction = m_st->sorptivity + wdepth;           // water depth included
   Theta = m_st->moistureDeficit * suction;

   Fold = m_pcs_this->GetNodeValue(msh_node, m_pcs_this->GetNodeValueIndex(
      "COUPLING"));
   F = Fold;

   do                                             // Newton iteration loop
   {
      Fiter = F;
      if (Fiter == 0)                             // avoids a = 0
         Fiter = 1.e-5;
      a = 1 - Theta / (Fiter + Theta);

      b = Fold - Fiter + Theta * log((Fiter + Theta) / (Fold + Theta))
         + conductivity * dt;                     // dt = timeStep
      F = Fiter + b / a;
   } while (fabs(F - Fiter) > 1.e-10);

   infiltration = (F - Fold) / dt;

   if (infiltration > m_st->rainfall)             // + wdepth / timestep )   // compare with available water
      infiltration = m_st->rainfall;              //  +  wdepth / timestep ;

   F = infiltration * dt + Fold;

   m_pcs_this->SetNodeValue(msh_node,
                                                  // output is cumulative infiltration
      m_pcs_this->GetNodeValueIndex("COUPLING") + 1, F);

   infiltration *= -area;
   value = infiltration;

}


/**************************************************************************
 FEMLib-Method:
 Task: Calculate coupling flux
 Mutual coupling of overland, Richards and groundwater flow
 Programing:
 01/2007 JOD Implementation
 09/2010 KR cleaned up code
 **************************************************************************/
#if !defined(USE_PETSC)
//&& !defined(NEW_EQS) // && defined(other parallel libs)//03~04.3012. WW
//#ifndef NEW_EQS                                   //WW. 06.11.2008
void GetCouplingNODValue(double &value, CSourceTerm* st, CNodeValue* cnodev)
{
	if((st->getProcessType() == FiniteElement::MASS_TRANSPORT ||
		      st->getProcessType() == FiniteElement::HEAT_TRANSPORT)
	 && st->pcs_type_name_cond == "LIQUID_FLOW")
	{
		GetCouplingNODValueConvectiveForm(value, st, cnodev->msh_node_number);
	}
	else if (st->getProcessType() == FiniteElement::GROUNDWATER_FLOW ||
      st->getProcessType() == FiniteElement::RICHARDS_FLOW ||
	  st->getProcessType() == FiniteElement::LIQUID_FLOW ||
      st->getProcessType() == FiniteElement::MULTI_PHASE_FLOW ||
      st->getProcessType() == FiniteElement::MASS_TRANSPORT ||
      st->getProcessType() == FiniteElement::HEAT_TRANSPORT)
	{
      GetCouplingNODValuePicard(value, st, cnodev);
	}
	else if (st->getProcessType() == FiniteElement::OVERLAND_FLOW)
	  GetCouplingNODValueNewton(value, st, cnodev);
	else
      std::cout << "Error in GetCouplingNODValue";
}
#endif

/**************************************************************************
 FEMLib-Method:
 Task: Calculate coupling flux for GetCouplingNODValue
 for the soil or groundwater flow case
 prerequisite: fluid data in mfp_vector[0] (2 times)
 Programing:
 01/2007 JOD Implementation
 10/2008 JOD overland node shifting for soil columns, averaging automatically 4.7.10
 12/2012 JOD Extension to TWO_PHASE_FLOW  5.3.07
 **************************************************************************/
#if !defined(USE_PETSC)
//&& !defined(NEW_EQS) // && defined(other parallel libs)//03~04.3012. WW
//#ifndef NEW_EQS                                   //WW. 06.11.2008
void GetCouplingNODValuePicard(double &value
	, CSourceTerm* m_st,
CNodeValue* cnodev)
{

   double relPerm, condArea;
   double h_this, h_cond, z_this, z_cond, h_cond_shifted, help;
   double gamma = 1;
   CRFProcess* m_pcs_this = NULL;
   CRFProcess* m_pcs_cond = NULL;  
   m_pcs_this = PCSGet(convertProcessTypeToString(m_st->getProcessType()));
   m_pcs_cond = PCSGet(m_st->pcs_type_name_cond);
   long mesh_node_number, mesh_node_number_conditional;

   if( (mesh_node_number = cnodev->msh_node_number) <  (long)m_pcs_this->m_msh->nod_vector.size())      
   {     // liquid 
	   mesh_node_number_conditional = cnodev->msh_node_number_conditional; 
	   if (m_st->getProcessType() == FiniteElement::RICHARDS_FLOW || m_st->getProcessType() == FiniteElement::LIQUID_FLOW)
          gamma = mfp_vector[0]->Density() * GRAVITY_CONSTANT;  // liquid pressure as a primary variable
       else if (m_st->getProcessType() == FiniteElement::MULTI_PHASE_FLOW)
	      gamma = -mfp_vector[0]->Density() * GRAVITY_CONSTANT; // capillary pressure as a primary variable
   }
   else
   {       // gas
	   mesh_node_number -= m_pcs_this->m_msh->nod_vector.size(); 
	   /*if( m_st->pcs_type_name_cond == "MULTI_PHASE_FLOW") 
	      mesh_node_number_conditional = cnodev->msh_node_number_conditional - m_pcs_this->m_msh->nod_vector.size(); 
	   else*/
	   mesh_node_number_conditional = cnodev->msh_node_number_conditional;
       gamma = mfp_vector[1]->Density() * GRAVITY_CONSTANT; // gas pressure as a primary variable
   }
   
   GetCouplingFieldVariables(m_pcs_this, m_pcs_cond, &h_this, &h_cond, &h_cond_shifted, &help,
      &z_this, &z_cond, m_st, cnodev, mesh_node_number, mesh_node_number_conditional, gamma);   // z_cond shifted for soil columns

   /////   relative interface permeability
   if(m_st->pcs_pv_name_cond == "PRESSURE2")// || m_st->pcs_type_name_cond == "OVERLAND_FLOW")
   { // gas
	    CRFProcess* m_pcs_overland = NULL;
		double head_overland;
        m_pcs_overland = PCSGet(m_st->pcs_type_name_cond);
        head_overland = m_pcs_overland->GetNodeValue(cnodev->msh_node_number_conditional, m_pcs_overland->GetNodeValueIndex("HEAD") + 1);
	    
		relPerm = m_st->GetRelativeInterfacePermeability(m_pcs_overland, head_overland, cnodev->msh_node_number_conditional);    // overland
	}
   else
   {  // liquid
     if (h_this < h_cond_shifted)  // flow direction from overland compartment
        relPerm = m_st->GetRelativeInterfacePermeability(m_pcs_cond, h_cond, cnodev->msh_node_number_conditional);
     else          // flow direction from subsurface compartment (, where now relPerm = 1)
	    relPerm = m_st->GetRelativeInterfacePermeability(m_pcs_this, h_this, cnodev->msh_node_number_conditional);
   }

		
    if( m_st->explicit_surface_water_pressure) { // put hyrostatic surface liquid pressure explicitly in coupling flux
		    h_cond = m_pcs_cond->GetNodeValue(mesh_node_number_conditional, 0);
	        h_cond_shifted =  h_cond - z_cond + z_this;
	}
  
   condArea = value * relPerm * m_st->getCoupLeakance();   // area * total interface permeability (m^2/s)        

   if (m_st->channel)                   // wetted perimeter, pcs_cond must be overland flow
      condArea *= m_st->channel_width + (h_cond - z_cond);
   ///// RHS part (a surface liquid pressure term) of coupling flux (m^3/s)                                            
   value = CalcCouplingValue(condArea, h_this, h_cond_shifted, z_cond, m_st);


   if (m_st->getProcessType() == FiniteElement::GROUNDWATER_FLOW && h_this < z_cond
      && m_st->pcs_type_name_cond == "OVERLAND_FLOW")
	 condArea = 0;   //  decoupled, only RHS
 

   if ( m_st->getProcessPrimaryVariable() == FiniteElement::PRESSURE 
		&& (mesh_node_number == cnodev->msh_node_number) ) // only for liquid phase
      condArea /= gamma; // pressure as a primary variable
   /////

#ifdef NEW_EQS
   // JOD 2018-5-17
   CSparseMatrix* A = NULL;
   A = m_pcs_this->get_eqs_new()->get_A();
   (*A)(cnodev->msh_node_number, cnodev->msh_node_number) += condArea;
#else
   MXInc(cnodev->msh_node_number, cnodev->msh_node_number, condArea);
#endif

   if(m_st->getProcessType() == FiniteElement::MULTI_PHASE_FLOW &&  m_st->pcs_pv_name_cond == "HEAD"   // 
		 && m_st->pcs_type_name_cond == "OVERLAND_FLOW" ) // ??? 
   {   // gas pressure term 

#ifdef NEW_EQS
	   // JOD 2018-5-17
	   (*A)(cnodev->msh_node_number, cnodev->msh_node_number+ m_pcs_this->m_msh->nod_vector.size()) -= condArea;
#else
		MXInc(cnodev->msh_node_number, cnodev->msh_node_number+ m_pcs_this->m_msh->nod_vector.size(), -condArea);
#endif
	   value  -= condArea * m_st->coup_given_value;
   }
  
   if (m_st->getProcessType() == FiniteElement::GROUNDWATER_FLOW)                   
   m_pcs_this->SetNodeValue( mesh_node_number,   // for GW water balance ply / pnt 
	   m_pcs_this->GetNodeValueIndex("FLUX") + 1, condArea); 


   if (m_st->no_surface_water_pressure)           // neglect hydrostatic surface liquid pressure
     value = 0; 

}
#endif


/**************************************************************************
 FEMLib-Method:
 Task: Source term for convective form of ADE
 Programing:
 03/2020 JOD Implementation

 Restrictions:
 	 takes mfp_vector[0]
 	 density model gets temperature and only that
 **************************************************************************/
#if !defined(USE_PETSC)
void GetCouplingNODValueConvectiveForm(double &value, CSourceTerm* m_st, const long mesh_node_number)//CNodeValue* cnodev)
{
    //const long mesh_node_number = cnodev->msh_node_number;
    long eq_index;
    long dom_node_index;
    
    //double  poro = 0.0;
    //int material_group;
    //long msh_ele;
    //size_t number_of_connected_elements;
    
    //get process
   CRFProcess* m_pcs_this = NULL;
   m_pcs_this = PCSGet(convertProcessTypeToString(m_st->getProcessType()));
   double nodal_val =  m_pcs_this->GetNodeValue(mesh_node_number, 1);

   eq_index = m_pcs_this->m_msh->nod_vector[mesh_node_number]->GetEquationIndex();
   //std::cout << temperature << ", ";
   //double density_test = mfp_vector[0]->Density();
   //std::cout << " mesh node number: " << mesh_node_number << "; Value: "<< value << '\n';
   //std::cout << " Equation Index: " << eq_index << '\n';

   if (mfp_vector[0]->get_flag_volumetric_heat_capacity())
       value *= mfp_vector[0]->get_volumetric_heat_capacity();
   else
       value *= mfp_vector[0]->Density(&nodal_val  // only first value of array set and as temperature
    		   ) * mfp_vector[0]->SpecificHeatCapacity(NULL, true);

#if defined(NEW_EQS)
#if defined(USE_MPI)  // JOD 2020-04-08
   //std::cout << "myrank: " << myrank << "; mesh node number: " << mesh_node_number << "; Value: " << value << '\n' ;
   CSparseMatrix* A = dom_vector[myrank]->get_eqs()->get_A();
   dom_node_index = dom_vector[myrank]->GetDOMNode(mesh_node_number);
   //std::cout << "mesh node number: " << mesh_node_number << "; dom_node_index: " << dom_node_index << '\n';
   (*A)(dom_node_index, dom_node_index) -= value;
#else // not USE_MPI
   // JOD 2020-3-25
   CSparseMatrix* A = m_pcs_this->get_eqs_new()->get_A();
   (*A)(mesh_node_number, mesh_node_number) -= value;
#endif
	
#else  // not  NEW_EQS
   MXInc(mesh_node_number, mesh_node_number, -value);
#endif

     value = 0;  // no right hand side term
}
#endif

/**************************************************************************
 FEMLib-Method:
 Task: Calculate coupling flux for GetCouplingNODValue
 for the overland flow case
 prerequisite: fluid data in mfp_vector[0]
 Programing:
 01/2007 JOD Implementation
 10/2008 JOD node shifting for soil columns 4.7.10
 12/2012 JOD coupling with TWO_PHASE_FLOW  5.3.07
 **************************************************************************/
#if !defined(USE_PETSC)
//&& !defined(NEW_EQS) // && defined(other parallel libs)//03~04.3012. WW
//#ifndef NEW_EQS                                   //WW. 06.11.2008
void GetCouplingNODValueNewton(double &value, CSourceTerm* m_st,
CNodeValue* cnodev)
{
   double relPerm, area, condArea, gamma = 0.0;
   double h_this, h_cond, z_this, z_cond, h_this_shifted, h_this_averaged;
   double epsilon = 1.e-7, value_jacobi, h_this_epsilon = 0.0,
      relPerm_epsilon, condArea_epsilon;          //OK411 epsilon as in pcs->assembleParabolicEquationNewton
   
   CRFProcess* m_pcs_cond = NULL;
   CRFProcess* m_pcs_this = NULL;
   m_pcs_this = PCSGet(convertProcessTypeToString(m_st->getProcessType()));
   m_pcs_cond = PCSGet(m_st->pcs_type_name_cond);
   ///// 
   if(m_st->pcs_type_name_cond == "RICHARDS_FLOW" || m_st->pcs_type_name_cond == "LIQUID_FLOW")
     gamma = mfp_vector[0]->Density() * GRAVITY_CONSTANT; // liquid pressure as a primary variable
   else if(m_st->pcs_type_name_cond == "MULTI_PHASE_FLOW")
     gamma = -mfp_vector[0]->Density() * GRAVITY_CONSTANT; // capillary pressure as a primary variable

   GetCouplingFieldVariables(m_pcs_this, m_pcs_cond, &h_this, &h_cond, &h_this_shifted,
                                // if soil column, z_this, h_this_averaged, h_this_shifted are shifted to the column 
      &h_this_averaged, &z_this, &z_cond, m_st, cnodev, cnodev->msh_node_number, cnodev->msh_node_number_conditional, gamma);

   /////   relative coupling interface permeability (calculate ALWAYS implicitly)
   if (h_this_shifted > h_cond)
   { // flow direction from the overland compartment
      relPerm = m_st->GetRelativeInterfacePermeability(m_pcs_this, h_this, cnodev->msh_node_number);
      relPerm_epsilon = m_st->GetRelativeInterfacePermeability(m_pcs_this, h_this + epsilon, // for jacobian
		  cnodev->msh_node_number);
   }
   else // flow direction from a subsurface compartment
      relPerm = relPerm_epsilon = 1; 

   if (m_st->node_averaging) 
   {  // h_this_epsilon for jacobian if multiple overland nodes are coupled to a soil column
      for (long i = 0; i < (long) cnodev->msh_node_numbers_averaging.size(); i++)
         if (cnodev->msh_node_numbers_averaging[i]
         == cnodev->msh_node_number)
            h_this_epsilon = h_this_averaged + epsilon
               * cnodev->msh_node_weights_averaging[i];
   } 
   else // h_this_epsilon for jacobian
      h_this_epsilon = h_this + epsilon;

   if (m_st->no_surface_water_pressure)  // neglect hydrostatic surface liquid pressure
      h_this_epsilon = z_this;

   if (m_st->explicit_surface_water_pressure)       
   {    // put hydrostatic surface liquid pressure explicitly in the coupling flux
	  h_this = m_pcs_this->GetNodeValue(cnodev->msh_node_number, 0);
      h_this_epsilon = h_this; // for jacobian

	   if (m_st->node_averaging)                  
      {               // shift overland node to soil column 
         h_this_shifted = h_this - m_pcs_this->m_msh->nod_vector[cnodev->msh_node_number]->getData()[2] + z_cond;
         h_this_averaged = 0;
         for (long i = 0; i
            < (long) cnodev->msh_node_numbers_averaging.size(); i++)
            h_this_averaged 
               += cnodev->msh_node_weights_averaging[i]
               * (m_pcs_this->GetNodeValue(
               cnodev->msh_node_numbers_averaging[i],
               0)
               - m_pcs_this->m_msh->nod_vector[cnodev->msh_node_numbers_averaging[i]]->getData()[2]);

         h_this_averaged += z_cond;
      
	   }

   } // end explicit
 
   area = value;                                 
   condArea = condArea_epsilon = area * m_st->getCoupLeakance();     
   condArea *= relPerm;
   condArea_epsilon *= relPerm_epsilon;

   if (m_st->channel)  
   {                           // wetted perimeter
      condArea *= m_st->channel_width + (h_this - z_this);
	  condArea_epsilon *= m_st->channel_width + (h_this_epsilon - z_this);
   }
   //////  coupling flux as a source term (m^3/s)   
   value = CalcCouplingValue(condArea, h_this_averaged, h_cond, z_cond, m_st);

   value_jacobi = -condArea_epsilon * (h_cond - h_this_epsilon) + value; // it's a Newton iteration
#ifndef NEW_EQS
   MXInc(cnodev->msh_node_number, cnodev->msh_node_number, value_jacobi 
      / epsilon);
   // else ...
#endif


   /////  output 
   m_pcs_this->SetNodeValue(cnodev->msh_node_number,
   // coupling flux (m/s)
      m_pcs_this->GetNodeValueIndex("COUPLING") + 1, -value/ area);
 
 
}
#endif

/**************************************************************************
 FEMLib-Method:
 Task: Calculate relative coupling permeability for GetCouplingNODValue(***).
 Programing: prerequisites: phase 0 in mfp
 06/2007 JOD Implementation
 10/2008 JOD include leakance 4.7.10
 10/2012 JOD extension to two-phase flow
 Gives 0 < relPerm < 1 dependent on liquid depth in the overland compartment
                       with reference to an interface thickness
 **************************************************************************/

double CSourceTerm::GetRelativeInterfacePermeability(CRFProcess* pcs, double head, long msh_node)
{

   double relPerm, sat, z;
  
   if(pcs->getProcessType() == FiniteElement::OVERLAND_FLOW)
   {                // flow direction from the overland compartment 

	 z = pcs->m_msh->nod_vector[msh_node]->getData()[2];
     sat =  (head - z) / std::max(1.e-6, st_rill_height); // liquid content in an interface
	       
     if( sat > 1 )
	    relPerm = 1; // liquid-covered surface (interface is filled)
     else if( sat <= 0 )
	    relPerm = 0;  // no liquid on the surface
     else
	    relPerm = pow(sat, 2*(1- sat)); // interface is partially filled by a liquid
		 
	 CRFProcess* m_pcs_multiPhase = NULL;  
     m_pcs_multiPhase = PCSGet("MULTI_PHASE_FLOW");
	
	 if(relPerm == 1 && pcs_pv_name_cond == "PRESSURE2") // surface closed for gas
	 {
        double capillaryPressure =  m_pcs_multiPhase->GetNodeValue(0, 3);
		
	    if( capillaryPressure > air_breaking_capillaryPressure || air_breaking == true ) 
		{	
	        relPerm *=  air_breaking_factor; // air-breaking
		 	air_breaking = true;
			
			if( capillaryPressure <  air_closing_capillaryPressure ) 
			{   // and  air-breaking (hysteresis)
			   air_breaking = false;
			}
		} // end air_breaking	
	 } // end relPerm = 1
	

	 if (pcs_pv_name_cond == "PRESSURE2")
	 {          // gas
		 relPerm = (1 - relPerm)* (1-  coup_residualPerm) +  coup_residualPerm;    
	 }
  } // end overland flow
  else
    relPerm = 1;  // flow direction not from the overland compartment 
	
  return relPerm;
}

/**************************************************************************
 FEMLib-Method:
 Task: Coupling of overland and soil flow by using water depth as soil boundary
 condition and flux term as overland source term according to
 Morita and Yen, J. Hydr. Eng. 184, 2002
 Programing: prerequisites: constant precipitation with assigned duration,
 phase = 0 in mfp, soil data in mmp_vetor[1] !!!!!
 06/2007 JOD Implementation
 **************************************************************************/
//#if !defined(NEW_EQS)
//&& !defined(USE_PETSC)                                   //WW. 06.11.2008
void GetCouplingNODValueMixed(double& value, CSourceTerm* m_st,
CNodeValue* cnodev)
{

   double cond1, cond0, pressure1, pressure0, bc_value, depth, gamma, sat,
      area;
   double leakance, deltaZ;
   int phase = 0;                                 // RESTRICTION for mfp !!!!!!!

   //WW CElem *m_ele = NULL;
   long msh_ele;
   int group, nidx;
   CRFProcess* m_pcs_cond = NULL;
   CRFProcess* m_pcs_this = NULL;
   m_pcs_this = PCSGet(convertProcessTypeToString (m_st->getProcessType()));
   m_pcs_cond = PCSGet(m_st->pcs_type_name_cond);

   area = value;
   leakance = m_st->getCoupLeakance();
   deltaZ = m_st->st_rill_height;
                                                  // phase  = 0 !!!!
   gamma = mfp_vector[0]->Density() * GRAVITY_CONSTANT;
   long msh_node_2nd;
   double const* const xyz_this (m_pcs_this->m_msh->nod_vector[cnodev->msh_node_number]->getData());
//   double y_this = m_pcs_this->m_msh->nod_vector[cnodev->msh_node_number]->Y();
//   double z_this = m_pcs_this->m_msh->nod_vector[cnodev->msh_node_number]->Z();

   msh_node_2nd = -1;                             //WW

   cond0 = leakance * deltaZ;
   cond1 = cond0;

   if (m_st->getProcessType () == FiniteElement::OVERLAND_FLOW)
   {

      ///// get number of second mesh node, provisional implementation
      double epsilon = 1.e-5;
//      double x_cond = m_pcs_cond->m_msh->nod_vector[cnodev->msh_node_number_conditional]->X();
//      double y_cond = m_pcs_cond->m_msh->nod_vector[cnodev->msh_node_number_conditional]->Y();
//      double z_cond = m_pcs_cond->m_msh->nod_vector[cnodev->msh_node_number_conditional]->Z();
      double const* const xyz_cond (m_pcs_cond->m_msh->nod_vector[cnodev->msh_node_number_conditional]->getData());

      for (size_t i = 0; i < m_pcs_cond->m_msh->nod_vector.size(); i++) {
			double const* const pnt_i(
					m_pcs_cond->m_msh->nod_vector[i]->getData());
			if (pnt_i[0] - xyz_cond[0] < epsilon) {
				if (pnt_i[1] - xyz_cond[1] < epsilon) {
					if (pnt_i[2] - (xyz_cond[2] - deltaZ) < epsilon) {
						msh_node_2nd = i;
					}
				}
			}
		}
      //////////////////////////

      nidx = m_pcs_cond->GetNodeValueIndex("PRESSURE1") + 1;

      pressure0 = m_pcs_cond->GetNodeValue(
         cnodev->msh_node_number_conditional, nidx);
      pressure1 = m_pcs_cond->GetNodeValue(msh_node_2nd, nidx);

                                                  // only one phase
      double gamma = mfp_vector[phase]->Density() * GRAVITY_CONSTANT;

      msh_ele
         = m_pcs_cond->m_msh->nod_vector[cnodev->msh_node_number_conditional]->getConnectedElementIDs()[0];
      //WW m_ele = m_pcs_cond->m_msh->ele_vector[msh_ele];
      group = m_pcs_cond->m_msh->ele_vector[msh_ele]->GetPatchIndex();

      //sat = mmp_vector[group]->SaturationCapillaryPressureFunction( -pressure0, 0);
      //cond0 *=  mmp_vector[group]->PermeabilitySaturationFunction(sat,0);

      sat = mmp_vector[group]->SaturationCapillaryPressureFunction(-pressure1);
      cond1 *= mmp_vector[group]->PermeabilitySaturationFunction(sat, phase);
      // use of relative permeability for second node (absolute perm. for top node !!!!)

      value = (pressure1 - pressure0 - deltaZ * gamma) * (cond0 + cond1)
         / (2* deltaZ * gamma);

      m_pcs_this->SetNodeValue(cnodev->msh_node_number,
         m_pcs_this->GetNodeValueIndex("COUPLING") + 1, -value);

      value *= area;
   }                                              // end overland
   else { // Richards
		///// get number of second mesh node, provisional implementation
		double epsilon = 1.e-5;
		for (size_t i = 0; i < m_pcs_this->m_msh->nod_vector.size(); i++) {
			double const* const pnt_i(
					m_pcs_this->m_msh->nod_vector[i]->getData());
			if (pnt_i[0] - xyz_this[0] < epsilon) {
				if (pnt_i[1] - xyz_this[1] < epsilon) {
					if (pnt_i[2] - (xyz_this[2] - deltaZ) < epsilon) {
						msh_node_2nd = i;
					}
				}
			}
		}
      //////////////////////////

      double inf_cap, supplyRate; //WW, rainfall;
      long
         bc_eqs_index =
         m_pcs_this->m_msh->nod_vector[cnodev->msh_node_number]->GetEquationIndex();
      double z_cond = m_pcs_cond->m_msh->nod_vector[cnodev->msh_node_number_conditional]->getData()[2];
      depth = std::max(0., m_pcs_cond->GetNodeValue(
         cnodev->msh_node_number_conditional,
         m_pcs_cond->GetNodeValueIndex("HEAD") + 1) - z_cond);

      nidx = m_pcs_this->GetNodeValueIndex("PRESSURE1") + 1;
      pressure0 = m_pcs_this->GetNodeValue(cnodev->msh_node_number, nidx);
      pressure1 = m_pcs_this->GetNodeValue(msh_node_2nd, nidx);

      msh_ele
         = m_pcs_this->m_msh->nod_vector[cnodev->msh_node_number]->getConnectedElementIDs()[0];
      //WW m_ele = m_pcs_this->m_msh->ele_vector[msh_ele];
      group = m_pcs_this->m_msh->ele_vector[msh_ele]->GetPatchIndex();

      //sat = mmp_vector[group]->SaturationCapillaryPressureFunction( -pressure0);
      //cond0 *=  mmp_vector[group]->PermeabilitySaturationFunction(sat,phase);

      sat = mmp_vector[group]->SaturationCapillaryPressureFunction(-pressure1);
      cond1 *= mmp_vector[group]->PermeabilitySaturationFunction(sat, phase);
      // use of relative permeability for second node (absolute perm. for top node !!!!)

      // calculate infiltration capacity
      /* //WW
      if (aktuelle_zeit < m_st->rainfall_duration)
         rainfall = m_st->rainfall;
      else
         rainfall = 0;
      */
      inf_cap = (depth + deltaZ - pressure1 / gamma) * (cond0 + cond1) / (2
         * deltaZ);
      supplyRate = m_st->rainfall;                //+ (depth ) / dt; // dt = timeStep

      m_pcs_this->SetNodeValue(cnodev->msh_node_number,
                                                  // update coupling variable for error estimation
         m_pcs_this->GetNodeValueIndex("COUPLING") + 1, inf_cap);

      if (inf_cap > supplyRate)
         bc_value = pressure1 - deltaZ * gamma + gamma * supplyRate * deltaZ
            * 2 / (cond0 + cond1);
      else
         bc_value = pressure1 - deltaZ * gamma + gamma * inf_cap * deltaZ
            * 2 / (cond0 + cond1);
      // bc_value = supplyRate * gamma * dt;
      /*
#ifndef NEW_EQS && !defined(USE_PETSC)
      MXRandbed(bc_eqs_index, bc_value, m_pcs_this->getEQSPointer()->b); //getEQSPointer. WW 
//else ...
#endif
*/
      value = 0;

   }                                              // end Richards

}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 02/2006 MB Implementation
 02/2006 WW Change argument
 09/2010 KR cleaned up code
 09/2010 TF commented out method
 **************************************************************************/
//double CSourceTermGroup::GetRiverNODValue(int i,CSourceTerm* m_st, long msh_node) //WW
//void GetRiverNODValue(double &value, CNodeValue* cnodev, const CSourceTerm* m_st) //WW
//{
//	double h;
//	double paraA(cnodev->node_parameterA); //HRiver
//	double paraB(cnodev->node_parameterB); //KRiverBed
//	double paraC(cnodev->node_parameterC); //WRiverBed
//	double paraD(cnodev->node_parameterD); //TRiverBed
//	double paraE(cnodev->node_parameterE); //BRiverBed
//	double NodeReachLength(value);
//	CRFProcess* m_pcs_this(NULL);
//
//	m_pcs_this = PCSGet(convertProcessTypeToString (m_st->getProcessType()));
//	//_________________________________________________________________
//	//paraA jetzt aus dem Prozess Overland Flow
//	if (m_st->isCoupled()) {
//
//		CRFProcess* m_pcs_cond = NULL;
//		m_pcs_cond = PCSGet(m_st->pcs_type_name_cond);
//
//		int nidx = m_pcs_cond->GetNodeValueIndex(m_st->pcs_pv_name_cond) + 1;
//		//WW    long node_cond = group_vector[i]->msh_node_number_conditional;
//		long node_cond = cnodev->msh_node_number_conditional; //WW
//		paraA = m_pcs_cond->GetNodeValue(node_cond, nidx);
//	}
//
//	double RiverConductance = paraB * paraC * NodeReachLength / (paraD - paraE);
//	int nidx1 = m_pcs_this->GetNodeValueIndex("HEAD") + 1;
//	h = m_pcs_this->GetNodeValue(cnodev->msh_node_number, nidx1);
//
//	if (h > paraD) { //HAquiver > BRiverBed
//		//q = (RiverConductance * HRiver)   -  (RiverConductance * HAquifer)
//		value = RiverConductance * paraA;
//		MXInc(cnodev->msh_node_number, cnodev->msh_node_number,	RiverConductance);
//	}
//	if (h < paraD) { //HAquiver < BRiverBed
//		//q = (RiverConductance * HRiver)   - (RiverConductance * BRiverBed)
//		value = RiverConductance * (paraA - paraD);
//	}
//	if (h == paraE)
//		value = 0.;
//	//_________________________________________________________________
//	//Safe Flux values
//	int nidxFLUX = m_pcs_this->GetNodeValueIndex("FLUX") + 1;
//	double flux = value / NodeReachLength; //fluxes in m^2/s !!
//	m_pcs_this->SetNodeValue(cnodev->msh_node_number, nidxFLUX, flux);
//}

/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 02/2006 MB Implementation
 02/2006 WW Change argument
 **************************************************************************/
//double CSourceTermGroup::GetCriticalDepthNODValue(CNodeValue* cnodev,CSourceTerm* m_st, long msh_node)
void GetCriticalDepthNODValue(double &value, CSourceTerm* m_st, long msh_node)
{
   double value_jacobi;
   double width, flowdepth, flowdepth3, flowdepth3_epsilon;
   long msh_ele;
   double epsilon = 1.e-7;                        // like in pcs->assembleParabolicEquationNewton

   CRFProcess* m_pcs_this = NULL;
   m_pcs_this = PCSGet(convertProcessTypeToString(m_st->getProcessType()));
   long nidx1 = m_pcs_this->GetNodeValueIndex("HEAD") + 1;
   flowdepth = m_pcs_this->GetNodeValue(msh_node, nidx1)
      - m_pcs_this->m_msh->nod_vector[msh_node]->getData()[2] - m_st->st_rill_height;

   if (flowdepth < 0.0)
   {
      value = 0;
      m_pcs_this->SetNodeValue(msh_node,
         m_pcs_this->GetNodeValueIndex("FLUX") + 0, -value);
   }
   else
   {
      flowdepth3 = MathLib::fastpow(flowdepth, 3);
      flowdepth3_epsilon = MathLib::fastpow(flowdepth + epsilon, 3);

      width = value;
      if (m_pcs_this->m_msh->GetMaxElementDim() == 1)
      {
         msh_ele
            = m_pcs_this->m_msh->nod_vector[msh_node]->getConnectedElementIDs()[0];
         int group = m_pcs_this->m_msh->ele_vector[msh_ele]->GetPatchIndex();
         width = mmp_vector[group]->overland_width;
      }

      value = -sqrt(GRAVITY_CONSTANT * flowdepth3) * width;

      value_jacobi = sqrt(GRAVITY_CONSTANT * flowdepth3_epsilon) * width
         + value;
      /*
#ifndef NEW_EQS && !defined(USE_PETSC)
                                                  // write source term into jacobi
      MXInc(msh_node, msh_node, value_jacobi / epsilon);
// else ...
#endif
*/
      m_pcs_this->SetNodeValue(msh_node,
         m_pcs_this->GetNodeValueIndex("FLUX") + 0, -value);

   }
}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 02/2006 MB JOD Implementation
 06/2007 JOD 2D case with slope in st-file
 **************************************************************************/
void GetNormalDepthNODValue(double &value, CSourceTerm* st, long msh_node)
{
   //WW  int AnzNodes = 0;
   //WW  double Haverage = 0;
   CRFProcess* pcs_this (PCSGet(st->getProcessType()));
   CFEMesh* mesh (pcs_this->m_msh);

   double value_for_jacobi, S_0;
   double epsilon = 1.e-7;                        // pcs->assembleParabolicEquationNewton !!!!!!!!!

   long msh_ele = mesh->nod_vector[msh_node]->getConnectedElementIDs()[0];
   CElem *m_ele = mesh->ele_vector[msh_ele];
   int group = mesh->ele_vector[msh_ele]->GetPatchIndex();
   double width = mmp_vector[group]->overland_width;
   double fric_coef = mmp_vector[group]->friction_coefficient;
   double slope_exp = mmp_vector[group]->friction_exp_slope;
   double depth_exp = mmp_vector[group]->friction_exp_depth;
   if (st->getNormalDepthSlope() == -1)           // TF: WARNING comparison of double via ==
   {
      if (mesh->GetMaxElementDim() > 1)
         std::cout << "!!!!! give slope for NORMAL DEPTH in st-file !!!!!"
            << "\n";

      double elementlength = sqrt(MathLib::sqrDist(m_ele->GetNode(1)->getData(), m_ele->GetNode(0)->getData()));
//    		  (MathLib::fastpow(m_ele->GetNode(1)->X()- m_ele->GetNode(0)->X(), 2)
//    		  + MathLib::fastpow(m_ele->GetNode(1)->Y()-m_ele->GetNode(0)->Y(), 2)
//    		  + MathLib::fastpow(m_ele->GetNode(1)->Z() - m_ele->GetNode(0)->Z(), 2));
      S_0 = (m_ele->GetNode(1)->getData()[2] - m_ele->GetNode(0)->getData()[2]) / elementlength;
      if (S_0 < 0)
         S_0 = -S_0;
   } else
   S_0 = st->getNormalDepthSlope();

   double flowdepth = pcs_this->GetNodeValue(msh_node, 1)
      - mesh->nod_vector[msh_node]->getData()[2] - st->st_rill_height;
   double flowdepth_epsilon = flowdepth + epsilon;
   if (flowdepth < 0.0)
   {
      flowdepth = 0.0;
      flowdepth_epsilon = 0.0;
   }

   double temp = width * fric_coef * pow(S_0, slope_exp);
   if (mmp_vector[group]->channel == 1)
   {
      value = -pow(flowdepth * width / (2 * flowdepth + width), depth_exp)
         * flowdepth * temp;
      value_for_jacobi = pow(flowdepth_epsilon * width / (2
         * flowdepth_epsilon + width), depth_exp) * flowdepth_epsilon
         * temp + value;
   }
   else
   {
      value = -pow(flowdepth, depth_exp + 1) * temp;
      value_for_jacobi = pow(flowdepth_epsilon, depth_exp + 1) * temp + value;
   }
   /*
#ifndef NEW_EQS && !defined(USE_PETSC)
                                                  // write source term into jacobi
   MXInc(msh_node, msh_node, value_for_jacobi / epsilon);
// else
#endif
*/
   pcs_this->SetNodeValue(msh_node, pcs_this->GetNodeValueIndex("FLUX")
      + 0, -value);
}
//#endif

/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 11/2007 JOD Implementation
 **************************************************************************/
void GetNODValue(double& value, CNodeValue* cnodev, CSourceTerm* st)
{
#if !defined(USE_PETSC)
	//&& !defined(NEW_EQS) // && defined(other parallel libs)//03~04.3012. WW
  //#ifndef NEW_EQS                                //WW. 06.11.2008
   if (st->isCoupled())
      GetCouplingNODValue(value, st, cnodev);
   else if (st->isAnalytical())
   {
      //WW      m_st_group->m_msh = m_msh;
                                                  //WW
      value = st->GetAnalyticalSolution(cnodev->msh_node_number);
      //WW         value = m_st_group->GetAnalyticalSolution(m_st,msh_node,(string)function_name[j]);
   }

   //	if (cnodev->node_distype == 5) // River Condition
   //	if (cnodev->getProcessDistributionType() == RIVER)
   //		GetRiverNODValue(value, cnodev, st); //MB
   //	if (cnodev->node_distype == 6) // CriticalDepth Condition
                                                  // CriticalDepth Condition
   if (cnodev->getProcessDistributionType() == FiniteElement::CRITICALDEPTH)
                                                  //MB
      GetCriticalDepthNODValue(value, st, cnodev->msh_node_number);
   //	if (cnodev->node_distype == 8) // NormalDepth Condition JOD
   if (cnodev->getProcessDistributionType() == FiniteElement::NORMALDEPTH)
                                                  //M
      GetNormalDepthNODValue(value, st, cnodev->msh_node_number);
#endif
   //	if (cnodev->node_distype == 10) // Philip infiltration JOD
   //		GetPhilipNODValue(value, st);
   //	if (cnodev->node_distype == 11) // Green_Ampt infiltration JOD
   if (cnodev->getProcessDistributionType() == FiniteElement::GREEN_AMPT)
      GetGreenAmptNODValue(value, st, cnodev->msh_node_number);

   //TN - Test flux with heat transfer coefficient
   if (st->getProcessPrimaryVariable() == FiniteElement::TEMPERATURE && st->getProcessDistributionType() == FiniteElement::TRANSFER_SURROUNDING) {
	   GetNODHeatTransfer(value, st, cnodev->geo_node_number);
	   //value = st->getTransferCoefficient()*cnodev->node_area*(cnodev->node_value - st->getValueSurrounding());
   }
   else if (st->getProcessPrimaryVariable() == FiniteElement::TEMPERATURE2 && st->getProcessDistributionType() == FiniteElement::TRANSFER_SURROUNDING) {
	   GetNODHeatTransfer(value, st, cnodev->geo_node_number);
	   //value = st->getTransferCoefficient()*cnodev->node_area*(cnodev->node_value - st->getValueSurrounding());
   }

}

/**************************************************************************
 FEMLib-Method:
 Task: Compute heat flux for heat transfer boundary condition
 Programing:
 04/2013 TN Implementation
 last modified:
 **************************************************************************/
void GetNODHeatTransfer(double& value, CSourceTerm* st, long geo_node){
   CRFProcess* m_pcs_this = NULL;
   double poro;

   //Get process type
   m_pcs_this = PCSGet(convertProcessTypeToString(st->getProcessType()));
   //Get Mesh
   CFEMesh* mesh (m_pcs_this->m_msh);
   
   //Get number of conneted elements
   size_t number_of_connected_elements = mesh->nod_vector[geo_node]->getConnectedElementIDs().size();

   poro = 0.0;
   double geo_area = 0.0;

   //loop over connected elements and get average porosity
   for (size_t i=0;i<number_of_connected_elements;i++){
	   long msh_ele = mesh->nod_vector[geo_node]->getConnectedElementIDs()[i];
	   int group = mesh->ele_vector[msh_ele]->GetPatchIndex();
	   poro += mmp_vector[group]->porosity;
	   geo_area += mmp_vector[group]->geo_area;
   }
   poro /= number_of_connected_elements;
   geo_area /= number_of_connected_elements;

   //if (mesh->isAxisymmetry() && mesh->GetMaxElementDim()!=1) //For axisymmetric 2D meshes geometry area is irrelevant
	  // geo_area = 1.0;
   
   

   //Get index of primary variable
   long nidx1 = m_pcs_this->GetNodeValueIndex(convertPrimaryVariableToString(st->getProcessPrimaryVariable())) + 1;

   //Get current primary variable value at that node
   double temp = m_pcs_this->GetNodeValue(geo_node, nidx1);

   //Find position of current node in st vectors
   size_t i;
   for (i=0; i<st->get_node_value_vectorArea().size(); i++){
	   if (geo_node == st->st_node_ids[i])
		   break;
   }

   value = st->getTransferCoefficient()*(st->getValueSurrounding() - temp);
   value *= st->get_node_value_vectorArea()[i]*geo_area;

   if (st->getProcessPrimaryVariable() == FiniteElement::TEMPERATURE2)
	   value *= (1.0 - poro);
   else if (st->getProcessPrimaryVariable() == FiniteElement::TEMPERATURE)
	   value *= poro;

}

/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 01/2005 OK Implementation
 last modified:
 **************************************************************************/
void STGroupDelete(std::string pcs_type_name, std::string pcs_pv_name)
{
   CSourceTermGroup* m_st_group = NULL;
   std::list<CSourceTermGroup*>::const_iterator p = st_group_list.begin();
   while (p != st_group_list.end())
   {
      m_st_group = *p;
      if ((m_st_group->pcs_type_name.compare(pcs_type_name) == 0)
         && (m_st_group->pcs_pv_name.compare(pcs_pv_name) == 0))
      {
         delete m_st_group;
         st_group_list.remove(m_st_group);
         return;
      }
      ++p;
   }
}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 11/2007 JOD Implementation
 modification:
 05/2010
 **************************************************************************/
void CSourceTermGroup::SetPLY(CSourceTerm* st, int ShiftInNodeVector)
{
	CGLPolyline* old_ply (GEOGetPLYByName(st->geo_name));
	if (old_ply) {
		std::vector<long> ply_nod_vector;
		std::vector<long> ply_nod_vector_cond;
		std::vector<double> ply_nod_val_vector;
		std::vector<double> ply_nod_vector_cond_length;

		double min_edge_length (m_msh->getMinEdgeLength());
		m_msh->setMinEdgeLength (old_ply->epsilon);
		m_msh->GetNODOnPLY(static_cast<const GEOLIB::Polyline*>(st->getGeoObj()), ply_nod_vector, true);
		m_msh->setMinEdgeLength (min_edge_length);  // reset value

		if (st->isCoupled() && (st->pcs_type_name_cond == "OVERLAND_FLOW" || st->pcs_type_name_cond2 == "OVERLAND_FLOW")  )
		// removed !!!!!		&& !st->isConnectedGeometry())  // if both combined connected nodes are set below
			SetPolylineNodeVectorConditional(st, ply_nod_vector, ply_nod_vector_cond);

	    if(st->hasThreshold()) // JOD 2018-02-20
	    {  // only point supported
		  st->msh_node_number_threshold = m_msh->GetNODOnPNT(
						static_cast<const GEOLIB::Point*>(st->geoInfo_threshold->getGeoObj()));
	    }

	    if(st->calculatedFromStorageRate()) // JOD 2018-02-22
	    {  // only point supported
		  st->storageRate.inlet_msh_node_numbers.push_back(m_msh->GetNODOnPNT(
						static_cast<const GEOLIB::Point*>(st->geoInfo_storageRateInlet->getGeoObj())));
		  st->storageRate.outlet_msh_node_numbers.push_back(m_msh->GetNODOnPNT(
						static_cast<const GEOLIB::Point*>(st->geoInfo_storageRateOutlet->getGeoObj())));
	    }

		SetPolylineNodeValueVector(st, ply_nod_vector, ply_nod_vector_cond, ply_nod_val_vector);

		if(st->verbosity)
		{
			std::cout << "\t" << ply_nod_vector.size() << " nodes with total value of " << std::accumulate(ply_nod_val_vector.begin(), ply_nod_val_vector.end(), 0.) << '\n';
			for(int i=0; i< ply_nod_vector.size(); i++)
			{
				std::cout << "\t\t" << i << ": " << ply_nod_vector[i] << " with value " << ply_nod_val_vector[i];
			  	if (st->isConnectedGeometry() &&
					  st->getConnectedGeometryCouplingType() != 2) // not borehole with given primary variable
					std::cout << " - connected to " << ply_nod_vector_cond[i];
				std::cout << '\n';
			}
		}

		if(st->scaling_verbosity && st->scaling_mode == 1)
		{
			std::cout << "Scaling mode: " << st->scaling_mode << '\n';
			for(long i=0; i < ply_nod_vector.size(); ++i)
			{
				std::cout << '\t' << ply_nod_vector[i] << ":\t" << ply_nod_val_vector[i] << '\n'; 
			}
		}

		if (st->distribute_volume_flux)   // 5.3.07 JOD
			DistributeVolumeFlux(st, ply_nod_vector, ply_nod_val_vector);

		st->st_node_ids.clear();
		st->st_node_ids.resize(ply_nod_vector.size());
		st->st_node_ids = ply_nod_vector;

		if (st->isConnectedGeometry() &&  // JOD 10/2018
				st->getConnectedGeometryCouplingType() != 2 // not borehole with given primary variable
				&& (st->pcs_type_name_cond != "OVERLAND_FLOW" &&  st->pcs_type_name_cond2 != "OVERLAND_FLOW"))
		{

			  st->SetPolylineNodeVectorConnected(ply_nod_vector, ply_nod_vector_cond);

			  if(st->average_mode == 1)
			  {
				  CRFProcess* m_pcs(PCSGet(pcs_type_name));

				  ply_nod_vector_cond_length.resize(ply_nod_vector_cond.size());
				  std::fill(ply_nod_vector_cond_length.begin(), ply_nod_vector_cond_length.end(), 1.);

				  if (m_msh->GetMaxElementDim() == 1)
					  FiniteElement::DomainIntegration(m_pcs, ply_nod_vector_cond, ply_nod_vector_cond_length);
				  else FiniteElement::EdgeIntegration(m_pcs->m_msh, ply_nod_vector_cond, ply_nod_vector_cond_length,
						  	  	  st->getProcessDistributionType(), st->getProcessPrimaryVariable(),
								  true, false, 0);//
								//bc->ignore_axisymmetry, st->isPressureBoundaryCondition(), st->scaling_mode);
			  }
			  else if(st->average_mode == -1)  // take nearest point
			  {
				CRFProcess* m_pcs(PCSGet(pcs_type_name));
				std::vector<long> ply_nod_vector_cond_new;
				for(std::size_t i = 0; i < ply_nod_vector.size(); ++i)
				{
			  		double distance_min = 1e10;
					int j_nearest = -1;

					for(std::size_t j = 0; j < ply_nod_vector_cond.size(); ++j)
					{
						const double distance_current = 
							(m_pcs->m_msh->nod_vector[ply_nod_vector[i]]->X() - m_pcs->m_msh->nod_vector[ply_nod_vector_cond[j]]->X()) *
							(m_pcs->m_msh->nod_vector[ply_nod_vector[i]]->X() - m_pcs->m_msh->nod_vector[ply_nod_vector_cond[j]]->X()) +
							(m_pcs->m_msh->nod_vector[ply_nod_vector[i]]->Y() - m_pcs->m_msh->nod_vector[ply_nod_vector_cond[j]]->Y()) *
							(m_pcs->m_msh->nod_vector[ply_nod_vector[i]]->Y() - m_pcs->m_msh->nod_vector[ply_nod_vector_cond[j]]->Y()) +
							(m_pcs->m_msh->nod_vector[ply_nod_vector[i]]->Z() - m_pcs->m_msh->nod_vector[ply_nod_vector_cond[j]]->Z()) *
							(m_pcs->m_msh->nod_vector[ply_nod_vector[i]]->Z() - m_pcs->m_msh->nod_vector[ply_nod_vector_cond[j]]->Z());
						if(distance_min > distance_current)
						{
							j_nearest = j;
						       	distance_min = distance_current;	
						}
					}
					ply_nod_vector_cond_new.push_back(ply_nod_vector_cond[j_nearest]);
				}				

				ply_nod_vector_cond = ply_nod_vector_cond_new;

			  }
		}
		

		if (st->isConstrainedST())
		{
			for (std::size_t i(0); i < st->st_node_ids.size(); i++)
			{
				st->_constrainedSTNodesIndices.push_back(-1);
				for (std::size_t i(0); i < st->getNumberOfConstrainedSTs(); i++)
					st->pushBackConstrainedSTNode(i, false);
			}
		}


		if (st->everyoneWithEveryone)  // JOD 8/2015   quick'n'dirty to test approach
			  {
		  	 	  std::cout << "\tConnect every node with every node of other surface\n";
				  std::vector<long>::iterator pos;
				  std::vector<double> ply_nod_val_vector_cond_original;


				  SetPolylineNodeValueVector(st, ply_nod_vector_cond, ply_nod_vector_cond /* nod used */, ply_nod_val_vector_cond_original);

		  		  const double total_val_cond_original = std::accumulate(ply_nod_val_vector_cond_original.begin(), ply_nod_val_vector_cond_original.end(), 0.);

	  	  		  if (st->verbosity && st->isConnectedGeometry() &&
			  		st->getConnectedGeometryCouplingType() != 2) // not borehole with given primary variable
					std::cout << "\t\t\t" << ply_nod_vector_cond.size() << " connected nodes with total value of " << total_val_cond_original << '\n';


				  int nod_vector_size = (int)ply_nod_vector.size();
				  int nod_vector_cond_size = (int)ply_nod_vector_cond.size();

				  for (int i = 0; i < nod_vector_size; i++)  // extend nod_vector
				  {
					  for (int j = 1; j < nod_vector_cond_size; j++)
					  {
						  pos = ply_nod_vector.begin() + i * nod_vector_cond_size + j;
						  ply_nod_vector.insert(pos, ply_nod_vector[i * nod_vector_cond_size + j - 1]);
					  }
				  }

				  for (int i = 1; i < nod_vector_size; i++)  // extend nod_vector_cond
				  {
					  for (int j = 0; j < nod_vector_cond_size; j++)
					  {
						  ply_nod_vector_cond.push_back(ply_nod_vector_cond[j]);
					  }
				  }

				  // extend nod_val_vector
				  std::vector<double> ply_nod_val_vector_original(ply_nod_val_vector);
				  ply_nod_val_vector.resize(nod_vector_size * nod_vector_cond_size);

				  for (int i = 0; i < nod_vector_size; i++)
				  {
					  for (int j = 0; j < nod_vector_cond_size; j++)
					  {
						  ply_nod_val_vector[i*nod_vector_cond_size + j] = ply_nod_val_vector_original[i] * ply_nod_val_vector_cond_original[j] / (total_val_cond_original);
					  }
				  }

		  		  if(st->verbosity)
		  		  {
		  			std::cout << "\tNow " << ply_nod_vector.size() << " nodes with total value of " << 
						std::accumulate(ply_nod_val_vector.begin(), ply_nod_val_vector.end(), 0.) << '\n';
		  		  	if(st->verbosity > 1)
		  				for (int i = 0; i < ply_nod_vector.size(); i++)
		  				{
				  			std::cout << "\t\t" << i << ": " << ply_nod_vector[i] << " connected to " << ply_nod_vector_cond[i] << " with value "<< ply_nod_val_vector[i] << '\n';
		  				}
		  		  }
			  }
		/////
		if(st->ogs_WDC != nullptr)  // for (3D) ATES with polyline wells
		{  	  // point for measurements
			std::vector<long> ply_nod_vector_well2;

			m_msh->setMinEdgeLength (old_ply->epsilon);
			m_msh->GetNODOnPLY(static_cast<const GEOLIB::Polyline*>(
					st->geoInfo_wellDoublet_well2_aquifer->getGeoObj()),
					ply_nod_vector_well2, true);
			m_msh->setMinEdgeLength (min_edge_length);  // reset value

			std::vector<double> well2_mesh_node_values;

			SetPolylineNodeValueVector(st, ply_nod_vector_well2,
					ply_nod_vector_well2, // cond (not used)
					well2_mesh_node_values);

			double total_value = std::accumulate(ply_nod_val_vector.begin(),
								ply_nod_val_vector.end(), 0.);
			for(auto& value: ply_nod_val_vector) value /= total_value;

			total_value = std::accumulate(well2_mesh_node_values.begin(),
					well2_mesh_node_values.end(), 0.);
			for(auto& value: well2_mesh_node_values) value /= total_value;



			st->ogs_WDC->set_doublet_mesh_nodes({
				std::vector<size_t>(ply_nod_vector.begin(), ply_nod_vector.end()),  // well1
				std::vector<size_t>(ply_nod_vector_well2.begin(), ply_nod_vector_well2.end()), // well2
				std::vector<size_t>(ply_nod_vector.begin(), ply_nod_vector.end()),
												// nod_val->msh_node_number  // heatExchanger
				std::vector<double>(ply_nod_val_vector.begin(), ply_nod_val_vector.end()),
												// well 1 - node area
				std::vector<double>(well2_mesh_node_values.begin(), well2_mesh_node_values.end()),
												// well 2 - node area
				std::vector<double>(ply_nod_val_vector.begin(), ply_nod_val_vector.end()),
																// heatExchanger - node area
			});

			CRFProcess* pcs_liquid = PCSGet(FiniteElement::LIQUID_FLOW);
			// warm well 1 - LIQUID_FLOW
			for(long i=0; i < ply_nod_vector.size(); ++i)
			 {
				 CNodeValue *nod_val_liquid_well (new CNodeValue());
				 nod_val_liquid_well->msh_node_number = ply_nod_vector[i];
				 nod_val_liquid_well->CurveIndex = st->CurveIndex;
				 nod_val_liquid_well->geo_node_number =
						 nod_val_liquid_well->msh_node_number - ShiftInNodeVector;
				 nod_val_liquid_well->node_value = st->geo_node_value * ply_nod_val_vector[i];
				 nod_val_liquid_well->tim_type_name = st->tim_type_name;

				 pcs_liquid->st_node_value.push_back(nod_val_liquid_well);
				 pcs_liquid->st_node.push_back(st);
			 }

			// cold well 2 - LIQUID_FLOW
			for(long i=0; i < ply_nod_vector_well2.size(); ++i)
			 {
				 CNodeValue *nod_val_liquid_well (new CNodeValue());
				 nod_val_liquid_well->msh_node_number = ply_nod_vector_well2[i];
				 nod_val_liquid_well->CurveIndex = st->CurveIndex;
				 nod_val_liquid_well->geo_node_number =
						 nod_val_liquid_well->msh_node_number - ShiftInNodeVector;
				 nod_val_liquid_well->node_value = -st->geo_node_value * well2_mesh_node_values[i];  // ! negative
				 nod_val_liquid_well->tim_type_name = st->tim_type_name;

				 pcs_liquid->st_node_value.push_back(nod_val_liquid_well);
				 pcs_liquid->st_node.push_back(st);
			 }


			std::copy (ply_nod_vector_well2.begin(), ply_nod_vector_well2.end(), std::back_inserter(ply_nod_vector));
			std::transform(well2_mesh_node_values.begin(), well2_mesh_node_values.end(), well2_mesh_node_values.begin(),
					[&](double& c){return -c;});
			std::copy (well2_mesh_node_values.begin(), well2_mesh_node_values.end(), std::back_inserter(ply_nod_val_vector));
		
	   }

		if(st->borehole_mode == 0)  // conductive with peaceman
		{
			double factor, radius_e;
			CRFProcess* m_pcs(PCSGet(pcs_type_name));

			for(size_t i=0; i < ply_nod_vector.size(); ++i)
			{
				CalculatePeaceman(st, m_pcs, ply_nod_vector[i],
						m_pcs->m_msh->nod_vector[ply_nod_vector[i]]->getConnectedElementOnPolyineIDs(
								ply_nod_vector, m_pcs->m_msh->ele_vector)
						, factor, radius_e);

				factor /= std::log(radius_e / st->borehole_data.radius);
				ply_nod_val_vector[i] *= factor;
				
         			// to use advective HEAT_TRANSPORT with conductive (peaceman) LIQUID_FLOW 
				if(m_pcs->ST_factor_kept.find(ply_nod_vector[i]) != m_pcs->ST_factor_kept.end())
                                       	m_pcs->ST_factor_kept[ply_nod_vector[i]] = m_pcs->ST_factor_kept[ply_nod_vector[i]] + factor;
                        	else
                               		m_pcs->ST_factor_kept[ply_nod_vector[i]] = factor;
			}
		}

		st->SetNodeValues(ply_nod_vector, ply_nod_vector_cond, ply_nod_val_vector, ply_nod_vector_cond_length, ShiftInNodeVector);
	} // end polyline
}


void CSourceTermGroup::SetPNT(CRFProcess* pcs, CSourceTerm* st,
const int ShiftInNodeVector)
{
 	 //05.2012. WW
 	 long node_id = m_msh->GetNODOnPNT(static_cast<const GEOLIB::Point*>(st->getGeoObj()));
 	 if(node_id < 0)
 	   return;
 	  CNodeValue *nod_val (new CNodeValue());

 	  // TF removed some checks - check validity of data while reading data

 	  nod_val->msh_node_number = m_msh->GetNODOnPNT(static_cast<const GEOLIB::Point*>(st->getGeoObj())) + ShiftInNodeVector;
 	  nod_val->CurveIndex = st->CurveIndex;
 	                                                 //WW
 	  nod_val->geo_node_number = nod_val->msh_node_number - ShiftInNodeVector;
 	  nod_val->node_value = st->geo_node_value;
 	  nod_val->tim_type_name = st->tim_type_name;

 	  if(st->isConnectedGeometry() && // JOD 2018-02-20
 	    st->getConnectedGeometryCouplingType() != 2)
 	  {
 	          nod_val->msh_node_number_conditional = m_msh->GetNODOnPNT(
 	       		   static_cast<const GEOLIB::Point*>(st->geoInfo_connected->getGeoObj()));
 	  }

 	  if(st->hasThreshold()) // JOD 2018-02-20
 	  {  // only point supported
 	         st->msh_node_number_threshold = m_msh->GetNODOnPNT(
 	       				static_cast<const GEOLIB::Point*>(st->geoInfo_threshold->getGeoObj()));
 	  }


 	  if(st->calculatedFromStorageRate()) // JOD 2018-02-22
 	  {  // only point supported
 	         st->storageRate.inlet_msh_node_numbers.push_back(m_msh->GetNODOnPNT(
 	       				static_cast<const GEOLIB::Point*>(st->geoInfo_storageRateInlet->getGeoObj())));
 	         st->storageRate.outlet_msh_node_numbers.push_back(m_msh->GetNODOnPNT(
 	       				static_cast<const GEOLIB::Point*>(st->geoInfo_storageRateOutlet->getGeoObj())));
 	          st->storageRate.inlet_msh_node_areas.push_back(1.);
 	          st->storageRate.outlet_msh_node_areas.push_back(1.);
 	          st->storageRate.inlet_totalArea = 1.;
 	          st->storageRate.outlet_totalArea = 1.;
 	  }

 	  if (st->getProcessDistributionType() == FiniteElement::CRITICALDEPTH)
 	  {
 	     //	if (st->dis_type_name.compare("CRITICALDEPTH") == 0) {
 	     nod_val->setProcessDistributionType (st->getProcessDistributionType());
 	     nod_val->node_area = 1.0;
 	     std::cout << "      - Critical depth" << "\n";
 	  }

 	  if (st->getProcessDistributionType() == FiniteElement::NORMALDEPTH)
 	  {
 	     nod_val->setProcessDistributionType (st->getProcessDistributionType());
 	     nod_val->node_area = 1.0;
 	         std::cout << "      - Normal depth" << "\n";
 	  }

 	  //	if (st->dis_type_name.compare("PHILIP") == 0) { // JOD
 	  //		nod_val->node_distype = 10;
 	  //		nod_val->node_area = 1.0;
 	  //	}
 	  // Added by CB
 	  if (st->getProcessDistributionType() == FiniteElement::CONSTANT_NEUMANN)
 	  {
 	    long msh_ele = m_msh->nod_vector[nod_val->msh_node_number]->getConnectedElementIDs()[0];
 	    if (mmp_vector[m_msh->ele_vector[msh_ele]->GetPatchIndex()]->GetGeoDimension() < 3)
 	      nod_val->node_value *= mmp_vector[m_msh->ele_vector[msh_ele]->GetPatchIndex()]->geo_area;
 	  }

 	  if (st->getProcessDistributionType() == FiniteElement::GREEN_AMPT)
 	  {
 	     nod_val->setProcessDistributionType (st->getProcessDistributionType());
 	     nod_val->node_area = 1.0;
 	         std::cout << "      - Green-Ampt" << "\n";
 	  }

 	  if (st->getProcessDistributionType() == FiniteElement::SYSTEM_DEPENDENT)
 	  {
 	     nod_val->setProcessDistributionType (st->getProcessDistributionType());
 	     pcs->compute_domain_face_normal = true;     //WW
 	     CElem* elem = NULL;
 	     CNode* cnode = NULL;                        //WW
 	     for (size_t i = 0; i < m_msh->ele_vector.size(); i++)
 	     {
 	        elem = m_msh->ele_vector[i];
 	        if (!elem->GetMark())
 	           continue;
 	        int nn = elem->GetNodesNumber(m_msh->getOrder());
 	        for (long j = 0; j < nn; j++)
 	        {
 	           cnode = elem->GetNode(j);             //WW
 	           if (cnode->GetIndex() == (size_t)st->geo_node_number)
 	              st->element_st_vector.push_back(i);
 	        }
 	     }
 	  }

 	  if (st->getProcessDistributionType() == FiniteElement::TRANSFER_SURROUNDING) { //TN - Belegung mit Fl�chenelementen
 	       	st->node_value_vectorArea.resize(1);
 	           st->node_value_vectorArea[0] = 1.0;
 	       	//nod_val->node_value = 0.0;
 	       	////Get process type
 	       	//CRFProcess* m_pcs_this = PCSGet(convertProcessTypeToString(st->getProcessType()));
 	       	////Get Mesh
 	       	//CFEMesh* mesh (m_pcs_this->m_msh);
 	       	//long msh_ele = mesh->nod_vector[nod_val->geo_node_number]->getConnectedElementIDs()[0];
 	        //   int group = mesh->ele_vector[msh_ele]->GetPatchIndex();

 	       	//st->node_value_vectorArea[0] = mmp_vector[group]->geo_area;
 	       	if (m_msh->isAxisymmetry() && m_msh->GetMaxElementDim() == 1)
 	       		st->node_value_vectorArea[0] *= m_msh->nod_vector[nod_val->geo_node_number]->X();
 	       			//2pi is mulitplicated during the integration process
 	  }

 	  if (st->getProcessDistributionType()==FiniteElement::RECHARGE)	//MW
 	  {
 	          nod_val->setProcessDistributionType (st->getProcessDistributionType());
 	  }

 	  st->st_node_ids.push_back(nod_val->geo_node_number);

 	  if (st->isConstrainedST())
 	  {
 	          st->_constrainedSTNodesIndices.push_back(-1);
 	          for (std::size_t i(0); i < st->getNumberOfConstrainedSTs(); i++)
 	       	   st->pushBackConstrainedSTNode(i,false);
 	  }

 	  if(st->ogs_WDC != nullptr)
 	  {  	  // point for measurements
 	        st->ogs_WDC->set_doublet_mesh_nodes({
 	       	 std::vector<size_t>(1, m_msh->GetNODOnPNT(
 	       			 static_cast<const GEOLIB::Point*>(
 	       					 st->geoInfo_wellDoublet_well1_aquifer->getGeoObj()))),  // well1
 	       			 std::vector<size_t>(1, m_msh->GetNODOnPNT(
 	       					 static_cast<const GEOLIB::Point*>(
 	       							 st->geoInfo_wellDoublet_well2_aquifer->getGeoObj()))),  // well2
 	       			  std::vector<size_t>{static_cast<size_t>(nod_val->msh_node_number)},  // heatExchanger
 	       				std::vector<double>(1, 1.),  // well 1 - node area
 	       				std::vector<double>(1, 1.), // well 2 - node area
 	       				std::vector<double>(1, 1.) // heatexchanger - node area
 	       			  });
 	         // liquid flow BCs
 	         CRFProcess* pcs_liquid = PCSGet(FiniteElement::LIQUID_FLOW);

 	         std::vector<long> liquidBC_mesh_nodes;
 	         std::vector<double> liquidBC_mesh_node_values;
 	         double total_value;

 	       	 // warm well 1
 	       	 if(st->geoInfo_wellDoublet_well1_liquidBC->getGeoType() == GEOLIB::POINT)
 	       	 {
 	       		 liquidBC_mesh_nodes.push_back(m_msh->GetNODOnPNT(static_cast<const GEOLIB::Point*>(
 	       							  st->geoInfo_wellDoublet_well1_liquidBC->getGeoObj())) + ShiftInNodeVector);
 	       		 liquidBC_mesh_node_values.resize(liquidBC_mesh_nodes.size(), 1.);
 	       	 }
 	       	 else if(st->geoInfo_wellDoublet_well1_liquidBC->getGeoType() == GEOLIB::POLYLINE)
 	       	 {
 	       		 m_msh->GetNODOnPLY(static_cast<const GEOLIB::Polyline*>(
 	       				 st->geoInfo_wellDoublet_well1_liquidBC->getGeoObj()), liquidBC_mesh_nodes, true);

 	       		 liquidBC_mesh_node_values.resize(liquidBC_mesh_nodes.size(), 1.);
 	       		 FiniteElement::EdgeIntegration(m_msh, liquidBC_mesh_nodes, liquidBC_mesh_node_values,
 	       				 st->getProcessDistributionType(), st->getProcessPrimaryVariable(), 
 	       				 st->ignore_axisymmetry, st->isPressureBoundaryCondition(), st->scaling_mode);
 	       		 total_value = std::accumulate(liquidBC_mesh_node_values.begin(), liquidBC_mesh_node_values.end(), 0.);
 	       		 for(auto& value: liquidBC_mesh_node_values) value /= total_value;
 	       	 }

 	       	 for(long i=0; i < liquidBC_mesh_nodes.size(); ++i)
 	       	 {
 	       		 CNodeValue *nod_val_liquid_well (new CNodeValue());
 	       		 nod_val_liquid_well->msh_node_number = liquidBC_mesh_nodes[i];
 	       		 nod_val_liquid_well->CurveIndex = st->CurveIndex;
 	       		 nod_val_liquid_well->geo_node_number = nod_val_liquid_well->msh_node_number - ShiftInNodeVector;
 	       		 nod_val_liquid_well->node_value = st->geo_node_value * liquidBC_mesh_node_values[i];
 	       		 nod_val_liquid_well->tim_type_name = st->tim_type_name;

 	       		 pcs_liquid->st_node_value.push_back(nod_val_liquid_well);
 	       		 pcs_liquid->st_node.push_back(st);
 	       	 }

 	       	 liquidBC_mesh_nodes.clear();
 	       	 //////// cold well 2 - copy from lines above except that node_value is negative
 	       	 if(st->geoInfo_wellDoublet_well2_liquidBC->getGeoType() == GEOLIB::POINT)
 	       	 {
 	       		 liquidBC_mesh_nodes.push_back(m_msh->GetNODOnPNT(static_cast<const GEOLIB::Point*>(
 	       							  st->geoInfo_wellDoublet_well2_liquidBC->getGeoObj())) + ShiftInNodeVector);
 	       		 liquidBC_mesh_node_values.resize(liquidBC_mesh_nodes.size(), 1.);
 	       	 }
 	       	 else if(st->geoInfo_wellDoublet_well2_liquidBC->getGeoType() == GEOLIB::POLYLINE)
 	       	 {
 	       		 m_msh->GetNODOnPLY(static_cast<const GEOLIB::Polyline*>(
 	       				 st->geoInfo_wellDoublet_well2_liquidBC->getGeoObj()), liquidBC_mesh_nodes, true);
 	       		 liquidBC_mesh_node_values.resize(liquidBC_mesh_nodes.size(), 1.);
 	       		 FiniteElement::EdgeIntegration(m_msh, liquidBC_mesh_nodes, liquidBC_mesh_node_values,
 	       				 st->getProcessDistributionType(), st->getProcessPrimaryVariable(), 
 	       				 st->ignore_axisymmetry, st->isPressureBoundaryCondition(), st->scaling_mode);
 	       		 total_value = std::accumulate(liquidBC_mesh_node_values.begin(), liquidBC_mesh_node_values.end(), 0.);
 	       		 for(auto& value: liquidBC_mesh_node_values) value /= total_value;
 	       	 }

 	       	 for(long i=0; i < liquidBC_mesh_nodes.size(); ++i)
 	       	 {
 	       		 CNodeValue *nod_val_liquid_well (new CNodeValue());
 	       		 nod_val_liquid_well->msh_node_number = liquidBC_mesh_nodes[i];
 	       		 nod_val_liquid_well->CurveIndex = st->CurveIndex;
 	       		 nod_val_liquid_well->geo_node_number = nod_val_liquid_well->msh_node_number - ShiftInNodeVector;
 	       		 nod_val_liquid_well->node_value = -st->geo_node_value * liquidBC_mesh_node_values[i];  // !!!!!
 	       		 nod_val_liquid_well->tim_type_name = st->tim_type_name;

 	       		 pcs_liquid->st_node_value.push_back(nod_val_liquid_well);
 	       		 pcs_liquid->st_node.push_back(st);
 	       	 }

 	         /*
 	         // warm well 1
 	         CNodeValue *nod_val_liquid_well1 (new CNodeValue());
 	         nod_val_liquid_well1->msh_node_number = m_msh->GetNODOnPNT(static_cast<const GEOLIB::Point*>(
 	       		  st->geoInfo_wellDoublet_well1_liquidBC->getGeoObj())) + ShiftInNodeVector;
 	         nod_val_liquid_well1->CurveIndex = st->CurveIndex;
 	         nod_val_liquid_well1->geo_node_number = nod_val_liquid_well1->msh_node_number - ShiftInNodeVector;
 	         nod_val_liquid_well1->node_value = st->geo_node_value;
 	         nod_val_liquid_well1->tim_type_name = st->tim_type_name;

 	         pcs_liquid->st_node_value.push_back(nod_val_liquid_well1);
 	         pcs_liquid->st_node.push_back(st);
 	         // cold well 2 - copy from lines above except node_value is negative
 	         CNodeValue *nod_val_liquid_well2 (new CNodeValue());
 	         nod_val_liquid_well2->msh_node_number = m_msh->GetNODOnPNT(static_cast<const GEOLIB::Point*>(
 	       		  st->geoInfo_wellDoublet_well2_liquidBC->getGeoObj())) + ShiftInNodeVector;
 	         nod_val_liquid_well2->CurveIndex = st->CurveIndex;
 	         nod_val_liquid_well2->geo_node_number = nod_val_liquid_well2->msh_node_number - ShiftInNodeVector;
 	         nod_val_liquid_well2->node_value = -st->geo_node_value;  // !!!!!
 	         nod_val_liquid_well2->tim_type_name = st->tim_type_name;

 	         pcs_liquid->st_node_value.push_back(nod_val_liquid_well2);
 	         pcs_liquid->st_node.push_back(st);
 	         */
 	  }

 	  if(st->ogs_contraflow != nullptr)  // JOD 2019-31-07
 	  {
 	          st->ogs_contraflow->add_node(nod_val->msh_node_number);

 	          GEOLIB::Point pnt(
 	       		   static_cast<const GEOLIB::Point*>(st->getGeoObj())->getData()[0],
 	       		   static_cast<const GEOLIB::Point*>(st->getGeoObj())->getData()[1],
 	       		   static_cast<const GEOLIB::Point*>(st->getGeoObj())->getData()[2]);

 	          std::vector<contra::SegmentData> segment_data_vec = st->ogs_contraflow->get_segment_data_vec();
 	          double z = pnt.getData()[2];
 	          for(int i=0; i < segment_data_vec.size(); ++i)
 	          {
 	       	   const int N = segment_data_vec[i].N;
 	       	   const double dz = segment_data_vec[i].L / N;

 	       	   for(int j=0; j<N; ++j)
 	       	   {
 	       		   z -= dz;
 	       		   pnt = GEOLIB::Point(
 	       				   static_cast<const GEOLIB::Point*>(st->getGeoObj())->getData()[0],
 	       				   static_cast<const GEOLIB::Point*>(st->getGeoObj())->getData()[1],
 	       				   z);
 	       		   st->ogs_contraflow->add_node(m_msh->GetNODOnPNT(static_cast<const GEOLIB::Point*>(&pnt)));
 	       	   }
 	          }
 	  }
	
	if(st->borehole_mode == 0)
	{
		double factor, radius_e;
	
		CalculatePeaceman(st, pcs, nod_val->geo_node_number,
			 pcs->m_msh->nod_vector[nod_val->geo_node_number]->getConnectedElementIDs()
				, factor, radius_e);
	
		factor /= std::log(radius_e / st->borehole_data.radius);
		nod_val->node_value *= factor;
	 	// to use advective HEAT_TRANSPORT with conductive (peaceman) LIQUID_FLOW 
	 	if(pcs->ST_factor_kept.find(nod_val->geo_node_number) != pcs->ST_factor_kept.end())
	 		pcs->ST_factor_kept[nod_val->geo_node_number] = pcs->ST_factor_kept[nod_val->geo_node_number] + factor;
	 	else
	 		pcs->ST_factor_kept[nod_val->geo_node_number] = factor;
	
	}
	
	pcs->st_node_value.push_back(nod_val);         //WW
   	pcs->st_node.push_back(st);                 //WW

}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 11/2007 JOD Implementation
 **************************************************************************/
void CSourceTermGroup::SetLIN(CRFProcess* m_pcs, CSourceTerm* m_st,
const int ShiftInNodeVector)
{
   (void)m_pcs;
   (void)m_st;
   (void)ShiftInNodeVector;
   /*OK411
    long number_of_nodes;
    vector<long>lin_nod_vector;
    vector<double>lin_nod_val_vector;
    CGLLine* m_lin = NULL;
    CGLPolyline* m_ply = NULL;
    long *nodes = NULL;
    m_lin = m_lin->GEOGetLine(m_st->geo_id);

    if(m_lin){
    double* coordinates;
   m_ply = new CGLPolyline;
   m_ply->point_vector.push_back(m_lin->m_point1);
   m_ply->point_vector.push_back(m_lin->m_point2);
   nodes = MSHGetNodesClose(&number_of_nodes, m_ply);//CC
   lin_nod_val_vector.resize(number_of_nodes);
   for(long i = 0; i < number_of_nodes; i++){
   lin_nod_val_vector[i] =  m_st->geo_node_value / number_of_nodes;
   coordinates = new double[3];
   coordinates[0] = GetNodeX(nodes[i]);
   coordinates[1] = GetNodeY(nodes[i]);
   coordinates[2] = GetNodeZ(nodes[i]);
   m_lin->nodes_coor_vector.push_back(coordinates);
   }
   //InterpolationAlongPolyline(m_polyline,node_value_vector);
   for(long i=0; i < number_of_nodes; i++){
   CNodeValue* m_nod_val = NULL;
   m_nod_val = new CNodeValue();
   m_nod_val->msh_node_number = -1;
   m_nod_val->msh_node_number = nodes[i]+ShiftInNodeVector;
   m_nod_val->geo_node_number = nodes[i];
   m_nod_val->node_value = lin_nod_val_vector[i];
   m_nod_val->CurveIndex = m_st->CurveIndex;
   //WW        group_vector.push_back(m_node_value);
   //WW        st_group_vector.push_back(m_st); //OK
   m_pcs->st_node_value.push_back(m_nod_val);  //WW
   m_pcs->st_node.push_back(m_st); //WW
   }
   lin_nod_val_vector.clear();
   m_ply->point_vector.clear();
   delete m_ply;
   }
   else
   cout << "Warning - CSourceTermGroup::Set: LIN not found" << '\n';
   */
}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing
 07/2005 OK Implementation based on CSourceTermGroup::Set
 modified
 07/2010 TF substituted GEOGetPLYByName
 **************************************************************************/
/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 11/2007 JOD Implementation
 **************************************************************************/
void CSourceTermGroup::SetDMN(CSourceTerm *m_st, const int ShiftInNodeVector)
{
   std::vector<long> dmn_nod_vector;
   std::vector<double> dmn_nod_val_vector;
   std::vector<long> dmn_nod_vector_cond;
   std::vector<double> dmn_nod_val_vector_cond;  // not used

   GEOGetNodesInMaterialDomain(m_msh, m_st->analytical_material_group,
      dmn_nod_vector, false);
   size_t number_of_nodes (dmn_nod_vector.size());
   dmn_nod_val_vector.resize(number_of_nodes);

   for (size_t i = 0; i < number_of_nodes; i++)
      dmn_nod_val_vector[i] = m_st->geo_node_value;

   m_st->SetNodeValues(dmn_nod_vector, dmn_nod_vector_cond,
      dmn_nod_val_vector, dmn_nod_val_vector_cond, ShiftInNodeVector);

}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 02/2008 JOD Implementation
 **************************************************************************/
void CSourceTermGroup::SetCOL(CSourceTerm *m_st, const int ShiftInNodeVector)
{
   long number_of_nodes;
   std::vector<long> col_nod_vector;
   std::vector<double> col_nod_val_vector;
   std::vector<long> col_nod_vector_cond;
   std::vector<double> col_nod_val_vector_cond;  // nod used

   long i = 0;
   if (m_st->geo_name == "BOTTOM")
      i = m_msh->getNumberOfMeshLayers() - 1;

   while (i < (long) m_msh->nod_vector.size())
   {
      col_nod_vector.push_back(i);
      i += m_msh->getNumberOfMeshLayers();
   }
   number_of_nodes = (long) col_nod_vector.size();
   col_nod_val_vector.resize(number_of_nodes);

   for (long i = 0; i < number_of_nodes; i++)
      col_nod_val_vector[i] = 1;

   m_st->SetSurfaceNodeVectorConditional(col_nod_vector, col_nod_vector_cond);

   m_st->SetNodeValues(col_nod_vector, col_nod_vector_cond,
      col_nod_val_vector, col_nod_val_vector_cond, ShiftInNodeVector);

}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 11/2007 JOD Implementation
 09/2010 WW  For the case that nodes are directly given
 **************************************************************************/
void CSourceTermGroup::SetSFC(CSourceTerm* m_st, const int ShiftInNodeVector)
{
   std::vector<long> sfc_nod_vector;
   std::vector<std::size_t> sfc_node_ids;
   std::vector<long> sfc_nod_vector_cond;
   std::vector<double> sfc_nod_val_vector;
   std::vector<double> sfc_nod_val_vector_cond;  // not used
   Surface* m_sfc = NULL;

   m_sfc = GEOGetSFCByName(m_st->geo_name);       //CC

   if (m_sfc)
   {
      std::cout << "Surface " << m_st->geo_name << '\n';
      SetSurfaceNodeVector( m_st->geo_name, m_st->getGeoObj(), sfc_node_ids);

      sfc_nod_vector.insert(sfc_nod_vector.begin(),
         sfc_node_ids.begin(), sfc_node_ids.end());
      if (m_st->isCoupled())
         m_st->SetSurfaceNodeVectorConditional(sfc_nod_vector,
            sfc_nod_vector_cond);
// CB JOD MERGE //

	  if (m_st->isConnectedGeometry() &&   // JOD 2/2015
			  m_st->getConnectedGeometryCouplingType() != 2) // not borehole with given primary variable
		  m_st->SetSurfaceNodeVectorConnected(sfc_nod_vector, sfc_nod_vector_cond);
	   if(m_st->hasThreshold()) // JOD 2018-03-7  - copyed from SetPnt
	   {  // only point supported
		   m_st->msh_node_number_threshold = m_msh->GetNODOnPNT(
						static_cast<const GEOLIB::Point*>(m_st->geoInfo_threshold->getGeoObj()));
	   }

	   if(m_st->calculatedFromStorageRate()) // JOD 2018-03-7  - copyed from SetPnt
	   {  // only point supported
		   if (m_st->geoInfo_storageRateInlet->getGeoType() == GEOLIB::POINT)
		   {
			   m_st->storageRate.inlet_msh_node_numbers.push_back(m_msh->GetNODOnPNT(
				static_cast<const GEOLIB::Point*>(m_st->geoInfo_storageRateInlet->getGeoObj())));
			   m_st->storageRate.outlet_msh_node_numbers.push_back(m_msh->GetNODOnPNT(
					static_cast<const GEOLIB::Point*>(m_st->geoInfo_storageRateOutlet->getGeoObj())));
			   m_st->storageRate.inlet_msh_node_areas.push_back(1.);
			   m_st->storageRate.outlet_msh_node_areas.push_back(1.);
			   m_st->storageRate.inlet_totalArea = 1.;
			   m_st->storageRate.outlet_totalArea = 1.;
		   }
		   if (m_st->geoInfo_storageRateInlet->getGeoType() == GEOLIB::SURFACE)
		   {
			   // inlet surface
			   sfc_node_ids.clear();
			   SetSurfaceNodeVector(m_st->storageRate.inlet_geometry_name, m_st->geoInfo_storageRateInlet->getGeoObj(), sfc_node_ids);
			   m_st->storageRate.inlet_msh_node_numbers.insert(m_st->storageRate.inlet_msh_node_numbers.begin(),
			         sfc_node_ids.begin(), sfc_node_ids.end());

			      sfc_nod_vector.insert(sfc_nod_vector.begin(),
			         sfc_node_ids.begin(), sfc_node_ids.end());
			   SetSurfaceNodeValueVector(m_st, GEOGetSFCByName(m_st->storageRate.inlet_geometry_name),
					   m_st->storageRate.inlet_msh_node_numbers, m_st->storageRate.inlet_msh_node_areas);

			   m_st->storageRate.inlet_totalArea = 0.;
			   for(int i=0; i < m_st->storageRate.inlet_msh_node_areas.size(); i++)
				   m_st->storageRate.inlet_totalArea += m_st->storageRate.inlet_msh_node_areas[i];

			   // outlet surface
			   sfc_node_ids.clear();
			   SetSurfaceNodeVector(m_st->storageRate.outlet_geometry_name, m_st->geoInfo_storageRateOutlet->getGeoObj(), sfc_node_ids);
			   m_st->storageRate.outlet_msh_node_numbers.insert(m_st->storageRate.outlet_msh_node_numbers.begin(),
			         sfc_node_ids.begin(), sfc_node_ids.end());

			   SetSurfaceNodeValueVector(m_st, GEOGetSFCByName(m_st->storageRate.outlet_geometry_name),
					   m_st->storageRate.outlet_msh_node_numbers, m_st->storageRate.outlet_msh_node_areas);

			   m_st->storageRate.outlet_totalArea = 0.;
			   for(int i=0; i < m_st->storageRate.outlet_msh_node_areas.size(); i++)
				   m_st->storageRate.outlet_totalArea += m_st->storageRate.outlet_msh_node_areas[i];
		   }
	   }

	//		m_st->SetDISType();
	SetSurfaceNodeValueVector(m_st, m_sfc, sfc_nod_vector, sfc_nod_val_vector);

	if(m_st->verbosity)
	{
		std::cout << "\t" << sfc_nod_vector.size() << " nodes with total value of " << std::accumulate(sfc_nod_val_vector.begin(), sfc_nod_val_vector.end(), 0.) << '\n';
		for(int i=0; i< sfc_nod_vector.size(); i++)
		{
			std::cout << "\t\t" << i << ": " << sfc_nod_vector[i] << " with value " << sfc_nod_val_vector[i];
		  	if (m_st->isConnectedGeometry() &&
				  m_st->getConnectedGeometryCouplingType() != 2) // not borehole with given primary variable
				std::cout << " - connected to " << sfc_nod_vector_cond[i];
			std::cout << '\n';
		}
	}


	  if (m_st->everyoneWithEveryone)  // JOD 8/2015 - modified 2022-02-23
	  {
		  std::cout << "\tConnect every node with every node of other surface\n";
		  std::vector<long>::iterator pos;
		  std::vector<double> sfc_nod_val_vector_cond_original;

		  SetSurfaceNodeValueVector(m_st, m_sfc, sfc_nod_vector_cond,
			  sfc_nod_val_vector_cond_original);
		  const double total_val_cond_original = std::accumulate(sfc_nod_val_vector_cond_original.begin(), sfc_nod_val_vector_cond_original.end(), 0.);

	  	  if (m_st->verbosity && m_st->isConnectedGeometry() &&
			  m_st->getConnectedGeometryCouplingType() != 2) // not borehole with given primary variable
			std::cout << "\t\t\t" << sfc_nod_vector_cond.size() << " connected nodes with total value of " << total_val_cond_original << '\n';

		  const int nod_vector_size_original = (int)sfc_nod_vector.size();
		  const int nod_vector_cond_size_original = (int)sfc_nod_vector_cond.size();

		  std::vector<double> sfc_nod_val_vector_original(sfc_nod_val_vector);

		  // examples
		  // nod_vector: 0, 1
		  // nod_vector_cond: 2, 3, 4
		  for (int i = 0; i < nod_vector_size_original; i++)  // extend nod_vector
		  {  // 0, 1 -> 0,0,0,1,1,1
			  for (int j = 1; j < nod_vector_cond_size_original; j++)
			  {
				  pos = sfc_nod_vector.begin() + i * nod_vector_cond_size_original + j;
				  sfc_nod_vector.insert(pos, sfc_nod_vector[i * nod_vector_cond_size_original + j - 1]);
			  }
		  }

		  for (int i = 1; i < nod_vector_size_original; i++)  // extend nod_vector_cond
		  {  // 2, 3, 4 -> 2, 3, 4, 2, 3, 4 
			  for (int j = 0; j < nod_vector_cond_size_original; j++)
			  {
				  sfc_nod_vector_cond.push_back(sfc_nod_vector_cond[j]);
			  }
		  }

		  // extend nod_val_vector
		  sfc_nod_val_vector.resize(sfc_nod_vector.size());

		  for (int i = 0; i < nod_vector_size_original; i++)
		  {
			  for (int j = 0; j < nod_vector_cond_size_original; j++)
			  {
				  sfc_nod_val_vector[i*nod_vector_cond_size_original + j] = sfc_nod_val_vector_original[i] * sfc_nod_val_vector_cond_original[j] / (total_val_cond_original);
			  }
		  }

		  if(m_st->verbosity)
		  {
		  	std::cout << "\tNow " << sfc_nod_vector.size() << " nodes with total value of " << std::accumulate(sfc_nod_val_vector.begin(), sfc_nod_val_vector.end(), 0.) << '\n';
		  	if(m_st->verbosity > 1)
		  		for (int i = 0; i < sfc_nod_vector.size(); i++)
		  		{
					  std::cout << "\t\t" << i << ": " << sfc_nod_vector[i] << " connected to " << sfc_nod_vector_cond[i] << " with value "<< sfc_nod_val_vector[i] << '\n';
		  		}
		  }
	  }  // end everyone with everyone

	  if (m_st->distribute_volume_flux)   // 5.3.07 JOD
		  DistributeVolumeFlux(m_st, sfc_nod_vector, sfc_nod_val_vector);

	  m_st->st_node_ids.clear();
	  m_st->st_node_ids.resize(sfc_nod_vector.size());
	  m_st->st_node_ids = sfc_nod_vector;

	  if (m_st->isConstrainedST())
	  {
		  for (std::size_t i(0); i < m_st->st_node_ids.size(); i++)
		  {
			  m_st->_constrainedSTNodesIndices.push_back(-1);
			  for (std::size_t i(0); i < m_st->getNumberOfConstrainedSTs(); i++)
				  m_st->pushBackConstrainedSTNode(i, false);
		  }
	  }

      m_st->SetNodeValues(sfc_nod_vector, sfc_nod_vector_cond,
         sfc_nod_val_vector, sfc_nod_val_vector_cond, ShiftInNodeVector);

   }                                              // end surface
}


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 11/2007 JOD Implementation
 **************************************************************************/
void CSourceTerm::SetNOD()
{

   std::vector<long> nod_vector;
   std::vector<long> nod_vector_cond;
   std::vector<double> nod_val_vector;
   std::vector<double> nod_val_vector_cond;
   int ShiftInNodeVector;

   nod_vector.push_back(msh_node_number);
   nod_vector_cond.push_back(msh_node_number);
   nod_val_vector.push_back(geo_node_value);
   nod_val_vector_cond.push_back(geo_node_value);

   /*nod_vector[0] = msh_node_number;
    nod_vector_cond[0] = msh_node_number;
    nod_val_vector[0] =geo_node_value;*/
   ShiftInNodeVector = 0;

   SetNodeValues(nod_vector, nod_vector_cond, nod_val_vector, nod_val_vector_cond,
      ShiftInNodeVector);

}


/**************************************************************************
 MSHLib-Method:
 Task:
 Programing:
 11/2007 JOD
 last modification:
 **************************************************************************/
//void CSourceTermGroup::SetPolylineNodeVector(CGLPolyline* m_ply,
//std::vector<long>&ply_nod_vector)
//{
//   if (m_ply->getType() == 100)                   // WW
//      m_msh->GetNodesOnArc(m_ply, ply_nod_vector);
//   else if (m_ply->getType() == 3)                // JOD
//      m_msh->GetNODOnPLY_XY(m_ply, ply_nod_vector);
//   else
//      m_msh->GetNODOnPLY(m_ply, ply_nod_vector);
//}


/**************************************************************************
 MSHLib-Method:
 Task:
 Programing:
 11/2007 JOD
 last modification:
 11/2021 JOD
 **************************************************************************/
void CSourceTermGroup::SetSurfaceNodeVector(const std::string geo_name, const GEOLIB::GeoObject* geo_obj,
		std::vector<size_t>&sfc_nod_vector)
{
	Surface* m_surface = GEOGetSFCByName(geo_name);
	if (m_surface->TIN)
	{
		std::vector<long> nodes_vector;
		m_msh->GetNODOnSFC_TIN(m_surface, nodes_vector);
		sfc_nod_vector.resize(nodes_vector.size());
		for(int i=0; i< nodes_vector.size(); ++i)
			sfc_nod_vector[i] = nodes_vector[i];
		std::cout << "\tTIN " << geo_name << ": " << nodes_vector.size() << " nodes" << std::endl;
	}
	else
	{
		const bool for_source = true;
		GEOLIB::Surface const* sfc(static_cast<const GEOLIB::Surface*> (geo_obj));
		m_msh->GetNODOnSFC(sfc, sfc_nod_vector, for_source);
	}

}

/**************************************************************************
 MSHLib-Method:
 Task:
 Programing:
 11/2007 JOD
 12/2012 JOD Extension to TWO_PHASE_FLOW 5.3.07
 last modification:
 **************************************************************************/
void CSourceTermGroup::SetPolylineNodeVectorConditional(CSourceTerm* st,
		std::vector<long>& ply_nod_vector, std::vector<long>& ply_nod_vector_cond)
{
   size_t assembled_mesh_node, number_of_nodes;

   if (st->node_averaging)
   {
      if (m_msh_cond)
      {
         if (pcs_type_name == "RICHARDS_FLOW" || pcs_type_name == "MULTI_PHASE_FLOW")
         {
            //				m_msh_cond->GetNODOnPLY(m_ply, ply_nod_vector_cond);
            m_msh_cond->GetNODOnPLY(static_cast<const GEOLIB::Polyline*>(st->getGeoObj()), ply_nod_vector_cond);
            number_of_nodes = ply_nod_vector_cond.size();
            assembled_mesh_node = 0;// ply_nod_vector[0]; // JOD  carefull !!! ply sometimes fails
            ply_nod_vector.resize(number_of_nodes);
            for (size_t i = 0; i < number_of_nodes; i++)
               ply_nod_vector[i] = assembled_mesh_node;
         }                                        // end richards / multi_phase
         else if (pcs_type_name == "OVERLAND_FLOW" || pcs_type_name == "GROUNDWATER_FLOW")                // JOD 4.10.01
         {
            number_of_nodes = ply_nod_vector.size();
            //				m_msh_cond->GetNODOnPLY(m_ply, ply_nod_vector_cond);
            m_msh_cond->GetNODOnPLY(static_cast<const GEOLIB::Polyline*>(st->getGeoObj()), ply_nod_vector_cond);
            assembled_mesh_node = 0; // ply_nod_vector_cond[0];  // JOD  carefull !!! ply sometimes fails
            ply_nod_vector_cond.resize(number_of_nodes);
            for (size_t i = 0; i < number_of_nodes; i++)
               ply_nod_vector_cond[i] = assembled_mesh_node;
         }                                        // end overland, groundwater
         else
            std::cout
               << "Warning in CSourceTermGroup::SetPolylineNodeVectorConditional - no area assembly for this process"
               << "\n";
      }                                           // end mesh_cond
      else
         std::cout << "Warning in CSourceTermGroup::SetPLY - no MSH_COND data"
            << "\n";
   }                                              // end area_assembly
   else
   {
      number_of_nodes = ply_nod_vector.size();
      ply_nod_vector_cond.resize(number_of_nodes);
      st->SetNOD2MSHNOD(ply_nod_vector, ply_nod_vector_cond);
   }                                              // end !area_assembly
}


/**************************************************************************
 MSHLib-Method:
 Task:
 Programing:
 11/2007 JOD
 last modification:
 **************************************************************************/
void CSourceTerm::InterpolatePolylineNodeValueVector(CGLPolyline* m_ply,
		std::vector<double>& Distribed, std::vector<double>& ply_nod_vector)
{
	for (long k = 0; k < (long) DistribedBC.size(); k++) {
		for (long l = 0; l < (long) m_ply->point_vector.size(); l++) {
			if (PointsHaveDistribedBC[k] == m_ply->point_vector[l]->id) {
				if (fabs(DistribedBC[k]) < MKleinsteZahl)
					DistribedBC[k] = 1.0e-20;
				m_ply->point_vector[l]->setPropert(Distribed[k]);
				break;
			}
		}
	}

	InterpolationAlongPolyline(m_ply, ply_nod_vector);
}

void CSourceTerm::InterpolatePolylineNodeValueVector(
		std::vector<double> const & nodes_as_interpol_points,
		std::vector<double>& node_values) const
{
	std::vector<double> interpolation_points;
	std::vector<double> interpolation_values;

	GEOLIB::Polyline const * ply (dynamic_cast<GEOLIB::Polyline const*>(this->getGeoObj()));

	for (size_t i(0); i < DistribedBC.size(); i++) {
		for (size_t j = 0; j < ply->getNumberOfPoints(); j++) {
			if ((size_t)(PointsHaveDistribedBC[i]) == ply->getPointID(j)) {
				interpolation_points.push_back (ply->getLength(j));
				if (fabs(DistribedBC[i]) < MKleinsteZahl)
					interpolation_values.push_back (1.0e-20);
				else
					interpolation_values.push_back (DistribedBC[i]);
				break;
			}
		}
	}

	MathLib::PiecewiseLinearInterpolation (interpolation_points, interpolation_values, nodes_as_interpol_points, node_values);
}



// 09/2010 TF
void CSourceTermGroup::SetPolylineNodeValueVector(CSourceTerm* st,
		std::vector<long> const & ply_nod_vector,
		std::vector<long> const & ply_nod_vector_cond,
		std::vector<double>& ply_nod_val_vector) const
{
	size_t number_of_nodes(ply_nod_vector.size());
	ply_nod_val_vector.resize(number_of_nodes);

	FiniteElement::DistributionType distype(st->getProcessDistributionType());

	// linear
	if (distype == FiniteElement::LINEAR || distype == FiniteElement::LINEAR_NEUMANN) {
		// fetch data for the linear interpolation
		GEOLIB::Polyline const* polyline (dynamic_cast<GEOLIB::Polyline const*>(st->getGeoObj()));
		if (polyline) {
			std::vector<double> nodes_as_interpol_points;
			m_msh->getPointsForInterpolationAlongPolyline (polyline, nodes_as_interpol_points);
			st->InterpolatePolylineNodeValueVector(nodes_as_interpol_points, ply_nod_val_vector);
		}
	} else if (distype == FiniteElement::SYSTEM_DEPENDENT) { //System Dependented YD
		CRFProcess* m_pcs(PCSGet(pcs_type_name));
		m_pcs->compute_domain_face_normal = true; //WW
		long no_face = (long) m_msh->face_vector.size();
		for (long i = 0; i < no_face; i++) {
			int node_on_line = 0;
			int no_vertex = m_msh->face_vector[i]->GetVertexNumber();
			for (long jj = 0; jj < no_vertex; jj++) {
				for (size_t kk = 0; kk < number_of_nodes; kk++) {
					if (ply_nod_vector[kk] == (m_msh->face_vector[i]->GetNodeIndex(jj)))
						node_on_line++;
				} // end nodes
			} // end vertices
			if (node_on_line == 2) st->element_st_vector.push_back(
					m_msh->face_vector[i]->GetOwner()->GetIndex());
		} // end faces
	} // end system dependent
   else if (distype == FiniteElement::FUNCTION) // 25.08.2011. WW
   {
      for (size_t i = 0; i < number_of_nodes; i++)
      {
         double const*const pnt (m_msh->nod_vector[ply_nod_vector[i]]->getData());
         ply_nod_val_vector[i] = st->dis_linear_f->getValue(pnt[0], pnt[1], pnt[2]);
      }
   }
	else //WW
	{
		for (size_t i = 0; i < number_of_nodes; i++) {
			ply_nod_val_vector[i] = st->geo_node_value;
			if (st->getProcessDistributionType() == FiniteElement::CONSTANT_GEO)
				ply_nod_val_vector[i] = st->geo_node_value / (double) number_of_nodes;
		}
	}
	if (distype == FiniteElement::CONSTANT_NEUMANN
			|| distype == FiniteElement::LINEAR_NEUMANN
			|| distype == FiniteElement::GREEN_AMPT
			|| distype==FiniteElement::RECHARGE)	//MW
	{
		if (m_msh->GetMaxElementDim() == 1 || st->assign_to_element_edge) // 1D  //WW MB
			FiniteElement::DomainIntegration(PCSGet(pcs_type_name), ply_nod_vector,
					ply_nod_val_vector);
		else FiniteElement::EdgeIntegration(m_msh, ply_nod_vector, ply_nod_val_vector,
				st->getProcessDistributionType(), st->getProcessPrimaryVariable(), 
				st->ignore_axisymmetry, st->isPressureBoundaryCondition(), st->scaling_mode);
	}

	if (distype == FiniteElement::CRITICALDEPTH
			|| distype == FiniteElement::NORMALDEPTH
			|| distype == FiniteElement::ANALYTICAL) {
		st->node_value_vectorArea.resize(number_of_nodes);
		for (size_t i = 0; i < number_of_nodes; i++)
			st->node_value_vectorArea[i] = 1.0; //Element width !
		FiniteElement::EdgeIntegration(m_msh, ply_nod_vector,
				st->node_value_vectorArea,
				st->getProcessDistributionType(), st->getProcessPrimaryVariable(), 
				st->ignore_axisymmetry, st->isPressureBoundaryCondition());
	}

	if (distype == FiniteElement::TRANSFER_SURROUNDING) { //TN - Belegung mit Fl�chenelementen
		st->node_value_vectorArea.resize(number_of_nodes);
		for (size_t i = 0; i < number_of_nodes; i++){
			st->node_value_vectorArea[i] = 1.0;
			ply_nod_val_vector[i] = 0.0;
		}

		if (m_msh->GetMaxElementDim() == 1) // 1D  //WW MB
			FiniteElement::DomainIntegration(PCSGet(pcs_type_name), ply_nod_vector,
					 st->node_value_vectorArea);
		else
			FiniteElement::EdgeIntegration(m_msh, ply_nod_vector, st->node_value_vectorArea,
					st->getProcessDistributionType(), st->getProcessPrimaryVariable(), 
					st->ignore_axisymmetry, st->isPressureBoundaryCondition());

		
		//CNode * a_node; //Fl�che wird hier im Knoten als patch area abgelegt
		//for (size_t i = 0; i < number_of_nodes; i++) 
		//{ 
		//	a_node = m_msh->nod_vector[ply_nod_vector[i]]; 
		//	a_node->patch_area  = st->node_value_vectorArea[i];
		//}
	}

	if (st->isCoupled() && st->node_averaging)
		AreaAssembly(st, ply_nod_vector_cond, ply_nod_val_vector);

	if (distype==FiniteElement::RECHARGE)	//MW
	{
		CRFProcess* m_pcs = PCSGet(pcs_type_name);
		MshEditor::sortNodesLexicographically(m_pcs->m_msh);

		st->setProcessDistributionType (st->getProcessDistributionType());
	}
}

/**************************************************************************
 MSHLib-Method:
 Task:
 Programing:
 11/2007 JOD
 12/2012 JOD Extension to TWO_PHASE_FLOW
 04/2020 JOD deactivated
 last modification:
 **************************************************************************/
void CSourceTermGroup::AreaAssembly(const CSourceTerm* const st,
const std::vector<long>& ply_nod_vector_cond,
std::vector<double>& ply_nod_val_vector) const
{
   if (pcs_type_name == "RICHARDS_FLOW" || pcs_type_name == "MULTI_PHASE_FLOW")
   {
      /*if (m_msh_cond->GetMaxElementDim() == 1)    // 1D  //WW MB
         DomainIntegration(m_msh_cond, ply_nod_vector_cond,
            ply_nod_val_vector);
      else
         st->EdgeIntegration(m_msh_cond, ply_nod_vector_cond,
            ply_nod_val_vector);
	    */
      double sum_node_value = 0;
      for (size_t i = 0; i < ply_nod_val_vector.size(); i++)
         sum_node_value += ply_nod_val_vector[i];
      for (size_t  i = 0; i < ply_nod_val_vector.size(); i++)
         ply_nod_val_vector[i] /= sum_node_value;
   }
}


/**************************************************************************
 MSHLib-Method:
 Task:
 Programing:
 11/2007 JOD
 last modification:
 **************************************************************************/
void CSourceTermGroup::SetSurfaceNodeValueVector(CSourceTerm* st,
		Surface* m_sfc, std::vector<long> const &sfc_nod_vector,
		std::vector<double>&sfc_nod_val_vector)
{
   // CRFProcess* m_pcs = NULL;
   // m_pcs = PCSGet(pcs_type_name);
   if (m_sfc == NULL)
	   throw std::runtime_error("Surface unknown");

   long number_of_nodes = (long) sfc_nod_vector.size();
   sfc_nod_val_vector.resize(number_of_nodes);

   for (long i = 0; i < number_of_nodes; i++){
       sfc_nod_val_vector[i] = st->geo_node_value;
       if (st->getProcessDistributionType() == FiniteElement::CONSTANT_GEO) // added by CB
         sfc_nod_val_vector[i] /= (double)number_of_nodes;
   }


   // KR & TF - case not used
   //	if (m_st->dis_type == 12) //To do. 4.10.06
   //		for (long i = 0; i < number_of_nodes; i++)
   //			sfc_nod_val_vector[i] = m_st->geo_node_value
   //					/ (double) number_of_nodes;

   //	if (st->dis_type == 2 || st->dis_type == 4) { // Piecewise linear distributed, polygon-wise WW
   if (st->getProcessDistributionType() == FiniteElement::LINEAR || st->getProcessDistributionType() == FiniteElement::LINEAR_NEUMANN)
   {
      CGLPolyline* m_ply = NULL;
      std::vector<CGLPolyline*>::iterator p =
         m_sfc->polyline_of_surface_vector.begin();
      p = m_sfc->polyline_of_surface_vector.begin();
      while (p != m_sfc->polyline_of_surface_vector.end())
      {
         m_ply = *p;
         long nPointsPly = (long) m_ply->point_vector.size();
         if (m_ply->point_vector.front() == m_ply->point_vector.back())
            nPointsPly -= 1;

         for (long k = 0; k < (long) st->DistribedBC.size(); k++)
         {
            for (long l = 0; l < nPointsPly; l++)
            {
               if (st->PointsHaveDistribedBC[k]
                  == m_ply->point_vector[l]->id)
               {
                  if (fabs(st->DistribedBC[k]) < MKleinsteZahl)
                     st->DistribedBC[k] = 1.0e-20;
                  m_ply->point_vector[l]->setPropert (st->DistribedBC[k]);
                  break;
               }                                  // end l
            }                                     // end k
         }                                        // end polyline
         // InterpolationAlongPolyline(m_polyline, node_value_vector);
         p++;
      }                                           // end while
   }                                              // end linear

   // neumann, Green-Ampt, Philip
   //	if (st->dis_type == 3 || st->dis_type == 4 || st->dis_type == 10
   //				|| st->dis_type == 11) {
   /*|| st->getProcessDistributionType() == PHILIP */
   if (st->getProcessDistributionType() == FiniteElement::CONSTANT_NEUMANN
		   || st->getProcessDistributionType() == FiniteElement::LINEAR_NEUMANN
		   || st->getProcessDistributionType() == FiniteElement::GREEN_AMPT
		   || st->getProcessDistributionType() == FiniteElement::RECHARGE)
   {
      if (m_msh->GetMaxElementDim() == 2)         // For all meshes with 1-D or 2-D elements
    	  FiniteElement::DomainIntegration(PCSGet(pcs_type_name), sfc_nod_vector, sfc_nod_val_vector);
      else if (m_msh->GetMaxElementDim() == 3)    // For all meshes with 3-D elements
    	  FiniteElement::FaceIntegration(m_msh, sfc_nod_vector, sfc_nod_val_vector, m_sfc,
        		 st->getProcessDistributionType(), st->_pcs->m_num->ele_gauss_points);
   }                                              // end neumann
  else if (st->getProcessDistributionType() == FiniteElement::FUNCTION) // 25.08.2011. WW
   {
      for (size_t j = 0; j < sfc_nod_vector.size(); j++)
      {
         double const*const pnt (m_msh->nod_vector[sfc_nod_vector[j]]->getData());
         sfc_nod_val_vector[j] = st->dis_linear_f->getValue(pnt[0], pnt[1], pnt[2]);
      }
   }

}

/**************************************************************************
 FEMLib-Method: Distributes source term [m^3/s] uniformly on SFC, PLY 
  use GEO_TYPE CONSTANT_NEUMANN 
 12/2012 5.3.07 JOD Implementation
 **************************************************************************/
void CSourceTermGroup::DistributeVolumeFlux(CSourceTerm* st, std::vector<long> const & nod_vector, 
		                      std::vector<double>& nod_val_vector)
{

	double area;
    std::vector<double> nod_val_vector_area;
	    
	area = 0;    
    st->geo_node_value = 1;

	if(st->getGeoType () == GEOLIB::POLYLINE)
	  SetPolylineNodeValueVector(st, nod_vector, nod_vector, nod_val_vector_area);
	else if(st->getGeoType () == GEOLIB::SURFACE)
	{
	   Surface* m_sfc = NULL;
       m_sfc = GEOGetSFCByName(st->geo_name); 	  
	   SetSurfaceNodeValueVector(st, m_sfc, nod_vector, nod_val_vector_area);
	}
	else
		std::cout << "GEO_TYPE not supported in CSourceTermGroup::DistributeVolumeFlux()" << "\n";

	for(int i=0; i<(int)nod_val_vector_area.size();i++)
		area +=nod_val_vector_area[i];
	if(area>0)
	{
		for(int i=0; i<(int)nod_val_vector.size();i++)
			nod_val_vector[i] /= area;
	}
	else
	{
		std::cout << "Using the geometric object " << st->geo_name << " does not find any nearby nodes. No ST applied there!" << "\n";
	}

}

/**************************************************************************
 MSHLib-Method:
 Task:
 Programing:
 11/2007 JOD
 last modification:
 **************************************************************************/
void CSourceTerm::SetSurfaceNodeVectorConditional(std::vector<long>&sfc_nod_vector,
std::vector<long>&sfc_nod_vector_cond)
{
   long number_of_nodes;
   number_of_nodes = (long) sfc_nod_vector.size();

   sfc_nod_vector_cond.resize(number_of_nodes);
   SetNOD2MSHNOD(sfc_nod_vector, sfc_nod_vector_cond);
}


/**************************************************************************
 MSHLib-Method:
 Task:
 Programing:
 11/2007 JOD
 last modification:
 **************************************************************************/
void CSourceTerm::SetNodeValues(const std::vector<long>& nodes, const std::vector<long>& nodes_cond,
		const std::vector<double>& node_values, const std::vector<double>& nodes_cond_length, const int& ShiftInNodeVector)
{
   CNodeValue *m_nod_val = NULL;
   size_t number_of_nodes (nodes.size());

   // Added by CB;  removed by JOD 2015-11-19 
   /*double geometry_area = 1.0;
   if (this->getProcessDistributionType() == FiniteElement::CONSTANT_NEUMANN)
   {
     //long msh_ele = fem_msh_vector[0]->nod_vector[nodes[0] + ShiftInNodeVector]->getConnectedElementIDs()[0];
     //if (mmp_vector[fem_msh_vector[0]->ele_vector[msh_ele]->GetPatchIndex()]->GetGeoDimension() < 3)
     //  geometry_area = mmp_vector[fem_msh_vector[0]->ele_vector[msh_ele]->GetPatchIndex()]->geo_area;
     if (mmp_vector[0]->GetGeoDimension() < 3)
       geometry_area = mmp_vector[0]->geo_area;

   }*/


   for (size_t i = 0; i < number_of_nodes; i++)
   {
      m_nod_val = new CNodeValue();
      m_nod_val->msh_node_number = nodes[i] + ShiftInNodeVector;
      m_nod_val->geo_node_number = nodes[i];
      m_nod_val->setProcessDistributionType (getProcessDistributionType());
      m_nod_val->node_value = node_values[i];
      m_nod_val->length = (geo_node_value == 0)? 0 : node_values[i] / geo_node_value;

      m_nod_val->scaling_node_group = scaling_node_group;
      // Added by CB;  removed by JOD 2015-11-19 
      //if (this->getProcessDistributionType() == FiniteElement::CONSTANT_NEUMANN)
      //  m_nod_val->node_value *= geometry_area;
      //m_nod_val->node_value *= geo_node_value;  // JOD 2020-05-29 ?????

      m_nod_val->CurveIndex = CurveIndex;
// CB JOD MERGE //
	
	  if (connected_geometry && // JOD 2/2015
   			connected_geometry_couplingType != 2) // not borehole mit given primary variable
	  {
		  m_nod_val->msh_node_number_conditional = nodes_cond[i];  // JOD 2021-12-10
		  m_nod_val->msh_vector_conditional = nodes_cond;
		  m_nod_val->msh_vector_conditional_length = nodes_cond_length;
	  }

    /**/
	  if (_coupled && nodes_cond.size() > 0)                               // JOD 4.7.10
      {
         m_nod_val->msh_node_number_conditional = nodes_cond[i];
                                                  // JOD 4.10.01
         if ((getProcessType() == FiniteElement::OVERLAND_FLOW || getProcessType() == FiniteElement::GROUNDWATER_FLOW)
        		 && node_averaging)
         {
            double weights = 0;
            for (size_t j = 0; j < number_of_nodes; j++)
            {
               m_nod_val->msh_node_numbers_averaging.push_back(nodes[j]);
               m_nod_val->msh_node_weights_averaging.push_back(
                  node_values[j]);
               weights += node_values[j];
            }
            for (size_t j = 0; j < number_of_nodes; j++)
               m_nod_val->msh_node_weights_averaging[j] /= weights;
         }
      }
      //WW        group_vector.push_back(m_node_value);
      //WW        st_group_vector.push_back(m_st); //OK
      //		if (getProcessDistributionType() == RIVER) {
      //			m_nod_val->node_value = node_value_vectorArea[i];
      //			m_nod_val->node_parameterA = node_value_vectorA[i];
      //			m_nod_val->node_parameterB = node_value_vectorB[i];
      //			m_nod_val->node_parameterC = node_value_vectorC[i];
      //			m_nod_val->node_parameterD = node_value_vectorD[i];
      //			m_nod_val->node_parameterE = node_value_vectorE[i];
      //		}
      //		if (dis_type == 6 || dis_type == 8 || dis_type == 9) // critical depth, normal depth, analytical
      if (getProcessDistributionType() == FiniteElement::CRITICALDEPTH
         || getProcessDistributionType() == FiniteElement::NORMALDEPTH
         || getProcessDistributionType() == FiniteElement::ANALYTICAL)
      {
         m_nod_val->node_value = node_value_vectorArea[i];
                                                  //CMCD bugfix on 4.9.06
         m_nod_val->node_area = node_value_vectorArea[i];
      }

      m_nod_val->setSTVectorIndex(i);
      m_nod_val->setSTVectorGroup(this->getSTVectorGroup());
      if(st_vector[m_nod_val->getSTVectorGroup()]->isConstrainedST())
         st_vector[m_nod_val->getSTVectorGroup()]->setConstrainedSTNodesIndex(m_nod_val->geo_node_number, m_nod_val->getSTVectorIndex());
      _pcs->st_node_value.push_back(m_nod_val);   //WW
      _pcs->st_node.push_back(this);              //WW

      if (_coupled && getProcessType() != FiniteElement::HEAT_TRANSPORT)  // !!!!!!!!!!
      {
    	  CRFProcess* m_pcs_OF = PCSGet("OVERLAND_FLOW");
    	  if(m_pcs_OF)
    	  {
    		  CNodeValue *m_nod_val2 = new CNodeValue();
    		  m_nod_val2->msh_node_number = nodes_cond[i];
    		  m_nod_val2->geo_node_number = nodes_cond[i];
    		  m_nod_val2->setProcessDistributionType (getProcessDistributionType());
    		  m_nod_val2->node_value = node_values[i];
    		  m_nod_val2->length = (geo_node_value == 0)? 0 : node_values[i] / geo_node_value;

    		  m_nod_val2->msh_node_number_conditional = nodes[i];

    		  m_pcs_OF->st_node_value.push_back(m_nod_val2);
    		  m_pcs_OF->st_node.push_back(this);
    	  }
      }
   }  // end nodes
}


// 09/2010 TF
//void CSourceTerm::SetNodeValues(const std::vector<size_t>& nodes, const std::vector<size_t>& nodes_cond,
//		const std::vector<double>& node_values, int ShiftInNodeVector)
//{
//   size_t number_of_nodes (nodes.size());
//
//   for (size_t i = 0; i < number_of_nodes; i++)
//   {
//      CNodeValue *m_nod_val = new CNodeValue();
//      m_nod_val->msh_node_number = nodes[i] + ShiftInNodeVector;
//      m_nod_val->geo_node_number = nodes[i];
//      m_nod_val->setProcessDistributionType (getProcessDistributionType());
//      m_nod_val->node_value = node_values[i];
//      m_nod_val->CurveIndex = CurveIndex;
//
//      if (_coupled)                               // JOD 4.7.10
//      {
//         m_nod_val->msh_node_number_conditional = nodes_cond[i];
//         if ((getProcessType() == OVERLAND_FLOW
//            || getProcessType() == GROUNDWATER_FLOW)
//            && node_averaging)                    // JOD 4.10.01
//         {
//            double weights = 0;
//            for (size_t j = 0; j < number_of_nodes; j++)
//            {
//               m_nod_val->msh_node_numbers_averaging.push_back(nodes[j]);
//               m_nod_val->msh_node_weights_averaging.push_back(
//                  node_values[j]);
//               weights += node_values[j];
//            }
//            for (size_t j = 0; j < number_of_nodes; j++)
//               m_nod_val->msh_node_weights_averaging[j] /= weights;
//         }
//      }
//
//      //		if (getProcessDistributionType() == RIVER) {
//      //			m_nod_val->node_value = node_value_vectorArea[i];
//      //			m_nod_val->node_parameterA = node_value_vectorA[i];
//      //			m_nod_val->node_parameterB = node_value_vectorB[i];
//      //			m_nod_val->node_parameterC = node_value_vectorC[i];
//      //			m_nod_val->node_parameterD = node_value_vectorD[i];
//      //			m_nod_val->node_parameterE = node_value_vectorE[i];
//      //		}
//      if (getProcessDistributionType() == FiniteElement::CRITICALDEPTH
//         || getProcessDistributionType() == FiniteElement::NORMALDEPTH
//         || getProcessDistributionType() == FiniteElement::ANALYTICAL)
//      {
//         m_nod_val->node_value = node_value_vectorArea[i];
//                                                  //CMCD bugfix on 4.9.06
//         m_nod_val->node_area = node_value_vectorArea[i];
//      }
//      _pcs->st_node_value.push_back(m_nod_val);   //WW
//      _pcs->st_node.push_back(this);              //WW
//   }                                              // end nodes
//}


/**************************************************************************
 MSHLib-Method:
 Task:
 Programing:
 11/2005 MB
 last modification:
 **************************************************************************/
//void CSourceTerm::SetNOD2MSHNOD(vector<long>&nodes,
//		vector<long>&conditional_nodes)
//{
//	CGLPoint* m_pnt = NULL;
//	long number;
//	CFEMesh* m_msh_cond = NULL;
//	CFEMesh* m_msh_this = NULL;
//
//	m_msh_cond = FEMGet(pcs_type_name_cond);
//	m_msh_this = FEMGet(convertProcessTypeToString(m_st->getProcessType()));
//	m_pnt = new CGLPoint;
//
//	for (long i = 0; i < (long) nodes.size(); i++) {
//		m_pnt->x = m_msh_this->nod_vector[nodes[i]]->X();
//		m_pnt->y = m_msh_this->nod_vector[nodes[i]]->Y();
//		m_pnt->z = m_msh_this->nod_vector[nodes[i]]->Z();
//
//		number = m_msh_cond->GetNODOnPNT(m_pnt);
//		conditional_nodes[i] = number;
//
//	}
//
//	delete m_pnt;
//
//}

void CSourceTerm::SetNOD2MSHNOD(std::vector<long>& nodes,
		std::vector<long>& conditional_nodes)
{
	CFEMesh* m_msh_cond(FEMGet(pcs_type_name_cond));
	CFEMesh* m_msh_this(FEMGet(convertProcessTypeToString(getProcessType())));

	for (size_t i = 0; i < nodes.size(); i++) {

		const double* coordinate;
		//coordinate = m_msh_this->nod_vector[nodes[i]]->getData();
		//const
		GEOLIB::Point pnt;//(coordinate);

	       if(pcs_type_name_cond == "OVERLAND_FLOW" || pcs_type_name_cond2 == "OVERLAND_FLOW")
			{
					coordinate = m_msh_this->nod_vector[nodes[i]]->getData();
					double coordinate2[3];
					coordinate2[0] = -coordinate[2];
					coordinate2[1] = 0.;
					coordinate2[2] = 10.;//coordinate[2];

					pnt = GEOLIB::Point(coordinate2);
			}
			else
			{
					coordinate = m_msh_this->nod_vector[nodes[i]]->getData();
					//const GEOLIB::Point
					pnt = GEOLIB::Point(coordinate);
			}


		//const GEOLIB::Point pnt(m_msh_this->nod_vector[nodes[i]]->getData());
//      pnt[0] = m_msh_this->nod_vector[nodes[i]]->X();
//      pnt[1] = m_msh_this->nod_vector[nodes[i]]->Y();
//      pnt[2] = m_msh_this->nod_vector[nodes[i]]->Z();

		conditional_nodes[i] = m_msh_cond->GetNODOnPNT(&pnt);
	}
}


void CSourceTerm::SetNOD2MSHNOD(const std::vector<size_t>& nodes,
		std::vector<size_t>& conditional_nodes) const
{
	CFEMesh* m_msh_cond(FEMGet(pcs_type_name_cond));
	CFEMesh* m_msh_this(FEMGet(convertProcessTypeToString(getProcessType())));

	for (size_t i = 0; i < nodes.size(); i++) {
		const GEOLIB::Point pnt(m_msh_this->nod_vector[nodes[i]]->getData());
//		pnt[0] = m_msh_this->nod_vector[nodes[i]]->X();
//		pnt[1] = m_msh_this->nod_vector[nodes[i]]->Y();
//		pnt[2] = m_msh_this->nod_vector[nodes[i]]->Z();

		conditional_nodes[i] = m_msh_cond->GetNODOnPNT(&pnt);
	}
}


/**************************************************************************
 GeoSys source term function:
 02/2009 WW Implementation
 **************************************************************************/
void CSourceTerm::DirectAssign(long ShiftInNodeVector)
{
   CRFProcess* m_pcs = PCSGet(convertProcessTypeToString(getProcessType()));

   if (getProcessDistributionType()==FiniteElement::CLIMATE) //NB for this type of ST, we assign a ST to each node on the Mesh surface (land surface)
   {
		std::vector<double> node_area_vec;
		MshEditor::getNodeAreas(m_pcs->m_msh, node_area_vec);
		const std::vector<GEOLIB::PointWithID*> &points ( MshEditor::getSurfaceNodes(*(m_pcs->m_msh)) );

		size_t nPoints (points.size());
		std::cout << points.size() << " nodes found on mesh surface. " << "\n";
			
		for (size_t i=0; i<nPoints; i++)
		{
			size_t node_id (points[i]->getID());
			CNodeValue* m_nod_val (new CNodeValue());
			m_nod_val->msh_node_number = node_id + ShiftInNodeVector;
			m_nod_val->geo_node_number = node_id;
			m_nod_val->setProcessDistributionType (getProcessDistributionType());
	   	    m_nod_val->node_value = std::numeric_limits<double>::min();  // values will be assigned in IncorporateSoureTerms (rf_pcs.cpp)
			m_nod_val->CurveIndex = CurveIndex;
			m_pcs->st_node_value.push_back(m_nod_val);
			m_pcs->st_node.push_back(this);
			m_pcs->m_msh->nod_vector[node_id]->patch_area = node_area_vec[node_id];
		}
		
		this->_distances = new MathLib::InverseDistanceInterpolation<GEOLIB::PointWithID*, GEOLIB::Station*>(points, this->_weather_stations);
   }
   else //NB this is the old version, where nodes were read from an separate input file
   {
   std::string line_string;
   std::string st_file_name;
   std::stringstream in;
   long n_index;
   double n_val;

   //========================================================================
   // File handling
   std::ifstream d_file(fname.c_str(), std::ios::in);
   //if (!st_file.good()) return;

   if (!d_file.good())
	   {
		  std::cout << "! Error in direct node source terms: Could not find file:!\n" << fname << "\n";
		  abort();
	   }
	   // Rewind the file
	   d_file.clear();
	   d_file.seekg(0L, std::ios::beg);
	   //========================================================================
	   while (!d_file.eof())
	   {
		  line_string = GetLineFromFile1(&d_file);
		  if (line_string.find("#STOP") != std::string::npos)
			 break;

		  in.str(line_string);
		  in >> n_index >> n_val;
		  in.clear();
		  //
		  CNodeValue* m_nod_val (new CNodeValue());
		  m_nod_val->msh_node_number = n_index + ShiftInNodeVector;
		  m_nod_val->geo_node_number = n_index;
		  m_nod_val->setProcessDistributionType (getProcessDistributionType());
		  m_nod_val->node_value = n_val;
		  m_nod_val->CurveIndex = CurveIndex;
		  m_pcs->st_node_value.push_back(m_nod_val);
		  m_pcs->st_node.push_back(this);
		  //
	   }                                              // eof
	}
}


/**************************************************************************
GeoSys source term function:
03/2010 WW Implementation
**************************************************************************/
std::string CSourceTerm::DirectAssign_Precipitation(double current_time)
{
   int i, size;
   double stepA, stepB, tim_val;

   long l, nbc_node, n_index, osize = 0;
   double n_val;

   std::string fileA, fileB;

   CRFProcess* m_pcs = NULL;
   CNodeValue *m_nod_val = NULL;
                                                  //PCSGet(pcs_type_name);
   m_pcs = PCSGet(convertProcessTypeToString (getProcessType()));

   if(start_pos_in_st<0)
      osize = (long)m_pcs->st_node.size();

   size = (int)precip_times.size();
   stepA = 0.;
   stepB = 0.;
   if(current_time < m_pcs->GetTimeStepping()->time_start||fabs(current_time - m_pcs->GetTimeStepping()->time_start)<DBL_MIN)
   {
      fileA = precip_files[0];
      stepB = -1.;
   }
   else if(current_time > precip_times[size-1]||fabs(current_time - precip_times[size-1])<DBL_MIN)
   {
      fileA = precip_files[size-1];
      stepB = -1.;
   }
   else
   {
      double step_b = DBL_MAX;
      double step_f = DBL_MAX;

      for(i=0; i<size; i++)
      {
         tim_val = precip_times[i];
         if(current_time>tim_val)
         {
            if((current_time-tim_val)<step_b)
            {
               step_b = current_time-tim_val;
               stepA = tim_val;
               fileA = precip_files[i];
            }
         }
         else
         {

            if((tim_val-current_time)<step_f)
            {
               step_f = tim_val-current_time;
               stepB = tim_val;
               fileB = precip_files[i];
            }
         }
      }
      if(fabs(stepA-current_time)<DBL_MIN)
         stepB = -1.;
      if(fabs(current_time-stepB)<DBL_MIN)
      {
         fileA = fileB;
         stepB = -1.;
      }
      if(fabs(stepA-stepB)<DBL_MIN)
         stepB = -1.;
   }

   fileA = FilePath + fileA;
   std::ifstream bin_sA(fileA.c_str(), std::ios::binary);
   std::ifstream bin_sB;

   if(!bin_sA.good())
   {
      std::cout<<"Could not find file "<< fileA<<'\n';
      exit(0);
   }
   bin_sA.setf(std::ios::scientific,std::ios::floatfield);
   bin_sA.precision(14);

   /*
   if(stepB>0.)
   {
      fileB = FilePath + fileB;
      bin_sB.open(fileB.c_str(), ios::binary);
      if(!bin_sB.good())
      {
          cout<<"Could not find file "<< fileB<<'\n';
          exit(0);
      }
      bin_sB.setf(ios::scientific,ios::floatfield);
   bin_sB.precision(14);
   bin_sB.read((char*)(&nbc_node), sizeof(nbc_node));
   }
   */

   bin_sA.read((char*)(&nbc_node), sizeof(nbc_node));

   double valA;                                   //, valB = 0.;
   for(l=0; l<nbc_node; l++)
   {

      bin_sA.read((char*)(&n_index), sizeof(n_index));
      bin_sA.read((char*)(&valA), sizeof(valA));
      /*
      if( stepB>0.)
      {
         bin_sB.read((char*)(&n_index), sizeof(n_index));
         bin_sB.read((char*)(&valB), sizeof(valB));
         n_val = valA + (current_time - stepA)*(valB-valA)/(stepB-stepA);
      }
      else
      */
      n_val = valA;

      //
      if(start_pos_in_st<0)
      {
         m_nod_val = new CNodeValue();
         m_pcs->st_node_value.push_back(m_nod_val);
         m_pcs->st_node.push_back(this);
      }
      else
         m_nod_val = m_pcs->st_node_value[l+start_pos_in_st];

      m_nod_val->msh_node_number = n_index ;
      m_nod_val->geo_node_number = n_index;
                                                  //node_distype = dis_type;
      m_nod_val->setProcessDistributionType (getProcessDistributionType());
      m_nod_val->node_value = n_val;
      m_nod_val->CurveIndex = CurveIndex;
      //
   }                                              //

   if(start_pos_in_st<0)
      start_pos_in_st = osize;

   bin_sA.close();
   //bin_sB.close();

   return fileA;
}


/**************************************************************************
 FEMLib-Method:
 Task: Analytical diffusion in matrix. Replaces matrix. See paper to be issued.
 Programing:
 11/2005 CMCD Implementation
 04/2006 Moved from CSourceTermGroup and changed the arguments
 last modification:
 04/2006 CMCD Updated
 **************************************************************************/
                                                  // , CSourceTerm *m_st)
double CSourceTerm::GetAnalyticalSolution(long location)
{
   int idx, n;
   int size, process_no;
   long i;
   long step, no_terms_included;
   double value, source, gradient, ref_value = 0.0;
   double timevalue;
   double fac = 1.0;
   double temp_time, temp_value;
   double pi = 3.1415926;
   double D = this->analytical_diffusion;
   double ne = this->analytical_porosity;
   double tort = this->analytical_tortousity;
   double Kd = this->analytical_linear_sorption_Kd;
   double rho = this->analytical_matrix_density;
   double Dtrans = (D * ne) / ((ne + Kd * rho) * tort);
   double Dsteady = D * ne / tort;
   double t0, tn, tnn, val1, val2, node_area;
   double tvol, vol, flux_area, tflux_area;
   double mass_solute_present, mass_to_remove;
   //WW  bool out = false;
   //WW  int dimension = this->analytical_material_group;
   std::string process;
   CRFProcess* m_pcs = NULL;
   m_pcs = PCSGet(convertProcessTypeToString(this->getProcessType()), convertPrimaryVariableToString(this->getProcessPrimaryVariable()));
   CFEMesh* m_msh = m_pcs->m_msh;                 //WW
   MeshLib::CElem* Ele = NULL;
   long node_number = location;                   //WW m_pcs->st_node_value[location]->msh_node_number;
   CNode* Node = m_msh->nod_vector[node_number];
   double area = m_pcs->st_node_value[location]->node_area;
   std::vector<double> time_history;
   std::vector<double> value_history;
   //Initialise
   time_history.clear();
   value_history.clear();
   t0 = tn = tnn = source = gradient = val1 = val2 = node_area = flux_area
      = 0.0;
   idx = m_pcs->GetNodeValueIndex(convertPrimaryVariableToString(this->getProcessPrimaryVariable()));
   value = m_pcs->GetNodeValue(node_number, idx);
   if (value < MKleinsteZahl)
      value = 0.0;
   timevalue = aktuelle_zeit;
   step = aktueller_zeitschritt;
   if (step < 10)
      step = 10;
   size = (int) analytical_processes.size();
   process_no = 0;
   //Domain or Polyline
   for (i = 0; i < size; i++)
   {
      if (analytical_processes[i] == convertPrimaryVariableToString(this->getProcessPrimaryVariable()))
      {
         if (this->getGeoType () == GEOLIB::POLYLINE)
         {
            if (this->getGeoName().compare(
               analytical_processes_polylines[i]) == 0)
               process_no = i;
         }
         //			if (this->geo_type_name.compare("DOMAIN") == 0)
         if (this->getGeoType () == GEOLIB::GEODOMAIN)
            process_no = i;
      }
   }
   //Identify process
   process_no *= 2;                               //first column time, second column value, hence two columns per process;

   //If time step require new calculation of source term then start
   if ((aktueller_zeitschritt - 1) % this->resolution == 0) {
      //Save data in a vector attached to the nodes
      this->SetNodePastValue(node_number, process_no, 0, timevalue);
      this->SetNodePastValue(node_number, process_no + 1, 0, value);

      //Recall historical data
      ref_value = this->GetNodePastValueReference(node_number, (process_no
         / 2));
      if ((size_t)step > this->number_of_terms)
         no_terms_included = this->number_of_terms;
      else
         no_terms_included = step;
      for (i = 0; i < no_terms_included; i++)
      {
         temp_time = this->GetNodePastValue(node_number, process_no, i);
         temp_value
            = (this->GetNodePastValue(node_number, process_no + 1, i)
            - ref_value);
         time_history.push_back(temp_time);
         value_history.push_back(temp_value);
      }

      //Calculate individual terms and sum for gradient
      for (i = no_terms_included - 1; i > 0; i--)
      {
         t0 = time_history[0];
         if (i == no_terms_included - 1)
            tn = (t0 - time_history[i]) + (time_history[i - 1]
               - time_history[i]);
         else
            tn = t0 - time_history[i + 1];
         tnn = t0 - time_history[i];
         val1 = 1 / (sqrt(pi * Dtrans * tn));
         val2 = 1 / (sqrt(pi * Dtrans * tnn));
         gradient += ((val2 - val1) * value_history[i]);
      }
      tn = t0 - time_history[1];
      tnn = 0;
      val1 = 1 / (sqrt(pi * Dtrans * tn));
      gradient -= (val1 * value_history[0]);

      //Area calculations
      mass_solute_present = 1.e99;                //Initially set value very high

      //Area for lines, triangles and quads in domain.
      //  if (area < 0) {//set in set source terms function, domain area = -1 to start with
      if (area < DBL_MIN) { // HS 04.2008
			tvol = 0.0;
			tflux_area = 0.0;
			for (i = 0; i < (int) Node->getConnectedElementIDs().size(); i++) {
				Ele = m_msh->ele_vector[Node->getConnectedElementIDs()[i]];
				vol = Ele->GetVolume(); //Assuming 1m thickness
				flux_area = Ele->GetFluxArea(); //Real thickness for a fracture
				n = Ele->GetVertexNumber();
				tvol += (vol / n);
				tflux_area += (flux_area / n);
			}
			node_area = tvol * 2.; //* 2 because the diffusion is in two direction perpendicular to the fracture
			mass_solute_present = tflux_area * tvol * value;
		}
		//Area for polylines
		else
		 node_area = area;

      //factor for conversion to energy for temperature if necessary
      fac = this->factor;
      source = gradient * Dsteady * fac * node_area;
      mass_to_remove = fabs(source) * dt;
      if (mass_to_remove > mass_solute_present)
      {
         source *= (mass_solute_present / mass_to_remove);
      }
      this->SetNodeLastValue(node_number, (process_no / 2), source);

   }                                              // If new source term calculation not required
   else
      source = this->GetNodeLastValue(node_number, (process_no / 2));
   return source;
}


/**************************************************************************
 FEMLib-Method:
 Task: master write function
 Programing:
 11/2005 CMCD Implementation, functions to access and write the time history data
 of nodes
 last modification:
 **************************************************************************/
void CSourceTerm::SetNodePastValue(long n, int idx, int pos, double value)
{
   bool endstepreached = false;
   pos = pos;                                     //WW

   //Check whether this is the first call
   CRFProcess* m_pcs = NULL;
   m_pcs = PCSGet(convertProcessTypeToString(getProcessType()), convertPrimaryVariableToString (getProcessPrimaryVariable()));
   if (!m_pcs)                                    //OK
   {
      std::cout << "Warning in SetNodePastValue - no PCS data" << '\n';
      return;
   }

   size_t size1 (m_pcs->nod_val_vector.size());

   //Create memory for the values
   if (node_history_vector.empty())
   {
      //WW     int number_of_terms = max_no_terms;
      for (size_t k = 0; k < size1; k++)
      {
         NODE_HISTORY *nh = new NODE_HISTORY;
         CreateHistoryNodeMemory(nh);
         node_history_vector.push_back(nh);
      }
      for (size_t k = 0; k < size1; k++)
      {
         for (size_t j = 0; j < _no_an_sol; j++)
            node_history_vector[k]->value_reference.push_back(-1.0);
      }
   }

   //Store the first set of values as reference values
   int flipflop = idx % 2;
   if (aktueller_zeitschritt == 1)
   {
      //if (size2 == idx)
      if (flipflop == 1)
         node_history_vector[n]->value_reference[(idx - 1) / 2] = value;
   }
   //  size2 = (int) node_history_vector[n]->value_reference.size();
   size_t steps = 10;
   if (aktueller_zeitschritt > steps)
      steps = aktueller_zeitschritt;
   if (_max_no_terms >= steps)
      steps = aktueller_zeitschritt;
   else
   {
      steps = _max_no_terms;
      endstepreached = true;
   }
   //Enter the value and push the other values back
   if (!endstepreached)
   {
      for (size_t k = steps - 1; k > 0; k--)
         node_history_vector[n]->value_store[idx][k]
            = node_history_vector[n]->value_store[idx][k - 1];
      node_history_vector[n]->value_store[idx][0] = value;
   }
   double cutvalue = 0.0;
   double nextvalue = 0.0;
   long no_steps_past_cutof = 0;
   if (endstepreached)
   {
      no_steps_past_cutof = aktueller_zeitschritt - _max_no_terms;
      cutvalue = node_history_vector[n]->value_store[idx][steps - 1];
      nextvalue = node_history_vector[n]->value_store[idx][steps - 2];
      node_history_vector[n]->value_store[idx][steps - 1] = (cutvalue
         * (double) (no_steps_past_cutof) + nextvalue)
         / ((double) no_steps_past_cutof + 1);
      for (size_t k = steps - 2; k > 0; k--)
         node_history_vector[n]->value_store[idx][k]
            = node_history_vector[n]->value_store[idx][k - 1];
      node_history_vector[n]->value_store[idx][0] = value;
   }
}


void CSourceTerm::SetNodeLastValue(long n, int idx, double value)
{
   size_t size (node_history_vector[n]->last_source_value.size());
   if (size == 0)
   {
      for (size_t i = 0; i < _no_an_sol; i++)
         node_history_vector[n]->last_source_value.push_back(0);
   }
   node_history_vector[n]->last_source_value[idx] = value;
}


double CSourceTerm::GetNodeLastValue(long n, int idx)
{
   double value = 0.0;
   value = node_history_vector[n]->last_source_value[idx];
   return value;
}


double CSourceTerm::GetNodePastValue(long n, int idx, int pos)
{
   double value;
   value = node_history_vector[n]->value_store[idx][pos];
   return value;
}


double CSourceTerm::GetNodePastValueReference(long n, int idx)
{
   double value;
   value = node_history_vector[n]->value_reference[idx];
   return value;
}


void CSourceTerm::CreateHistoryNodeMemory(NODE_HISTORY* nh)
{
   size_t s_col = _no_an_sol * 2;
   size_t s_row = number_of_terms;

   nh->value_store = new double*[s_col];
   for (size_t i = 0; i < s_col; i++)
      nh->value_store[i] = new double[s_row];

   for (size_t k = 0; k < s_col; k++)
   {
      for (size_t j = 0; j < s_row; j++)
         nh->value_store[k][j] = 0.0;
   }
}


void CSourceTerm::DeleteHistoryNodeMemory()
{
   size_t size (node_history_vector.size());
   size_t s_row = _no_an_sol * 2;

   if (size > 0)
   {
      for (size_t j = 0; j < size; j++)
      {
         for (size_t i = 0; i < s_row; i++)
            delete node_history_vector[j]->value_store[i];
         delete node_history_vector[j]->value_store;
      }
      node_history_vector.clear();
   }
}


///**************************************************************************
//FEMLib-Method:
//07/2007 OK Implementation
//modifications:
//05/2010 TF restructured a little bit
//**************************************************************************/
//CSourceTerm* STGet(const std::string& pcs_name, const std::string& geo_type_name, const std::string& geo_name)
//{
//  for(size_t i=0;i<st_vector.size();i++) {
//    if((st_vector[i]->pcs_type_name.compare(pcs_name)==0) &&
//       (st_vector[i]->geo_type_name.compare(geo_type_name)==0) &&
//       (st_vector[i]->getGeoName().compare(geo_name)==0))
//      return st_vector[i];
//
////    if((st_vector[i]->pcs_type_name.compare(pcs_name)==0) &&
////    		(st_vector[i]->geo_type_name.compare(geo_type_name)==0) &&
////    		geo_type_name.compare ("POINT") == 0 &&
////    		(st_vector[i]->getGeoObjIdx() == compare(geo_name)==0))
////    	return st_vector[i];
//  }
//  return NULL;
//}

/**************************************************************************
 FEMLib-Method: 4.7.10 shift and average field variables
 10/2008 JOD Implementation
 12/2012 JOD extension to two-phase flow  5.3.07
 gives heads for liquids and pressures for gases
 overland / groundwater nodes are shifted to a soil column if required
 **************************************************************************/
void GetCouplingFieldVariables( CRFProcess* m_pcs_this,  CRFProcess* m_pcs_cond, double* h_this, double* h_cond,
double* h_shifted, double* h_averaged, double* z_this, double* z_cond,
CSourceTerm* m_st, CNodeValue* cnodev, long msh_node_number, long msh_node_number_cond, double gamma)
{

   int nidx, nidx_cond;

   if (m_st->getProcessType() == FiniteElement::HEAT_TRANSPORT
	   || m_st->getProcessType() == FiniteElement::MASS_TRANSPORT)
   {
	   *z_this = 0;
	   *z_cond = 0;
   }
   else
   {
	   *z_this = m_pcs_this->m_msh->nod_vector[msh_node_number]->getData()[2];
	   *z_cond  = m_pcs_cond->m_msh->nod_vector[msh_node_number_cond]->getData()[2];
   }
   nidx = m_pcs_this->GetNodeValueIndex (convertPrimaryVariableToString (m_st->getProcessPrimaryVariable())) + 1;
   nidx_cond = m_pcs_cond->GetNodeValueIndex(m_st->pcs_pv_name_cond) + 1;

   *h_this = m_pcs_this->GetNodeValue(msh_node_number, nidx);

   if(m_st->pcs_pv_name_cond == "GIVEN_HEIGHT" || m_st->pcs_pv_name_cond == "PRESSURE2")
   {
	 *h_cond = *h_shifted = m_st->coup_given_value; // coupled to fixed pressure head / or atmospheric pressure
	 *z_this = *z_cond = 0;
	 return;
   }
   else if(m_st->pcs_pv_name_cond == "GIVEN_PRESSURE")
  {
	 *h_cond = *h_shifted = m_st->coup_given_value / gamma; // coupled to fixed pressure head / or atmospheric pressure
	 *z_this = *z_cond = 0;
	 return;
   }
   else
     *h_cond = m_pcs_cond->GetNodeValue(msh_node_number_cond, nidx_cond); // coupled to a liquid

   if(m_st->pcs_pv_name_cond == "PRESSURE1") {
	
	   if(m_st->pcs_type_name_cond == "MULTI_PHASE_FLOW")
	   { // capillary pressure as a primary variable in the coupled process 
	     double gasPressure =  m_pcs_cond->GetNodeValue(msh_node_number_cond, 3);
	     *h_cond -= (gasPressure - m_st->coup_given_value);
      }
   }

   if (m_st->getProcessType() == FiniteElement::OVERLAND_FLOW 
	   || m_st->getProcessType() == FiniteElement::GROUNDWATER_FLOW)
   {
      if (m_st->node_averaging)                 
      { // average pressures from overland / groundwater (soil column coupled to multiple nodes)
         *h_shifted = *h_this - *z_this + *z_cond;  // shift overland / groundwater node to soil column
         *h_averaged = 0;
         for (long i = 0; i
            < (long) cnodev->msh_node_numbers_averaging.size(); i++)
            *h_averaged
               += cnodev->msh_node_weights_averaging[i]
               * (m_pcs_this->GetNodeValue(
               cnodev->msh_node_numbers_averaging[i],
               nidx)
               - m_pcs_this->m_msh->nod_vector[cnodev->msh_node_numbers_averaging[i]]->getData()[2]);

         *h_averaged += *z_cond; // convert to head
         *z_this = *z_cond;
      }                                           // end averaging
      else                                        // no averaging
      {   // no column, no shifting
         *h_shifted = *h_this;
         *h_averaged = *h_this;
      }                                           // end no averaging

      if (m_st->pcs_pv_name_cond == "PRESSURE1")
      {       // convert to head
		// if (m_st->pcs_type_name_cond == "MULTI_PHASE_FLOW")
          // *h_cond /= -gamma;
		 //else
		   *h_cond /= gamma;
         
		 *h_cond += *z_cond;
      }
      if (m_st->pcs_type_name_cond == "GROUNDWATER_FLOW")
         h_cond = std::max(h_cond, z_this);       //groundwater level might not reach surface  
   }          // end overland flow
   else      
   {     // subsurface flow
      if (m_st->pcs_pv_name_cond == "PRESSURE1")  // JOD 4.10.01
      { // convert to head
          //if (m_st->pcs_type_name_cond == "MULTI_PHASE_FLOW")
          // *h_cond /= -gamma;
		 //else
		   *h_cond /= gamma;

         *h_cond += *z_cond;
      }
      if (m_st->node_averaging)                   
      {             // shift overland / groundwater node to the soil column
         *h_shifted = *h_cond - *z_cond;
         *h_shifted += *z_this;
         *z_cond = *z_this;
      }                                           // end averaging
      else
         *h_shifted = *h_cond;  // no column, no shifting

      if (m_st->getProcessPrimaryVariable() == FiniteElement::PRESSURE)
      {   // convert pressure into head
         *h_this /=  gamma;   
         *h_this += *z_this;
      }
   }                                              // end subsurface flow

}

/**************************************************************************
 FEMLib-Method: 4.7.10
 10/2008 JOD Implementation
 **************************************************************************/
double CalcCouplingValue(double factor, double h_this, double h_cond,
double z_cond, CSourceTerm* m_st)
{

   if (m_st->getProcessType() == FiniteElement::OVERLAND_FLOW)
   {
      if (m_st->no_surface_water_pressure)        
         return factor * (h_cond - z_cond); // neglect hydrostatic surface liquid pressure
      else
         return factor * (h_cond - h_this);
   }                                              // Richards' & groundwater flow
   else
   {
      if (m_st->getProcessType() == FiniteElement::GROUNDWATER_FLOW)
      {
         if (h_this < z_cond && m_st->pcs_type_name_cond == "OVERLAND_FLOW")
            return factor * (h_cond - z_cond); // decoupled
         else
            return factor * h_cond;
      }
      else if (m_st->getProcessType() == FiniteElement::LIQUID_FLOW)
    	  return factor * h_cond;
      else
    	  return factor * (h_cond - z_cond);  // Richards', two-phase flow (liquid & gas)
   }
}


/**************************************************************************
MSHLib-Method:
Task:
Programing:
02/2015 JOD
last modification: 10 / 2018 use GeoObj
**************************************************************************/
void CSourceTerm::SetSurfaceNodeVectorConnected(std::vector<long>&sfc_nod_vector,
	std::vector<long>&sfc_nod_vector_cond)
{
	std::vector<std::size_t> sfc_node_cond_ids;
	CFEMesh* m_msh(FEMGet(convertProcessTypeToString(getProcessType())));

	GEOLIB::Surface const& sfc_connected(
			*(dynamic_cast<GEOLIB::Surface const*>(geoInfo_connected->getGeoObj()))
	      );

	m_msh->GetNODOnSFC(&sfc_connected, sfc_node_cond_ids, true);
	sfc_nod_vector_cond.insert(sfc_nod_vector_cond.begin(),
	    		  sfc_node_cond_ids.begin(), sfc_node_cond_ids.end());

	// !!!!! now something to guarantee match of nodes

}

/**************************************************************************
MSHLib-Method:
Task:
Programing:  same as CSourceTerm::SetSurfaceNodeVectorConnected()
10/2018 JOD
last modification:
**************************************************************************/
void CSourceTerm::SetPolylineNodeVectorConnected(std::vector<long>&ply_nod_vector,
	std::vector<long>&ply_nod_vector_cond)
{
	std::vector<std::size_t> ply_node_cond_ids;
	CFEMesh* m_msh(FEMGet(convertProcessTypeToString(getProcessType())));

	GEOLIB::Polyline const& ply_connected(
				*(dynamic_cast<GEOLIB::Polyline const*>(geoInfo_connected->getGeoObj()))
		      );

	m_msh->GetNODOnPLY(&ply_connected, ply_node_cond_ids, true);
	ply_nod_vector_cond.insert(ply_nod_vector_cond.begin(),
		    		  ply_node_cond_ids.begin(), ply_node_cond_ids.end());

	// !!!!! now something to guarantee match of nodes
}

/**************************************************************************
FEMLib-Method:
Task:	Incorporate a direct connection between two geometries
which is not represented by mesh elements
Programing:
02/2015 SB Implementation
02/2015 JOD4cST Modification to use MKL PARDSIO solver 
**************************************************************************/

void IncorporateConnectedGeometries(double &value, CNodeValue* cnodev, CSourceTerm* m_st, CRFProcess* m_pcs)
{
	// get data from input *.st

	double alpha;

	if(m_st->GetBoreholeMode() == 1)  // advective heat transport
	{

      	if(mfp_vector[0]->get_flag_volumetric_heat_capacity()) // 1st mfp property for liquid !!!
      		alpha = mfp_vector[0]->get_volumetric_heat_capacity();
      	else
			alpha = mfp_vector[0]->SpecificHeatCapacity() * mfp_vector[0]->Density();

		if(m_st->isCoupled())
		{
            CRFProcess* m_pcs_flow = PCSGet(m_st->pcs_type_name_cond);

            CRFProcess* m_pcs_flow2 =  (m_st->pcs_type_name_cond2 == "")? PCSGet(m_st->pcs_type_name_cond):  PCSGet(m_st->pcs_type_name_cond2);


			if(m_pcs_flow)
			{
				const double gamma =   (m_st->pcs_type_name_cond2 == "OVERLAND_FLOW")? 9810.: 1.;
double pressure_cond;
if(m_st->pcs_type_name_cond2 == "OVERLAND_FLOW")
{
	pressure_cond = m_pcs_flow2->GetNodeValue(cnodev->msh_node_number_conditional, 1) - 2*  m_pcs_flow2->m_msh->nod_vector[cnodev->msh_node_number_conditional]->getData()[2] +
			m_pcs_flow2->m_msh->nod_vector[cnodev->msh_node_number_conditional]->getData()[0];
	pressure_cond *= 9810;
}
else
pressure_cond = 	m_pcs_flow2->GetNodeValue(cnodev->msh_node_number_conditional, 1);


				//long no = cnodev->msh_node_number;
				//long no2 = cnodev->msh_node_number_conditional;
				//std::cout << no <<std::endl;
				if(m_pcs_flow->ST_factor_kept.find(cnodev->msh_node_number) != m_pcs_flow->ST_factor_kept.end())
				{
					const double fluid_flux = m_pcs_flow->ST_factor_kept[cnodev->msh_node_number] *  // without area to node (value)
						(pressure_cond - m_pcs_flow->GetNodeValue(cnodev->msh_node_number, 1)) ;

					if(m_st->getConnectedGeometryVerbosity() > 1)
					{
						std::cout << "\t\t\tNode " << cnodev->msh_node_number << " with pressure " << 
						m_pcs_flow->GetNodeValue(cnodev->msh_node_number, 1) << " connected to\n\t\t\tnode " << 
				      		cnodev->msh_node_number_conditional << " with pressure " << 
							pressure_cond<< '\n';
	                	std::cout << "\t\t\tFluid volumetric heat capacity: " << alpha << '\n';
	                	std::cout << "\t\t\tFluid flux: " << fluid_flux << '\n';
					}

					alpha *= fluid_flux;
				}
				else
        	        		throw std::runtime_error("Error in IncorporateConnectedGeometry for borehole - Factor for LIQUID_FLOW not kept");
			}
        	       	else
        	        	throw std::runtime_error("Error in IncorporateConnectedGeometry for borehole - No LIQUID_FLOW");


		}
		else  // not coupled to LIQUID_FLOW - source / sink term used with $KEEP_VALUES
		{
               		CRFProcess* m_pcs_liquid = PCSGet("LIQUID_FLOW");

        	        if(m_pcs_liquid->ST_values_kept.find(cnodev->msh_node_number) != m_pcs_liquid->ST_values_kept.end())
        	       	{
        	       		alpha *= m_pcs_liquid->ST_values_kept[cnodev->msh_node_number];

        	             	if(m_st->getConnectedGeometryVerbosity() > 0)
        	               		std::cout << " - Q_fluid: " << m_pcs_liquid->ST_values_kept[cnodev->msh_node_number] << '\n';
        	              
        	              	if(fabs(m_pcs_liquid->ST_values_kept[cnodev->msh_node_number]) < 1.e-10)  // switch off ST - no connection
				{
					alpha = 0.;
				}
        	     	}
        	       	else
        	         	throw std::runtime_error("Error in IncorporateConnectedGeometry for borehole - No LIQUID_FLOW");
		}
	}
	else
		alpha = m_st->connected_geometry_exchange_term; // unit [1/m]


	double alpha_value = alpha * value;               // value is area to node (distance between connected nodes can be added into preprocessing SetST())

	// determine downwind node part I - select now msh_node_number obtained from input *.gli in SetST()   - nodes are fixed for modes 0, 1 - if mode 2, nodes will be rearranged according to velocity later on)

	const long ToNode = (alpha_value > 0) ?  cnodev->msh_node_number : cnodev->msh_node_number_conditional;
	long FromNode = (alpha_value > 0) ? cnodev->msh_node_number_conditional : cnodev->msh_node_number; 

//std::cout << "value: " <<value << ", to " << ToNode << " from " << FromNode << '\n';
	if(m_st->getConnectedGeometryCouplingType() == 2) // borehole with given primary variable
		FromNode = -1;

	// now we have all data
	m_pcs->IncorporateNodeConnectionSourceTerms(FromNode, ToNode, alpha_value, m_st, value);
}

/**************************************************************************
FEMLib-Method:
Task:	Switch source term off if primary variable of transport process
        passes an upper or lower threshold
Programing:
01/2018 JOD implementation  - restricted to LIQUID_FLOW with HEAT_TRANSPORT
							 temperature taken either explicit or implicit
**************************************************************************/

double CSourceTerm::CheckThreshold(const double &value, const CNodeValue* cnodev) const
{
        CRFProcess* m_pcs = PCSGet(threshold.process);  // HEAT_TRANSPORT
        double distance, result, running_value;
        int sign;
    	long node_number;
    	if(threshold_geometry == true)
    		node_number = msh_node_number_threshold;  // only for point now
    	else
    		node_number = cnodev->msh_node_number;

    	running_value = m_pcs->GetNodeValue(node_number, threshold.scheme);  // temperature

    	if(threshold.type == Threshold::lower)
                sign = 1;
        else if(threshold.type == Threshold::upper)
                sign = -1;
        else
        {
                throw runtime_error("SourceTerm::CheckThreshold - Threshold type must be 1 or 2");
        }

        if(m_pcs == NULL)
        {
                throw runtime_error("SourceTerm::CheckThreshold - Process unknown");
        }
        else
        {
                distance = sign * (running_value - threshold.value);  // positive if source-term should be on, else negative

                if(threshold.scheme == Threshold::_explicit)  // explicit
                { // hard switch
                        if(distance > 0)
                                result = value;  // source-term on (completely)
                        else
                                result = 0;  // source-term off
                }
                else if(threshold.scheme == Threshold::_implicit)  // implicit
                { // soft switch by smoothing - generate S-curve between threshold.value and threshold.value + threshold.delta
                        double relativeDistance = distance / threshold.delta;

                        if(relativeDistance < 0)
                        {
                                result = 0;  // source-term off
                        }
                        else if(relativeDistance > 1)
                                result = value;  // source-term on (completely)
                        else
                        {
                        	result = pow(relativeDistance, 2*(1-relativeDistance)) * value;  // source-term partly on
                        }
                }
                else
                {
                        throw runtime_error("SourceTerm::CheckThreshold - Threshold scheme must be 0 (explicit) or 1 (implicit)");
                }
        }

        if (threshold.verbosity > 1)
        {
        	std::cout << "    Source term with ";
        	if(threshold.type == Threshold::lower)
        		std::cout << "lower ";
        	else if(threshold.type == Threshold::upper)
        		std::cout << "upper ";
        	std::cout << "threshold temperature of " << threshold.value << "\n";
        	std::cout << " \tTemperature at reference node " << node_number << " : " << running_value << "\n";
        }
        if (threshold.verbosity > 0)
        {
        	std::cout << "\tSource term at node " << cnodev->msh_node_number << " : " << result << "\n";
        }
        return result;

}


/**************************************************************************
FEMLib-Method:
Task:	Calculates an LIQUID_FLOW flow for a given power value
        Temperature from HEAT_TRANSPORT
        can be extended to MASS_TRANSPORT
Programing:
02/2018 JOD implementation
**************************************************************************/

double CSourceTerm::CalculateFromStorageRate(const double &value, const CNodeValue* cnodev) const
{
		CRFProcess* m_pcs = PCSGet(getProcessType());
        CRFProcess* m_pcs_transport = PCSGet(storageRate.process);
        double densityInlet, specCapacityInlet;
        double densityOutlet, specCapacityOutlet;
        double primValsInlet[3]; // for density and specific capacity
        double primValsOutlet[3];
		double factor, result;

        if(m_pcs == NULL)
        {
                throw runtime_error("SourceTerm::CalculateFromStorageRate - Process unknown");
                return 0;
        }

        primValsInlet[0] = 0.;  // pressure
        primValsInlet[1] = 0.;  // temperature, process must be HEAT_TRANSPORT
        //primValsInlet[2] : concentration in density calculation and saturation in capacity calculation !!!!!

        // averaging over surface, one node in case of point
        for(int i=0; i < storageRate.inlet_msh_node_numbers.size(); i++)
        {
        	primValsInlet[0] += m_pcs->GetNodeValue(storageRate.inlet_msh_node_numbers[i], 1)
        			* storageRate.inlet_msh_node_areas[i];
        	primValsInlet[1] += m_pcs_transport->GetNodeValue(storageRate.inlet_msh_node_numbers[i], 1)
        			* storageRate.inlet_msh_node_areas[i];
        }
        primValsInlet[0] /= storageRate.inlet_totalArea;
        primValsInlet[1] /= storageRate.inlet_totalArea;

        densityInlet = mfp_vector[0]->Density(primValsInlet);
        specCapacityInlet = mfp_vector[0]->SpecificHeatCapacity(primValsInlet);

        primValsOutlet[0] = 0.;  // pressure
        primValsOutlet[1] = 0.;  // temperature, process must be HEAT_TRANSPORT
        //primValsInlet[2] : concentration in density calculation and saturation in capacity calculation !!!!!

        // averaging over surface, one node in case of point
        for(int i=0; i < storageRate.outlet_msh_node_numbers.size(); i++)
        {
        	primValsOutlet[0] += m_pcs->GetNodeValue(storageRate.outlet_msh_node_numbers[i], 1)
        			* storageRate.outlet_msh_node_areas[i];
        	primValsOutlet[1] += m_pcs_transport->GetNodeValue(storageRate.outlet_msh_node_numbers[i], 1)
        			* storageRate.outlet_msh_node_areas[i];
        }
        primValsOutlet[0] /= storageRate.outlet_totalArea;
        primValsOutlet[1] /= storageRate.outlet_totalArea;


        densityOutlet = mfp_vector[0]->Density(primValsOutlet);
        specCapacityOutlet = mfp_vector[0]->SpecificHeatCapacity(primValsOutlet);

        factor = densityInlet * specCapacityInlet * primValsInlet[1] - densityOutlet * specCapacityOutlet * primValsOutlet[1];

        if(fabs(factor) <std::numeric_limits<double>::epsilon())
        	return 0;  // no division by zero - maximum flow rate treated below

        result = storageRate.inputValue / factor;
        if(result > storageRate.absMaximum)
        	result = storageRate.absMaximum;
        else if(result < -storageRate.absMaximum)
        	result = -storageRate.absMaximum;

        if(storageRate.verbosity > 0)
        {
        	std::cout << "\tlux from storage rate at node " << cnodev->msh_node_number << "\n";
        	std::cout << " \t\tSource term value: " << result << " (at time "
        			<< aktuelle_zeit  << ")"<<  "\n";
        }
        if(storageRate.verbosity > 1)
        {
        	std::cout << "\t\t\tWarm well / pipe   - node " << storageRate.inlet_msh_node_numbers[0] << "\n";
        	std::cout << "\t\t\t\tdensity\t: " << densityInlet << "\n";
        	std::cout << "\t\t\t\tspecific capacity\t: " << specCapacityInlet << "\n";
        	std::cout << "\t\t\t\ttemperature\t: " << primValsInlet[1] << "\n";
        	std::cout << "\t\t\tCold well / pipe   - node " << storageRate.outlet_msh_node_numbers[0] << "\n";
        	std::cout << "\t\t\t\tdensity\t: " << densityOutlet << "\n";
        	std::cout << "\t\t\t\tspecific capacity\t: " << specCapacityOutlet << "\n";
        	std::cout << "\t\t\t\ttemperature\t: " << primValsOutlet[1] << "\n";
        	std::cout << "\t\t\tThermal flux (input) : " << factor * result << "\n";
        }

        return value * result;
}


/**************************************************************************
FEMLib-Method:
Task:
Programing:
06/2018 JOD implementation
// passes temperatures and volumetic heat capacities to WDC and returns WDC result

!!! use constant (e.g. volumetric) heat capacity, when using polylines for warm and cold well,
               dependency on pressure and temeprature not considered in NNNC
heat capacity (pressure, temperature) when point heat source /sink term for heat exchanger
only heat capacity and density of fluid is considered (mfp_vector[0])

**************************************************************************/
double CSourceTerm::apply_wellDoubletControl(double value,
				const CNodeValue* cnodev, const double& aktuelle_zeit, CRFProcess* m_pcs)
{
	if(!m_pcs)
		throw std::runtime_error("WellDoubletControl - No PCS");

	const int ndx1 = 1;  // implicit - take new values for capacity calculations

	std::vector<size_t> heatExchanger_aquifer_mesh_nodes, upwind_aquifer_mesh_nodes;
	std::vector<double> heatExchanger_aquifer_mesh_nodes_area_fraction, upwind_aquifer_mesh_nodes_area_fraction;

	int operation_type = ogs_WDC->get_aquifer_mesh_nodes(aktuelle_zeit, wdc_flag_extract_and_reinject,
				heatExchanger_aquifer_mesh_nodes, heatExchanger_aquifer_mesh_nodes_area_fraction,
				upwind_aquifer_mesh_nodes, upwind_aquifer_mesh_nodes_area_fraction);
	// operation_type 1: storing, -1: retrieving
	if(wdc_flag_extract_and_reinject // polyline approach, no connection by elements
		&& m_pcs->getProcessType() == FiniteElement::HEAT_TRANSPORT)
	{
		if(operation_type * value < 0)  // storing with injection at warm well and retrieving with injection at cold well
			return 0;
		if(value < 0)  // -1 is only flag for retrieval: set value positive
			value = -value;
	}

	ogs_WDC->set_heat_exchanger_mesh_nodes(heatExchanger_aquifer_mesh_nodes, heatExchanger_aquifer_mesh_nodes_area_fraction);

	/*std::cout << "heat_exchanger\n";
	for(int i=0; i< heatExchanger_aquifer_mesh_nodes.size(); ++i)
	{
		std::cout << heatExchanger_aquifer_mesh_nodes[i] << " " << heatExchanger_aquifer_mesh_nodes_area_fraction[i] << "\n";
	}

	std::cout << "unpwind\n";
	for(int i=0; i< upwind_aquifer_mesh_nodes.size(); ++i)
	{
		std::cout << upwind_aquifer_mesh_nodes[i] << " " << upwind_aquifer_mesh_nodes_area_fraction[i] << "\n";
	}*/

	const CRFProcess* const m_pcs_liquid = (m_pcs->getProcessType() ==
			FiniteElement::LIQUID_FLOW)? m_pcs : PCSGet("LIQUID_FLOW");
	CRFProcess* m_pcs_heat = (m_pcs->getProcessType() ==
			FiniteElement::HEAT_TRANSPORT)? m_pcs : PCSGet("HEAT_TRANSPORT");

	// heat capacity(pressure, temperature) here considered
	double variables_heatExchanger[3] {
		//ogs_WDC->get_extremum(m_pcs_liquid, ndx1, ogs_WDC->get_doublet_mesh_nodes().heatExchanger),
		//ogs_WDC->get_extremum(m_pcs_heat, ndx1, ogs_WDC->get_doublet_mesh_nodes().heatExchanger)
		m_pcs_liquid->GetWeightedAverageNodeValue(heatExchanger_aquifer_mesh_nodes, // pressure
				heatExchanger_aquifer_mesh_nodes_area_fraction, ndx1),  // !?
		m_pcs_heat->GetWeightedAverageNodeValue(heatExchanger_aquifer_mesh_nodes,   // temperature
				heatExchanger_aquifer_mesh_nodes_area_fraction, ndx1) };
	double variables_upwindAquifer[3] {
		m_pcs_liquid->GetWeightedAverageNodeValue(upwind_aquifer_mesh_nodes,  // pressure
				upwind_aquifer_mesh_nodes_area_fraction, ndx1),  // !?
		m_pcs_heat->GetWeightedAverageNodeValue(upwind_aquifer_mesh_nodes,   // temperature
				upwind_aquifer_mesh_nodes_area_fraction, ndx1) };

	double volumetric_heat_capacity_heatExchanger = (mfp_vector[0]->get_flag_volumetric_heat_capacity()) ?
			mfp_vector[0]->get_volumetric_heat_capacity() :
			mfp_vector[0]->SpecificHeatCapacity(&variables_heatExchanger[1]) * mfp_vector[0]->Density(&variables_heatExchanger[1]);
	double volumetric_heat_capacity_upwindAquifer = (mfp_vector[0]->get_flag_volumetric_heat_capacity()) ?
			mfp_vector[0]->get_volumetric_heat_capacity() :
			mfp_vector[0]->SpecificHeatCapacity(&variables_upwindAquifer[1]) * mfp_vector[0]->Density(&variables_upwindAquifer[1]);

	double result = ogs_WDC->call_WDC(m_pcs,
			{ //ogs_WDC->get_extremum(m_pcs_heat, ndx1, ogs_WDC->get_doublet_mesh_nodes().heatExchanger),  // extremum temperature at heat exchanger
				m_pcs_heat->GetWeightedAverageNodeValue(heatExchanger_aquifer_mesh_nodes,
						heatExchanger_aquifer_mesh_nodes_area_fraction, ndx1),
				m_pcs_heat->GetWeightedAverageNodeValue(upwind_aquifer_mesh_nodes,
					upwind_aquifer_mesh_nodes_area_fraction, ndx1),  // temperature in upwind aquifer
				volumetric_heat_capacity_heatExchanger,
				volumetric_heat_capacity_upwindAquifer },
				heatExchanger_aquifer_mesh_nodes);

	// NNNC
	if((ogs_WDC->get_parameter_list().begin()->indicator == 0 
		|| ogs_WDC->get_parameter_list().begin()->indicator == 1 
		|| ogs_WDC->get_parameter_list().begin()->indicator == 2)  // ST / ST
			&& wdc_flag_extract_and_reinject 
			&& m_pcs->getProcessType() == FiniteElement::HEAT_TRANSPORT 
			&& (result > 1.e-10 || result < -1.e-10)  // use NNNC and connect wells only if storage is active (storing or retrieving)
			)
	{
		//connected_geometry_mode = 1;  // downwind fixed - is set in Read function
		connected_geometry_verbose_level = 0;

		for(size_t i=0; i < upwind_aquifer_mesh_nodes.size(); ++i)
		{
			for(size_t j=0; j < heatExchanger_aquifer_mesh_nodes.size(); ++j)
			{
				double value_connection = value;
				// constant (volumetric) heat capacity
				const double factor = value * heatExchanger_aquifer_mesh_nodes_area_fraction[j] *
								upwind_aquifer_mesh_nodes_area_fraction[i] *
								volumetric_heat_capacity_upwindAquifer *
								fabs(ogs_WDC->get_WellDoubletControl()->get_result().Q_W);

				m_pcs_heat->IncorporateNodeConnectionSourceTerms(
															upwind_aquifer_mesh_nodes[i],
															heatExchanger_aquifer_mesh_nodes[j],
															factor,
															this, value_connection);

				if(connected_geometry_couplingType == 0)
				{  // set RHS values
					m_pcs_heat->add_to_RHS(
#if defined(USE_MPI)
							dom_vector[myrank]->GetDOMNode(heatExchanger_aquifer_mesh_nodes[j]),

#else
							heatExchanger_aquifer_mesh_nodes[j],
#endif
							value_connection);

					double convective_term = -factor / volumetric_heat_capacity_upwindAquifer;
					GetCouplingNODValueConvectiveForm(convective_term, this, heatExchanger_aquifer_mesh_nodes[j]);
				}
			}
		}
	}

	return value * result;
}

// JOD 2019-7-30
double CSourceTerm::apply_contraflow(const double &value, const double& aktuelle_zeit, CRFProcess* m_pcs, double* eqs_rhs)
{
	std::vector<long> nodes_vec = ogs_contraflow->get_nodes_vec();
	stru3::DVec T_s(nodes_vec.size());
	for(std::size_t i=0; i < nodes_vec.size(); ++i)
	{
		T_s[i] = m_pcs->GetNodeValue(nodes_vec[i], 1);
	}

	//std::cout << "\napply contraflow\n";
	contra::Result result;
	if(ogs_contraflow->call_contraflow(aktuelle_zeit, T_s, result))
	{
		std::vector<contra::SegmentData> segment_data_vec = ogs_contraflow->get_segment_data_vec();

		 //for(int i=0; i < nodes_vec.size(); ++i) std::cout << nodes_vec[i] << " ";
		 //std::cout << std::endl;
		 //for(int i=0; i < nodes_vec.size(); ++i) std::cout << result.T_out[i] << " ";
		 //std::cout << std::endl;
		 //std::cout << "Resistances: " << result.resistances_vec[0].R_0_Delta << " " << result.resistances_vec[0].R_1_Delta << std::endl;  // resistances vec due to segmetns
		 //std::cout << "value: " << " " << value << "\n";
		int j = 0;  // segment index
		int k = 0;  // node in segment index
		double L_ele = 0.;
		int N = segment_data_vec[0].N;

		for(std::size_t i=0; i < nodes_vec.size(); ++i)
		{
			if(k == 1)
			{
				L_ele = segment_data_vec[j].L/N;
			}
			else if(k == N)
			{
				k = 0;
				if(j < int(segment_data_vec.size()-1))
				{
					L_ele /=2;  // inside BHE
					++j;
				}
				else
					L_ele = 0.;  // end of BHE
			}

			if(k == 0)
			{
				N = segment_data_vec[j].N;
				L_ele += segment_data_vec[j].L/(2*N);

			}

#ifdef NEW_EQS
			CSparseMatrix* A = NULL;
		   A = m_pcs->get_eqs_new()->get_A();
		   (*A)(nodes_vec[i], nodes_vec[i]) += value * L_ele* ( (1. / result.resistances_vec[j].R_0_Delta) + 1. / result.resistances_vec[j].R_1_Delta);
#else
		   MXInc(nodes_vec[i], nodes_vec[i],  value * L_ele* ( (1. / result.resistances_vec[j].R_0_Delta) + 1. / result.resistances_vec[j].R_1_Delta));
#endif
			// not implemented for PETSc
		   eqs_rhs[nodes_vec[i]] += value * L_ele * ((result.T_in[i] / result.resistances_vec[j].R_0_Delta) + result.T_out[i] / result.resistances_vec[j].R_1_Delta);

		   //double flux = L_ele * (( (result.T_in[i] / result.resistances_vec[j].R_0_Delta) + (result.T_out[i] / result.resistances_vec[j].R_1_Delta)
		//	- ( (1. / result.resistances_vec[j].R_0_Delta) + 1. / result.resistances_vec[j].R_1_Delta)) * T_s[i]);

		   ++k;
		}
	} // end contraflow true (no shutin)

	return 0.;
}


// JOD 2021-12-06
void CSourceTerm::CalculateScalingForNode(const CNodeValue* const cnodev,
		const long& msh_node,  // of local processor in case of MPI
		const CFEMesh* const m_msh, const double& value,
		// the following are updated
		std::vector<scaling_type>& scaling_vec,
		std::map<int, double>& scaling_total_source_term_vector,
		std::map<int, double>& scaling_vec_sum)
{

	scaling_type scaling;

	scaling.node_number = cnodev->msh_node_number;
#if defined(USE_MPI)
	scaling.node_number_local = msh_node;
#endif
	scaling.group_number = GetScalingNodeGroup();
	scaling.keep_value = keep_values;
	scaling.verbosity = scaling_verbosity;

	// calculate a scaling factor for the node
	double scaling_factor = 0;
	std::vector<size_t> elements_connected = m_msh->nod_vector[cnodev->msh_node_number]->getConnectedElementIDs();
	for (long i = 0; i < elements_connected.size(); ++i)
	{
		// permeability
		const CElem* ele = m_msh->ele_vector[elements_connected[i]];
		const double perm = mmp_vector[ele->GetPatchIndex()]->permeability_tensor[0];  // x-direction
		// viscosity - dependent on temperature, nothing else
		CFluidProperties* mfp = MFPGet("LIQUID");
		const CRFProcess* pcs_heat = PCSGet("HEAT_TRANSPORT");
		double temperature = 0;

		const int number_of_nodes = ele->getNodeIndices().Size();
		for (int j = 0; j < number_of_nodes; ++j)
		{
			temperature += pcs_heat->GetNodeValue(ele->GetNodeIndex(j), 1); // implizit
		}
		temperature /= number_of_nodes;

		double prim_values[3];
		prim_values[1] = temperature;

		const double visc = mfp->Viscosity(prim_values);  // use with viscosity model 3
		//
		scaling_factor += perm / visc;
	}

	scaling_factor *=  cnodev->length / elements_connected.size();

	// store the scaling_factor, increment the sum of scalin_factor and source term value for the group (same CSourceTerm instance)
	if(scaling_vec_sum.find(GetScalingNodeGroup()) != scaling_vec_sum.end())
		scaling_vec_sum[GetScalingNodeGroup()] += scaling_factor;
	else
		scaling_vec_sum[GetScalingNodeGroup()] = scaling_factor;

	if(scaling_total_source_term_vector.find(GetScalingNodeGroup()) != scaling_total_source_term_vector.end())
		scaling_total_source_term_vector[GetScalingNodeGroup()] += value;
	else
		scaling_total_source_term_vector[GetScalingNodeGroup()] = value;

	scaling.node_value = scaling_factor;
	scaling_vec.push_back(scaling);
}


// JOD 2022/2 Implementation
void CalculatePeaceman(const CSourceTerm* const m_st, CRFProcess* m_pcs, const long& node_number,
		const std::vector<size_t>& elements_connected, double& factor, double& radius_e)
{

	CElem* face = new CElem(1);
	CFiniteElementStd* fem = new CFiniteElementStd(m_pcs, 21); // 2D: X, Y,  m_pcs->m_msh->GetCoordinateFlag() // coord flag:
	face->SetFace();

	double sum = 0.;
	double sum_ln = 0.;
	factor = 0.;
	bool take_elements_above = false;  // to avoid that face occurs twice
	bool take_elements_below = false;
	int number_of_taken_elements = 0;

	for (size_t j = 0; j < elements_connected.size(); ++j)
	{
		// for r_e
		bool neglect_element = false;

		CElem* elem = m_pcs->m_msh->ele_vector[elements_connected[j]];
					//if (!elem->GetMark())   // !!! do not care about deactivation of elements
					//	continue;
		for(size_t k = 0; k< elem->GetVertexNumber(); ++k)
		{
			if(m_pcs->m_msh->nod_vector[elem->GetNode(k)->GetEquationIndex()]->Z() > m_pcs->m_msh->nod_vector[node_number]->Z() + 1e-5)
			{
				if(take_elements_below)
				{
					neglect_element = true;
					continue;
				}
				else
					take_elements_above = true;
				}

			if(m_pcs->m_msh->nod_vector[elem->GetNode(k)->GetEquationIndex()]->Z() < m_pcs->m_msh->nod_vector[node_number]->Z() - 1e-5)
			{
				if(take_elements_above)
				{
					neglect_element = true;
					continue;
				}
				else
					take_elements_below = true;
            }
		}

		if(neglect_element)
			continue;

		number_of_taken_elements++;

		if(elem->GetDimension() == 2)
		{
			elem->SetOrder(m_pcs->m_msh->getOrder());
			elem->ComputeVolume();
			//fem->ConfigElement(elem, m_pcs->m_num->ele_gauss_points, false);
		}
		else if(elem->GetDimension() == 3)  // pris & hex
		{
			int nfaces, nfn, nodesFace[8];

			CNode* e_node;

			nfaces = elem->GetFacesNumber();
			elem->SetOrder(m_pcs->m_msh->getOrder()); // ?

			int facenumber = -1;
			for (int k = 0; k < nfaces; k++)
			{
				nfn = elem->GetElementFaceNodes(k, nodesFace);
				int counter = 0;
				for (int l = 0; l < nfn; l++)
				{
					e_node = elem->GetNode(nodesFace[l]);
					if(fabs(m_pcs->m_msh->nod_vector[node_number]->Z() - e_node->Z()) < 1.e-5)
						counter++;
				}
				if(counter == nfn)
				{
					facenumber = k;
					break;
				}
			}
			if(facenumber == -1)
				throw std::runtime_error("Error in borehole ST - Element face not found");

			face->SetFace(elem, facenumber);
			face->SetOrder(m_pcs->m_msh->getOrder());
			face->ComputeVolume();
			elem = face;
		}
		else
			throw std::runtime_error("Error in Borehole ST: elements must be 2 or 3D");

		fem->ConfigElement(elem, m_pcs->m_num->ele_gauss_points, false);
		fem->SetMemory();
		fem->SetMaterial();
		fem->CalcLaplace(true);

		Math_Group::Matrix* laplace = fem->get_Laplace();


		if(m_st->getConnectedGeometryVerbosity() > 1)
			laplace->Write();

		CNode* wellNode = NULL;
		int well_ndx;
		for(size_t k = 0; k< elem->GetVertexNumber(); ++k)
			if(elem->GetNode(k)->GetEquationIndex() == node_number)
			{
				wellNode = elem->GetNode(k);
				well_ndx = k;
			}

		for(size_t k = 0; k< elem->GetVertexNumber(); ++k)
		{
			CNode* node = elem->GetNode(k);
			if(node->GetEquationIndex() !=  wellNode->GetEquationIndex() &&
					std::fabs(node->Z() - wellNode->Z()) < 1.e-5)
			{
				const double distance = std::sqrt(
						(wellNode->X() - node->X()) * (wellNode->X() - node->X()) +
						(wellNode->Y() - node->Y()) * (wellNode->Y() - node->Y()));
						// +
						// (wellNode->Z() - node->Z()) * (wellNode->Z() - node->Z()));  // !!! mesh must be horizontal
				const double entry = (*laplace)(well_ndx, k);

				if(m_st->getConnectedGeometryVerbosity() > 1)
					std::cout << "\tNodes " << node_number << " "<< node->GetEquationIndex() << " with distance " <<  distance <<
					" (" << well_ndx<< ", "<< k << "): " <<  entry << std::endl;

				if(distance > 1e-5)
				{
					sum -= entry;
					sum_ln -= entry * std::log(distance);
				}
			}
		}

		if(m_st->get_borehole_modified_aquifer_parameter() == -1)
		{
			/////////////////////////////////////////////////////////////////
			// factor
			const int group = m_pcs->m_msh->ele_vector[elements_connected[j]]->GetPatchIndex();

			switch(m_st->getProcessType())
			{
				case  FiniteElement::LIQUID_FLOW:

					factor += mmp_vector[group]->PermeabilityTensor(group)[0] /  // x-direction
								mfp_vector[0]->Viscosity();  // 1st mfp instance is LIQUID
					break;

				case  FiniteElement::HEAT_TRANSPORT:
				{
					const int dimen = m_pcs->m_msh->GetCoordinateFlag() / 10;
	
					double heat_conductivity_solid[9];
					msp_vector[group]->HeatConductivityTensor(dimen, heat_conductivity_solid, group, j);
					const double porosity = mmp_vector[group]->porosity_model_values[0];  // !!!

					factor += porosity * mfp_vector[0]->HeatConductivity() + // 1st mfp instance is LIQUID
							(1-porosity) * heat_conductivity_solid[0];  // one entry
					break;
				}
				default:
					throw std::runtime_error("Error in Borehole ST - PCS unknown or not supported");
			}
		}


	}  // elements_connected


	if(m_st->get_borehole_modified_aquifer_parameter() != -1)
		factor = 6.283185307179586 * m_st->get_borehole_modified_aquifer_parameter();
	else
		factor *= 6.283185307179586 / number_of_taken_elements;
	

	radius_e = std::exp( (sum_ln - 6.283185307179586) / sum  ) ;  // peaceman

	if(m_st->getConnectedGeometryVerbosity())
		std::cout << "\tPeaceman r_e: " << radius_e  << " at node " << node_number << std::endl;

	delete face;
	delete fem;

}

