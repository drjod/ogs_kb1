/**
 * \file FEM/Output.cpp
 * 05/04/2011 LB Refactoring: Moved from rf_out_new.h
 *
 * Implementation of Output class
 */

// ** INCLUDES **
#include "Configure.h"
#include "Output.h"

#include <fstream>
#include <iostream>
#include <string>

#include <cfloat> // DBL_EPSILON

#include "Configure.h"

#include "FEMIO/GeoIO.h"
#include "GEOObjects.h"
#include "StringTools.h"
#include "fem_ele_std.h"
#include "files0.h"
#include "makros.h"
#include "mathlib.h"
#include "msh_lib.h"
#include "fem_ele.h"
#include "problem.h"
#include "rf_msp_new.h"
#include "rf_pcs.h"
#include "rf_random_walk.h"
#include "rf_tim_new.h"
#include "vtk.h"

// MathLib
#include "MathTools.h"
#include "matrix_class.h" // JOD 2014-11-10
#include "LegacyVtkInterface.h" // JOD 8/2015

#include "mathlib.h"
#include "fem_ele.h"
#include "tools.h"
#include "FileTools.h"

#include "logger.h"
bool flag_block_output_of_initial_values = false;

#include "rf_st_new.h" /// JOD 2018-4-6 to use FaceIntegration (surface-averaged output)

extern size_t max_dim;                            //OK411 todo

#ifdef CHEMAPP
#include "eqlink.h"
#endif

#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL)
#include "par_ddc.h"
#endif
#ifdef SUPERCOMPUTER
// kg44 this is usefull for io-buffering as endl flushes the buffer
#define endl '\n'     // Introduced by WW. LB super bad programming style: this breaks platform independet IO
#define MY_IO_BUFSIZE 4096
#endif // SUPERCOMPUTER
#ifdef GEM_REACT
#include "rf_REACT_GEM.h"
#endif // GEM_REACT

using MeshLib::CFEMesh;
using MeshLib::CElem;
using MeshLib::CEdge;
using MeshLib::CNode;

using namespace std;

COutput::COutput() :
	GeoInfo(GEOLIB::GEODOMAIN), ProcessInfo(), _id(0), out_amplifier(0.0),
	m_msh(NULL), nSteps(-1), _new_file_opened(false), dat_type_name("TECPLOT")
{
	tim_type_name = "TIMES";
	m_pcs = NULL;
	vtk = NULL; //NW
	tecplot_zone_share = false; // 10.2012. WW
	tecplot_datapack_block = 0; // 10.2014 BW
	_ignore_axisymmetry = false;
	VARIABLESHARING = false;	//BG
#if defined(USE_PETSC) || defined(USE_MPI) //|| defined(other parallel libs)//01.3014. WW
	int_disp = 0;
	offset = 0;
	domain_output_counter = 0;
#endif
}

COutput::COutput(size_t id) :
	GeoInfo(GEOLIB::GEODOMAIN), ProcessInfo(), _id(id), out_amplifier(0.0),
	m_msh(NULL), nSteps(-1), _new_file_opened(false), dat_type_name("TECPLOT")
{
	tim_type_name = "TIMES";
	m_pcs = NULL;
	vtk = NULL; //NW
	tecplot_zone_share = false; // 10.2012. WW
	tecplot_datapack_block = 0; // 10.2014 BW
	VARIABLESHARING = false;	//BG
#if defined(USE_PETSC) || defined(USE_MPI) //|| defined(other parallel libs)//01.3014. WW
	int_disp = 0;
	domain_output_counter = 0;
#endif
}
#if defined(USE_PETSC) || defined(USE_MPI) //|| defined(other parallel libs)//03.3012. WW
void COutput::setMPI_Info(const int rank, const int size, std::string rank_str)
{
	mrank = rank;
	msize = size;
	mrank_str = rank_str;
}  
#endif

/*!
   Create the instance of class CVTK
   04.2012. WW
 */
void COutput::CreateVTKInstance(void)
{
#if defined(USE_PETSC) || defined(USE_MPI) //|| defined(other parallel libs)//03.3012. WW
  vtk = new CVTK(mrank, mrank_str);
#else
   vtk = new CVTK();
#endif

}
void COutput::init()
{
	/*if (getProcessType () == FiniteElement::INVALID_PROCESS)
	{
		std::cerr <<
		"COutput::init(): could not initialize process pointer (process type INVALID_PROCESS) and appropriate mesh"
		          << "\n";
		std::cerr <<
		"COutput::init(): trying to fetch process pointer using msh_type_name ... " <<
		"\n";
		if(msh_type_name.size() > 0)
		{
			_pcs = PCSGet(msh_type_name);
			if (_pcs)
				std::cerr << " successful" << "\n";
			else
			{
				std::cerr << " failed" << "\n";
				exit (1);
			}
		}
		else
			std::cerr << " failed" << "\n";
	}
	*/  // removed by JOD 8/2015 since WARNIMG for common case
	m_msh = FEMGet(convertProcessTypeToString(getProcessType()));

    setInternalVarialbeNames(m_msh); //NW

	for (size_t j = 0; j < _nod_value_vector.size(); j++) { // 8/2015 JOD for DELTA_ output
		if (_nod_value_vector[j].find("DELTA") == 0) {
			_pcs = GetPCS(_nod_value_vector[j]);
			_pcs->StoreInitialValues(_nod_value_vector[j]);
		}
	}


    // For binary output of the domain data
#if defined(USE_PETSC) // || defined(other solver libs)//01.3014. WW
    if( (getGeoType() == GEOLIB::GEODOMAIN) || (dat_type_name.compare("BINARY") != 0 ) )
    {
      //dat_type_name = "BINARY";
       setDataArrayDisp();
    }  
#endif

}

COutput::~COutput()
{
	mmp_value_vector.clear();             //OK

	if (this->vtk != NULL)
		delete vtk;               //NW
}

const std::string& COutput::getGeoName () const
{
	return geo_name;
}

/**************************************************************************
   FEMLib-Method:
   Task: OUT read function
   Programing:
   06/2004 OK Implementation
   07/2004 WW Remove old files
   11/2004 OK string streaming by SB for lines
   03/2005 OK PCS_TYPE
   12/2005 OK DIS_TYPE
   12/2005 OK MSH_TYPE
   08/2008 OK MAT
   06/2010 TF formated, restructured, signature changed, use new GEOLIB data structures
   09/2010 TF signature changed, removed some variables
**************************************************************************/
ios::pos_type COutput::Read(std::ifstream& in_str,
                            const GEOLIB::GEOObjects& geo_obj, const std::string& unique_geo_name)
{
	std::string line_string;
	bool new_keyword = false;
	ios::pos_type position;
	bool new_subkeyword = false;
	std::string tec_file_name;
	ios::pos_type position_line;
	bool ok = true;
	std::stringstream in;
	string name;
	ios::pos_type position_subkeyword;

	// Schleife ueber alle Phasen bzw. Komponenten
	while (!new_keyword)
	{
		position = in_str.tellg();
		if (new_subkeyword)
			in_str.seekg(position_subkeyword, ios::beg);
		new_subkeyword = false;
		// SB new input		in_str.getline(buffer,MAX_ZEILE);
		// SB new input         line_string = buffer;
		line_string.clear();
		line_string = GetLineFromFile1(&in_str);
		if (line_string.size() < 1)
			break;

		if (Keyword(line_string))
			return position;

		// subkeyword found
		if (line_string.find("$NOD_VALUES") != string::npos)
		{
			while ((!new_keyword) && (!new_subkeyword))
			{
				position_subkeyword = in_str.tellg();
				//SB input with comments  in_str >> line_string>>ws;
				line_string = GetLineFromFile1(&in_str);
				if (line_string.find("#") != string::npos)
					return position;
				if (line_string.find("$") != string::npos)
				{
					new_subkeyword = true;
					break;
				}
				if (line_string.size() == 0)
					break;  //SB: empty line
				in.str(line_string);
				in >> name;
				//_alias_nod_value_vector.push_back(name);
                _nod_value_vector.push_back(name);
				in.clear();
			}

			continue;
		}
		//--------------------------------------------------------------------
		// subkeyword found //MX
		if (line_string.find("$PCON_VALUES") != string::npos)
		{
			while ((!new_keyword) && (!new_subkeyword))
			{
				position_subkeyword = in_str.tellg();
				in_str >> line_string >> ws;
				if (line_string.find("#") != string::npos)
					return position;
				if (line_string.find("$") != string::npos)
				{
					new_subkeyword = true;
					break;
				}
				if (line_string.size() == 0)
					break;
				_pcon_value_vector.push_back(line_string);
			}
			continue;
		}

		//--------------------------------------------------------------------
		// subkeyword found
		if (line_string.find("$ELE_VALUES") != string::npos)
		{
			ok = true;
			while (ok)
			{
				position_line = in_str.tellg();
				in_str >> line_string;
				if (SubKeyword(line_string))
				{
					in_str.seekg(position_line, ios::beg);
					ok = false;
					continue;
				}
				if (Keyword(line_string))
					return position;
				_ele_value_vector.push_back(line_string);
				in_str.ignore(MAX_ZEILE, '\n');
			}
			/*
			   // Commented by WW
			   // Remove files
			   tec_file_name = file_base_name + "_domain_ele" + TEC_FILE_EXTENSION;
			   remove(tec_file_name.c_str());
			 */
			continue;
		}
		//-------------------------------------------------------------------- // Added 03.2010 JTARON
		// subkeyword found
		if (line_string.find("$RWPT_VALUES") != string::npos)
		{
			while ((!new_keyword) && (!new_subkeyword))
			{
				position_subkeyword = in_str.tellg();
				line_string = GetLineFromFile1(&in_str);
				if (line_string.find("#") != string::npos)
					return position;
				if (line_string.find("$") != string::npos)
				{
					new_subkeyword = true;
					break;
				}
				if (line_string.size() == 0)
					break;  //SB: empty line
				in.str(line_string);
				in >> name;
				_rwpt_value_vector.push_back(name);
				in.clear();
			}
			continue;
		}

		//subkeyword found
		if (line_string.find("$GEO_TYPE") != string::npos)
		{
			FileIO::GeoIO::readGeoInfo (this, in_str, geo_name, geo_obj, unique_geo_name);
			continue;
		}

		// subkeyword found
		if (line_string.find("$TIM_TYPE") != string::npos)
		{
			while ((!new_keyword) && (!new_subkeyword))
			{
				position_subkeyword = in_str.tellg();
				in_str >> line_string;
				if (line_string.size() == 0) //SB
					break;
				if (line_string.find("#") != string::npos)
				{
					new_keyword = true;
					break;
				}
				if (line_string.find("$") != string::npos)
				{
					new_subkeyword = true;
					break;
				}
				if (line_string.find("STEPS") != string::npos)
				{
					in_str >> nSteps;
					tim_type_name = "STEPS"; //OK
					break; //kg44 I guess that was missing..otherwise it pushes back a time_vector!
				}
				// JT 2010, reconfigured (and added RWPT)... didn't work
				if (line_string.find("STEPPING") != string::npos)
				{
					double stepping_length, stepping_end, stepping_current;
					in_str >> stepping_length >> stepping_end;
					stepping_current = stepping_length;
					while (stepping_current <= stepping_end)
					{
						time_vector.push_back(stepping_current);
						//						rwpt_time_vector.push_back(stepping_current);
						stepping_current += stepping_length;
					}
				}
				else
					time_vector.push_back(strtod(line_string.data(), NULL));
				//					rwpt_time_vector.push_back(strtod(line_string.data(), NULL));
				in_str.ignore(MAX_ZEILE, '\n');
			}
			continue;
		}

		// subkeyword found
		if (line_string.find("$DAT_TYPE") != string::npos)
		{
			in_str >> dat_type_name;
			if (dat_type_name == "CONTENT") // JOD 2/2015
			{
				in_str >> mmp_index;
			}
			if (dat_type_name == "VOLUME") // JOD 2020-1-16
			{
						mmp_index = -2;  // flag for volume calculation
						in_str >> domainIntegration_lowerThreshold >> domainIntegration_upperThreshold; // JOD 2020-1-15
			}
			in_str.ignore(MAX_ZEILE, '\n');
			continue;
		}

		// Coordinates of each node as well as connection list is stored only for the first time step; BG: 05/2011
        if (line_string.find("$VARIABLESHARING") != string::npos)
        {
	       this->VARIABLESHARING = true;
		   continue;
        }

		// subkeyword found
		if (line_string.find("$AMPLIFIER") != string::npos)
		{
			in_str >> out_amplifier;
			in_str.ignore(MAX_ZEILE, '\n');
			continue;
		}

		// subkeyword found
		if (line_string.find("$PCS_TYPE") != string::npos)
		{
			std::string tmp_pcs_type_name;
			in_str >> tmp_pcs_type_name;
			setProcessType(FiniteElement::convertProcessType(tmp_pcs_type_name));
			in_str.ignore(MAX_ZEILE, '\n');
			/* // Comment by WW
			   // Remove files
			   tec_file_name = pcs_type_name + "_" + "domain" + "_line" + TEC_FILE_EXTENSION;
			   remove(tec_file_name.c_str());
			   tec_file_name = pcs_type_name + "_" + "domain" + "_tri" + TEC_FILE_EXTENSION;
			   remove(tec_file_name.c_str());
			   tec_file_name = pcs_type_name + "_" + "domain" + "_quad" + TEC_FILE_EXTENSION;
			   remove(tec_file_name.c_str());
			   tec_file_name = pcs_type_name + "_" + "domain" + "_tet" + TEC_FILE_EXTENSION;
			   remove(tec_file_name.c_str());
			   tec_file_name = pcs_type_name + "_" + "domain" + "_pris" + TEC_FILE_EXTENSION;
			   remove(tec_file_name.c_str());
			   tec_file_name = pcs_type_name + "_" + "domain" + "_hex" + TEC_FILE_EXTENSION;
			   remove(tec_file_name.c_str());
			 */
			continue;
		}

		// subkeyword found
		if (line_string.find("$DIS_TYPE") != string::npos)
		{
			std::string dis_type_name;
			in_str >> dis_type_name;
			setProcessDistributionType (FiniteElement::convertDisType(dis_type_name));
			in_str.ignore(MAX_ZEILE, '\n');
			continue;
		}

		// subkeyword found
		if (line_string.find("$MSH_TYPE") != string::npos)
		{
			in_str >> msh_type_name;
			in_str.ignore(MAX_ZEILE, '\n');
			continue;
		}

		//OK
		if (line_string.find("$MMP_VALUES") != string::npos)
		{
			ok = true;
			while (ok)
			{
				position_line = in_str.tellg();
				in_str >> line_string;
				if (SubKeyword(line_string))
				{
					in_str.seekg(position_line, ios::beg);
					ok = false;
					continue;
				}
				if (Keyword(line_string))
					return position;
				mmp_value_vector.push_back(line_string);
				in_str.ignore(MAX_ZEILE, '\n');
			}
			continue;
		}

		//OK
		if (line_string.find("$MFP_VALUES") != string::npos)
		{
			ok = true;
			while (ok)
			{
				position_line = in_str.tellg();
				in_str >> line_string;
				if (SubKeyword(line_string))
				{
					in_str.seekg(position_line, ios::beg);
					ok = false;
					continue;
				}
				if (Keyword(line_string))
					return position;
				mfp_value_vector.push_back(line_string);
				in_str.ignore(MAX_ZEILE, '\n');
			}

			continue;
		}
		// For teplot zone share. 10.2012. WW
		if (line_string.find("$TECPLOT_ZONE_SHARE") != string::npos)
		{
			tecplot_zone_share = true; 
			continue;
		}
		if (line_string.find("$IGNORE_AXISYMMETRY") != string::npos)  // JOD 2018-08-17  for content output
		{
			_ignore_axisymmetry = true;
			std::cout << "\tIgnore axisymmetry in content output " << mmp_index << '\n';
			continue;
		}
		if (line_string.find("$APPEND_DATA") != string::npos)  // JOD 2020-05-08
		{
			_new_file_opened = true;  // declared as already deleted
			flag_block_output_of_initial_values = true;  // no output for step 0
			std::cout << "\tAppend data (File not deleted)\n";
			logger.block_deletion();
			continue;
		}
		// For tecplot block datapack format 10.2014 BW
		if (line_string.find("$TECPLOT_DATAPACK_BLOCK") != string::npos)
		{
			in_str >> tecplot_datapack_block;
			if (tecplot_datapack_block == 1)
				std::cout << "Warning: Only Block Data Packing Type of Tecplot Output Will be Genereated!!" << '\n';
			in_str.ignore(MAX_ZEILE, '\n');
			continue;
		}
	}
	return position;
}

/**************************************************************************
   FEMLib-Method:
   Task: write function
   Programing:
   06/2004 OK Implementation
   01/2005 OK Extensions
   12/2005 OK DIS_TYPE
   12/2005 OK MSH_TYPE
   12/2008 NW DAT_TYPE
   05/2009 OK bug fix STEPS
**************************************************************************/
void COutput::Write(fstream* out_file)
{
	//--------------------------------------------------------------------
	// KEYWORD
	*out_file << "#OUTPUT" << "\n";
	//--------------------------------------------------------------------
	// PCS_TYPE
	*out_file << " $PCS_TYPE" << "\n" << "  ";
	*out_file << convertProcessTypeToString(getProcessType()) << "\n";
	//--------------------------------------------------------------------
	// NOD_VALUES
	*out_file << " $NOD_VALUES" << "\n";
	size_t nod_value_vector_size(_nod_value_vector.size());
	for (size_t i = 0; i < nod_value_vector_size; i++)
		*out_file << "  " << _nod_value_vector[i] << "\n";
	//--------------------------------------------------------------------
	// ELE_VALUES
	*out_file << " $ELE_VALUES" << "\n";
	size_t ele_value_vector_size (_ele_value_vector.size());
	for (size_t i = 0; i < ele_value_vector_size; i++)
		*out_file << "  " << _ele_value_vector[i] << "\n";
	//--------------------------------------------------------------------
	// GEO_TYPE
	*out_file << " $GEO_TYPE" << "\n";
	*out_file << "  ";
	*out_file << getGeoTypeAsString() << " " << geo_name << "\n";
	//--------------------------------------------------------------------
	// TIM_TYPE
	*out_file << " $TIM_TYPE" << "\n";
	if (tim_type_name == "STEPS")
		*out_file << "  " << tim_type_name << " " << nSteps << "\n";
	else
	{
		size_t time_vector_size (time_vector.size());
		for (size_t i = 0; i < time_vector_size; i++)
			*out_file << "  " << time_vector[i] << "\n";
	}

	// DIS_TYPE
	//	if (_dis_type_name.size() > 0) {
	//		*out_file << " $DIS_TYPE" << "\n";
	//		*out_file << "  ";
	//		*out_file << _dis_type_name << "\n";
	//	}
	if (getProcessDistributionType() != FiniteElement::INVALID_DIS_TYPE)
	{
		*out_file << " $DIS_TYPE" << "\n";
		*out_file << "  ";
		*out_file << convertDisTypeToString (getProcessDistributionType()) << "\n";
	}

	// MSH_TYPE
	if (msh_type_name.size() > 0)
	{
		*out_file << " $MSH_TYPE" << "\n";
		*out_file << "  ";
		*out_file << msh_type_name << "\n";
	}
	//--------------------------------------------------------------------
	// DAT_TYPE
	*out_file << " $DAT_TYPE" << "\n";
	*out_file << "  ";
	*out_file << dat_type_name << "\n";
	//--------------------------------------------------------------------
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   08/2004 OK Implementation
   03/2005 OK MultiMSH
   08/2005 WW Changes for MultiMSH
   12/2005 OK VAR,MSH,PCS concept
   07/2007 NW Multi Mesh Type
**************************************************************************/
void COutput::NODWriteDOMDataTEC()
{
	int te = 0;
	string eleType;
	string tec_file_name;
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL)
	char tf_name[10];
	std::cout << "Process " << myrank << " in WriteDOMDataTEC" << "\n";
#endif
	//----------------------------------------------------------------------
	// Tests
	//OK4704
	if((_nod_value_vector.size() == 0) && (mfp_value_vector.size() == 0))
		return;
	//......................................................................
	// MSH
	//m_msh = FEMGet(pcs_type_name);
	//  m_msh = GetMSH();
	if(!m_msh)
	{
		cout << "Warning in COutput::NODWriteDOMDataTEC() - no MSH data" << "\n";
		return;
	}
	//======================================================================
	vector<int> mesh_type_list;           //NW
	if (m_msh->getNumberOfLines() > 0)
		mesh_type_list.push_back(1);
	if (m_msh->getNumberOfQuads() > 0)
		mesh_type_list.push_back(2);
	if (m_msh->getNumberOfHexs () > 0)
		mesh_type_list.push_back(3);
	if (m_msh->getNumberOfTris () > 0)
		mesh_type_list.push_back(4);
	if (m_msh->getNumberOfTets () > 0)
		mesh_type_list.push_back(5);
	if (m_msh->getNumberOfPrisms () > 0)
		mesh_type_list.push_back(6);
	if (m_msh->getNumberOfPyramids() > 0)
		mesh_type_list.push_back(7);

	// Output files for each mesh type
	//NW
	for (int i = 0; i < (int)mesh_type_list.size(); i++)
	{
		te = mesh_type_list[i];
		//----------------------------------------------------------------------
		// File name handling
		tec_file_name = file_base_name + "_" + "domain";
		if(msh_type_name.size() > 0) // MultiMSH
			tec_file_name += "_" + msh_type_name;
		if(getProcessType() != FiniteElement::INVALID_PROCESS) // PCS
			tec_file_name += "_" + convertProcessTypeToString(getProcessType());
		//======================================================================
		switch (te)               //NW
		{
		case 1:
			tec_file_name += "_line";
			eleType = "QUADRILATERAL";
			break;
		case 2:
			tec_file_name += "_quad";
			eleType = "QUADRILATERAL";
			break;
		case 3:
			tec_file_name += "_hex";
			eleType = "BRICK";
			break;
		case 4:
			tec_file_name += "_tri";
			eleType = "QUADRILATERAL";
			break;
		case 5:
			tec_file_name += "_tet";
			eleType = "TETRAHEDRON";
			break;
		case 6:
			tec_file_name += "_pris";
			eleType = "BRICK";
			break;
		case 7:
			tec_file_name += "_pyra";
			eleType = "BRICK";
			break;
		}
		/*
		   if(m_msh->msh_no_line>0)
		   {
		      tec_file_name += "_line";
		      eleType = "QUADRILATERAL";
		     te=1;
		   }
		   else if (m_msh->msh_no_quad>0)
		   {
		      tec_file_name += "_quad";
		      eleType = "QUADRILATERAL";
		   te=2;
		   }
		   else if (m_msh->msh_no_hexs>0)
		   {
		   tec_file_name += "_hex";
		   eleType = "BRICK";
		   te=3;
		   }
		   else if (m_msh->msh_no_tris>0)
		   {
		   tec_file_name += "_tri";
		   //???Who was this eleType = "TRIANGLE";
		   eleType = "QUADRILATERAL";
		   te=4;
		   }
		   else if (m_msh->msh_no_tets>0)
		   {
		   tec_file_name += "_tet";
		   eleType = "TETRAHEDRON";
		   te=5;
		   }
		   else if (m_msh->msh_no_pris>0)
		   {
		   tec_file_name += "_pris";
		   eleType = "BRICK";
		   te=6;
		   }
		 */
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL)
		sprintf(tf_name, "%d", myrank);
		tec_file_name += "_" + string(tf_name);
		std::cout << "Tecplot filename: " << tec_file_name << "\n";
#endif
#if defined(USE_PETSC)  //|| defined(other parallel libs)//03.3012. WW
		tec_file_name += "_"+mrank_str;
		std::cout << "Tecplot filename: " << tec_file_name << "\n";
#endif
		tec_file_name += TEC_FILE_EXTENSION;
		//WW
		if (!_new_file_opened)
			remove(tec_file_name.c_str());

		fstream tec_file (tec_file_name.data(),ios::app | ios::out);
		tec_file.setf(ios::scientific,ios::floatfield);
		tec_file.precision(12);
		if (!tec_file.good())
			return;
#ifdef SUPERCOMPUTER
		// kg44 buffer the output
		char mybuf1 [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
		tec_file.rdbuf()->pubsetbuf(mybuf1,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
#endif
		//
		WriteTECHeader(tec_file,te,eleType);
		WriteTECNodeData(tec_file);

		// 08.2012. WW
        if(tecplot_zone_share)
		{
	     	if(!_new_file_opened)  
	           WriteTECElementData(tec_file,te);
		}
		else
		{
            WriteTECElementData(tec_file,te);
		}
        _new_file_opened = true;
		tec_file.close();         // kg44 close file
		//--------------------------------------------------------------------
		// tri elements
		// ***** 07/2010 TF commented out block since the global variable is always zero
		//    if(msh_no_tris>0){
		//    //string tec_file_name = pcs_type_name + "_" + "domain" + "_tri" + TEC_FILE_EXTENSION;
		//#ifdef SUPERCOMPUTER
		//// buffer the output
		//      char sxbuf1[MY_IO_BUFSIZE*MY_IO_BUFSIZE];
		//#endif
		//      string tec_file_name = file_base_name + "_" + "domain" + "_tri" + TEC_FILE_EXTENSION;
		//      fstream tec_file1 (tec_file_name.data(),ios::app|ios::out);
		//      tec_file1.setf(ios::scientific,ios::floatfield);
		//      tec_file1.precision(12);
		//      if (!tec_file1.good()) return;
		//#ifdef SUPERCOMPUTER
		//      tec_file1.rdbuf()->pubsetbuf(sxbuf1,MY_IO_BUFSIZE*MY_IO_BUFSIZE);
		//#endif
		//      //OK  tec_file1.clear();
		//      //OK  tec_file1.seekg(0L,ios::beg);
		//      WriteTECHeader(tec_file1,4,"TRIANGLE");
		//      WriteTECNodeData(tec_file1);
		//      WriteTECElementData(tec_file1,4);
		//      tec_file1.close(); // kg44 close file
		//    }
		//--------------------------------------------------------------------
		// quad elements
		// ***** 07/2010 TF commented out block since the global variable is always zero
		//    if(msh_no_quad>0){
		//      //string tec_file_name = pcs_type_name + "_" + "domain" + "_quad" + TEC_FILE_EXTENSION;
		//#ifdef SUPERCOMPUTER
		//      char sxbuf2[MY_IO_BUFSIZE*MY_IO_BUFSIZE];
		//#endif
		//      string tec_file_name = file_base_name + "_" + "domain" + "_quad" + TEC_FILE_EXTENSION;
		//      fstream tec_file (tec_file_name.data(),ios::app|ios::out);
		//
		//      tec_file.setf(ios::scientific,ios::floatfield);
		//      tec_file.precision(12);
		//      if (!tec_file.good()) return;
		//#ifdef SUPERCOMPUTER
		//      tec_file.rdbuf()->pubsetbuf(sxbuf2,MY_IO_BUFSIZE*MY_IO_BUFSIZE);
		//#endif
		//      WriteTECHeader(tec_file,2,"QUADRILATERAL");
		//      WriteTECNodeData(tec_file);
		//      WriteTECElementData(tec_file,2);
		//      tec_file.close(); // kg44 close file
		//    }
		//--------------------------------------------------------------------
		// tet elements
		// ***** 07/2010 TF commented out block since the global variable is always zero
		//    if(msh_no_tets>0){
		//      //string tec_file_name = pcs_type_name + "_" + "domain" + "_tet" + TEC_FILE_EXTENSION;
		//#ifdef SUPERCOMPUTER
		//      char sxbuf3[MY_IO_BUFSIZE*MY_IO_BUFSIZE];
		//#endif
		//
		//      string tec_file_name = file_base_name + "_" + "domain" + "_tet";
		//
		//#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL)
		//      sprintf(tf_name, "%d", myrank);
		//      tec_file_name += "_" + string(tf_name);
		//#endif
		//
		//      tec_file_name += TEC_FILE_EXTENSION;
		//
		//      fstream tec_file (tec_file_name.data(),ios::app|ios::out);
		//
		//      tec_file.setf(ios::scientific,ios::floatfield);
		//      tec_file.precision(12);
		//      if (!tec_file.good()) return;
		//#ifdef SUPERCOMPUTER
		//      tec_file.rdbuf()->pubsetbuf(sxbuf3,MY_IO_BUFSIZE*MY_IO_BUFSIZE);
		//#endif
		//
		//      WriteTECHeader(tec_file,5,"TETRAHEDRON");
		//      WriteTECNodeData(tec_file);
		//      WriteTECElementData(tec_file,5);
		//      tec_file.close(); // kg44 close file
		//    }
		//--------------------------------------------------------------------
		//    // pris elements
		// ***** 07/2010 TF commented out block since the global variable is always zero
		//    if(msh_no_pris>0){
		//      //string tec_file_name = pcs_type_name + "_" + "domain" + "_pris" + TEC_FILE_EXTENSION;
		//#ifdef SUPERCOMPUTER
		//        char sxbuf4[MY_IO_BUFSIZE*MY_IO_BUFSIZE];
		//#endif
		//      string tec_file_name = file_base_name + "_" + "domain" + "_pris" + TEC_FILE_EXTENSION;
		//      fstream tec_file (tec_file_name.data(),ios::app|ios::out);
		//
		//      tec_file.setf(ios::scientific,ios::floatfield);
		//      tec_file.precision(12);
		//      if (!tec_file.good()) return;
		//#ifdef SUPERCOMPUTER
		//      tec_file.rdbuf()->pubsetbuf(sxbuf4,MY_IO_BUFSIZE*MY_IO_BUFSIZE);
		//#endif
		//
		//      WriteTECHeader(tec_file,6,"BRICK");
		//      WriteTECNodeData(tec_file);
		//      WriteTECElementData(tec_file,6);
		//      tec_file.close(); // kg44 close file
		//    }
		//--------------------------------------------------------------------
		// hex elements
		// ***** 07/2010 TF commented out block since the global variable is always zero
		//    if(msh_no_hexs>0){
		//      //string tec_file_name = pcs_type_name + "_" + "domain" + "_hex" + TEC_FILE_EXTENSION;
		//#ifdef SUPERCOMPUTER
		//        char sxbuf5[MY_IO_BUFSIZE*MY_IO_BUFSIZE];
		//#endif
		//
		//      string tec_file_name = file_base_name + "_" + "domain" + "_hex" + TEC_FILE_EXTENSION;
		//      fstream tec_file (tec_file_name.data(),ios::app|ios::out);
		//
		//
		//      tec_file.setf(ios::scientific,ios::floatfield);
		//      tec_file.precision(12);
		//      if (!tec_file.good()) return;
		//#ifdef SUPERCOMPUTER
		//      tec_file.rdbuf()->pubsetbuf(sxbuf5,MY_IO_BUFSIZE*MY_IO_BUFSIZE);
		//#endif
		//      WriteTECHeader(tec_file,3,"BRICK");
		//      WriteTECNodeData(tec_file);
		//      WriteTECElementData(tec_file,3);
		//      tec_file.close(); // kg44 close file
		//    }
	}
}

/*
   Set data array displacement for parallel output
  
    WW  12.2013
*/
#if defined(USE_PETSC) // || defined(other solver libs)//01.2014. WW
void COutput::setDataArrayDisp()
{
   //   MPI_Barrier (MPI_COMM_WORLD);
   //
   int *i_cnt;
   int *i_disp;
   int *i_recv;

   i_cnt =  new int[msize];
   i_disp = new int[msize];
   i_recv = new int[msize];

   for(int i=0; i<msize; i++)
   {
     i_cnt[i] = 1; 
     i_disp[i] = i;
   }

   int size_local =  fem_msh_vector[0]->getNumNodesLocal();

   MPI_Allgatherv(&size_local, 1, MPI_INT, i_recv, i_cnt, i_disp,
                   MPI_INT, MPI_COMM_WORLD);  
  
   int_disp = 0;
   for(int i=0; i<mrank; i++)
   {
       int_disp += i_recv[i]; 
   }

   delete [] i_cnt;
   delete [] i_disp;
   delete [] i_recv;
   //   MPI_Barrier (MPI_COMM_WORLD);
} 

/*
    Write variable informations of the domain
  
    WW  12.2013
*/
void COutput::NODDomainWriteBinary_Header()
{

   if(mrank != 0)
      return;

   if( dat_type_name.compare("BINARY") != 0 )
     return;

   string file_name;

   file_name = file_base_name +  "_" + convertProcessTypeToString(getProcessType()) + "_domain_" + "node_value_header.txt";
   std::cout << "Name of the header file: " << file_name << "\n";


   ofstream os (file_name.data(), ios::trunc | ios::out);
   if (!os.good())
   {
	  return;
   }

   os << msize << "\n";
  
   m_pcs =  GetPCS();

   os  << domain_output_counter  <<  "\n";
 
   const size_t num_prim_unknowns = m_pcs->GetPrimaryVNumber();
   const size_t num_2nd_unknowns = m_pcs->GetSecondaryVNumber();

   os <<  convertProcessTypeToString(getProcessType()) << "\n";
   os << num_prim_unknowns + num_2nd_unknowns <<  "\n";

   for(size_t i=0; i < num_prim_unknowns; i++)
   {
      os <<  m_pcs->GetPrimaryVName(i) << " ";
   } 
   for(size_t i=0; i < num_2nd_unknowns; i++)
   {
       os << m_pcs->GetSecondaryVName(i) << " ";
   } 
   os << "\n";
   
   // Write number of unknowns
   os << m_pcs->m_msh->getNumNodesGlobal()  << "\n";

   os.close();
}

/*
  WW 08.2013
*/
void COutput::NODDomainWriteBinary()
{
   string file_name;

   file_name = file_base_name + "_" + convertProcessTypeToString(getProcessType()) + "_domain_variables" + ".bin";
   std::cout << "Name of the binary file for node and element data: " << file_name << "\n";

   domain_output_counter++;  

   if(!_new_file_opened)
   {
      remove(file_name.c_str());      
   }

   m_pcs =  GetPCS();
  
   MPI_Barrier (MPI_COMM_WORLD);

   MPI_Offset offset_new;
   MPI_File fh;
   int rc = 0;
   
   if(!_new_file_opened)
   {
      rc = MPI_File_open(MPI_COMM_WORLD, &file_name[0], MPI_MODE_WRONLY | MPI_MODE_CREATE,  MPI_INFO_NULL, &fh);
      offset = 0;
   }
   else
   {
       rc = MPI_File_open(MPI_COMM_WORLD, &file_name[0], MPI_MODE_WRONLY | MPI_MODE_APPEND,  MPI_INFO_NULL, &fh);
   }
   
   if (rc ) 
   {	   
       MPI_Finalize();
       cout<<"Cannot open "<<file_name<<"does not exist." <<"\n";
       exit(0);
   }
 
   //MPI_File_get_position( fh, &offset ); 
   // Write time and remember the number of processes#
   string ftype = "native";
  
   offset_new = offset + mrank*sizeof(double); 
   MPI_File_set_view(fh, offset_new, MPI_DOUBLE, MPI_DOUBLE,  &ftype[0], MPI_INFO_NULL);
   MPI_File_write(fh, &_time, 1, MPI_DOUBLE, MPI_STATUS_IGNORE); //_all
   offset += msize*sizeof(double); 

   const size_t num_prim_unknowns = m_pcs->GetPrimaryVNumber();
   const size_t num_2nd_unknowns = m_pcs->GetSecondaryVNumber();
   // Write unknowns
   size_t n_unknowns = 0; 
   n_unknowns = m_pcs->m_msh->getNumNodesLocal();
   const int nn = m_pcs->m_msh->getNumNodesGlobal();

   // Write primary unknowns
   for(size_t i=0; i < num_prim_unknowns; i++)
   {
      double *node_values = m_pcs->getNodeValue_per_Variable(2*i + 1);
      offset_new = offset + int_disp*sizeof(double);
      MPI_File_set_view(fh, offset_new, MPI_DOUBLE, MPI_DOUBLE,  &ftype[0], MPI_INFO_NULL);
      MPI_File_write(fh, node_values, n_unknowns, MPI_DOUBLE, MPI_STATUS_IGNORE); //_all
      offset += nn * sizeof(double);  
   } 

   // Write secondary unknowns
   for(size_t i=0; i < num_2nd_unknowns; i++)
   {
       double *node_values = m_pcs->getNodeValue_per_Variable(2*num_prim_unknowns + i);
       offset_new = offset + int_disp*sizeof(double);
       MPI_File_set_view(fh, offset_new, MPI_DOUBLE, MPI_DOUBLE,  &ftype[0], MPI_INFO_NULL);
       MPI_File_write(fh, node_values, n_unknowns, MPI_DOUBLE, MPI_STATUS_IGNORE); //_all
       offset += nn * sizeof(double);  
   } 

   MPI_File_sync( fh ) ; 
   MPI_Barrier( MPI_COMM_WORLD ) ;
   MPI_File_sync( fh ) ; 
   MPI_File_close(&fh);
   _new_file_opened = true;
}
#endif //  end of USE_PETSC

/**************************************************************************
   FEMLib-Method:
   Programing:
   08/2004 OK Implementation
   08/2004 WW Output node variables by their names given in .out file
   03/2005 OK MultiMSH
   08/2005 WW Correction of node index
   12/2005 OK Mass transport specifics
   OK ??? too many specifics
**************************************************************************/
void COutput::WriteTECNodeData(fstream &tec_file)
{
	const size_t nName(_nod_value_vector.size());
	double val_n = 0.;                    //WW
	int nidx, nidx_dm[3];
	vector<int> NodeIndex(nName);
	string nod_value_name;                //OK
	CNode *node = NULL;
	CRFProcess *deform_pcs = NULL; // 23.01.2012. WW. nulltpr

	int timelevel;
	//	m_msh = GetMSH();
	CRFProcess* m_pcs_out = NULL;
	// MSH
	for (size_t k = 0; k < nName; k++)
	{
		m_pcs = PCSGet(_nod_value_vector[k], true);
		if (m_pcs != NULL)
		{
            NodeIndex[k] = m_pcs->GetNodeValueIndex(_nod_value_vector[k],true); // JT Latest.
            if(     (m_pcs->getProcessType() == FiniteElement::DEFORMATION)
                 || (m_pcs->getProcessType() == FiniteElement::DEFORMATION_DYNAMIC)
                 ||  (m_pcs->getProcessType() == FiniteElement::DEFORMATION_FLOW)
                 ||  (m_pcs->getProcessType() == FiniteElement::DEFORMATION_H2)			
               )
            {
                deform_pcs = m_pcs;
            }
         }
    }

	if (deform_pcs)    // 23.01.2012. WW.
	{
		nidx_dm[0] = deform_pcs->GetNodeValueIndex("DISPLACEMENT_X1") + 1;
		nidx_dm[1] = deform_pcs->GetNodeValueIndex("DISPLACEMENT_Y1") + 1;
		if (max_dim > 1)
			nidx_dm[2] = deform_pcs->GetNodeValueIndex("DISPLACEMENT_Z1") + 1;
		else
			nidx_dm[2] = -1;
	}
	// 08.2012. WW
	bool out_coord = true;
	if(tecplot_zone_share && _new_file_opened)
       out_coord = false;
	for (size_t j = 0; j < m_msh->GetNodesNumber(false); j++)
	{
       node = m_msh->nod_vector[j];  // 23.01.2013. WW
       const size_t n_id = node->GetIndex(); 

	   if(out_coord) // 08.2012. WW
	   {
		  // XYZ
		  const double *x = node->getData(); // 23.01.2013. WW

		  // Amplifying DISPLACEMENTs
	      if (deform_pcs)  // 23.01.2012. WW.
		  {
	         for (size_t i = 0; i <max_dim+1; i++)
		        tec_file << x[i] + out_amplifier * m_pcs->GetNodeValue(n_id, nidx_dm[i]) << " ";
	         for (size_t i = max_dim+1; i<3; i++)
		        tec_file << x[i] << " ";
		  }
		  else
		  {
	         for (size_t i = 0; i < 3; i++)
		        tec_file << x[i] << " ";
		  }
	   }  
		// NOD values
		// Mass transport
		//     if(pcs_type_name.compare("MASS_TRANSPORT")==0){
		if (getProcessType() == FiniteElement::MASS_TRANSPORT)
			for (size_t i = 0; i < _nod_value_vector.size(); i++)
			{
				std::string nod_value_name = _nod_value_vector[i];
				for (size_t l = 0; l < pcs_vector.size(); l++)
				{
					m_pcs = pcs_vector[l];
					//					if (m_pcs->pcs_type_name.compare("MASS_TRANSPORT") == 0) {
					if (m_pcs->getProcessType () == FiniteElement::MASS_TRANSPORT)
					{
						timelevel = 0;
						for (size_t m = 0;
						     m < m_pcs->nod_val_name_vector.size(); m++)
							if (m_pcs->nod_val_name_vector[m].compare(
							            nod_value_name) == 0)
							{
								m_pcs_out = PCSGet(FiniteElement::MASS_TRANSPORT,
								                   nod_value_name);
								if (!m_pcs_out)
									continue;
								if (timelevel == 1)
								{
									nidx =
									        m_pcs_out->
									        GetNodeValueIndex(
									                nod_value_name)
									        +
									        timelevel;
									tec_file <<
									m_pcs_out->GetNodeValue(n_id, nidx) << " ";
								}
								timelevel++;
							}
					}
				}
			}
		else
		{
			for (size_t k = 0; k < nName; k++)
			{
				m_pcs = GetPCS(_nod_value_vector[k]);
				if (m_pcs != NULL) { //WW

					if (NodeIndex[k] > -1) {
						if (_nod_value_vector[k].find("DELTA") == 0) // JOD 2014-11-10 
							val_n = m_pcs->GetNodeValue(n_id, 1) - m_pcs->GetNodeValue(n_id, NodeIndex[k]);
					    else 
							val_n = m_pcs->GetNodeValue(n_id, NodeIndex[k]); //WW
						tec_file << val_n << " ";
						// WTP  - reactivated by JOD 2015-11-19 to match output with benchmark reference files - if SATURATION1 requested, SATURATION2 is added
						// WTP  - and deactivated again by wtp 19.11.2015
						//if ((m_pcs->type == 1212 || m_pcs->type == 42)
						//	&& _nod_value_vector[k].find("SATURATION") != string::npos) //WW
						//	tec_file << 1. - val_n << " ";
						////
					}
				}
			}
			//OK4704
			for (size_t k = 0; k < mfp_value_vector.size(); k++)
				//tec_file << MFPGetNodeValue(m_msh->nod_vector[j]->GetIndex(),mfp_value_vector[k]) << " "; //NB
				tec_file << MFPGetNodeValue(n_id,
				                            mfp_value_vector[k], atoi(&mfp_value_vector[k][mfp_value_vector[k].size() - 1]) - 1) << " ";  //NB: MFP output for all phases
		}
		tec_file << "\n";
	}
	_new_file_opened = true;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   08/2004 OK Implementation
   03/2005 OK MultiMSH
   08/2005 WW Wite for MultiMSH
   12/2005 OK GetMSH
   07/2007 NW Multi Mesh Type
**************************************************************************/
void COutput::WriteTECElementData(fstream &tec_file,int e_type)
{
	for (size_t i = 0; i < m_msh->ele_vector.size(); i++)
	{
		if (!m_msh->ele_vector[i]->GetMark())
			continue;
		//NW
		if (m_msh->ele_vector[i]->GetElementType() == e_type)
			m_msh->ele_vector[i]->WriteIndex_TEC(tec_file);
	}
}

/**************************************************************************
   FEMLib-Method:
   Programing:
   08/2004 OK Implementation
   08/2004 WW Header by the names gives in .out file
   03/2005 OK MultiMSH
   04/2005 WW Output active elements only
   08/2005 WW Output by MSH
   12/2005 OK GetMSH
**************************************************************************/
void COutput::WriteTECHeader(fstream &tec_file,int e_type, string e_type_name)
{
	// MSH
	//	m_msh = GetMSH();

	//OK411
	size_t no_elements = 0;
	const size_t mesh_ele_vector_size (m_msh->ele_vector.size());
	for (size_t i = 0; i < mesh_ele_vector_size; i++)
		if (m_msh->ele_vector[i]->GetMark())
			if (m_msh->ele_vector[i]->GetElementType() == e_type)
				no_elements++;
	//--------------------------------------------------------------------
	// Write Header I: variables
	CRFProcess* pcs = NULL;               //WW
	const size_t nName (_nod_value_vector.size());
	tec_file << "VARIABLES  = \"X\",\"Y\",\"Z\"";
	for (size_t k = 0; k < nName; k++)
	{
		tec_file << ", \"" << _nod_value_vector[k] << "\"";
		//-------------------------------------WW
		// WTP  - reactivated by JOD 2015-11-19 to match output with benchmark reference files - if SATURATION1 requested, SATURATION2 is added
		//pcs = GetPCS(_nod_value_vector[k]);
		//if (pcs != NULL)
		//	if ((pcs->type == 1212 ||
		//	     pcs->type == 42 ) && _nod_value_vector[k].find("SATURATION")
		//	    != string::npos)
		//		tec_file << ", SATURATION2";
		//-------------------------------------WW
	}
	const size_t mfp_value_vector_size (mfp_value_vector.size());
	for (size_t k = 0; k < mfp_value_vector_size; k++)
		//NB
		tec_file << ", \"" << mfp_value_vector[k] << "\"";

	// PCON
	//MX
	const size_t nPconName (_pcon_value_vector.size());
	for (size_t k = 0; k < nPconName; k++)
		//MX
		tec_file << ", " << _pcon_value_vector[k] << "";
	tec_file << "\n";

	//--------------------------------------------------------------------
	// Write Header II: zone
	tec_file << "ZONE T=\"";
	tec_file << _time << "s\", ";
	//OK411
	tec_file << "N=" << m_msh->GetNodesNumber(false) << ", ";
	tec_file << "E=" << no_elements << ", ";
	tec_file << "F=" << "FEPOINT" << ", ";
	tec_file << "ET=" << e_type_name;
	// JOD 2020-3-20 from BW - data accuracy for each variable
	tec_file << ", DT=(DOUBLE,DOUBLE,DOUBLE"; // BW, for the accuracy of the coordinates
	for (size_t k = 0; k < nName; k++) // BW, for the nodal variables, hard coded as SINGLE, i.e. 6 digits
	{
		tec_file << ",SINGLE";
	}
	tec_file << ")\n";
	//--------------------------------------------------------------------
    // Write Header III: solution time			; BG 05/2011
    tec_file << "STRANDID=1, SOLUTIONTIME=";
    tec_file << _time;			// << "s\"";
    tec_file << "\n";

    //--------------------------------------------------------------------
    // Write Header IV: Variable sharing		; BG 05/2011
    if (this->VARIABLESHARING == true)
    {
	   //int timestep = this->getNSteps;
      //if (this->
    }
	//
    if(_new_file_opened && tecplot_zone_share)  // 08.2012. WW
	{
        tec_file <<"VARSHARELIST=([1-3]=1)"<<"\n";
        tec_file <<"CONNECTIVITYSHAREZONE=1"<<"\n";
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   09/2004 OK Implementation
   01/2006 OK VAR,PCS,MSH concept
**************************************************************************/
void COutput::ELEWriteDOMDataTEC()
{
	//----------------------------------------------------------------------
	if(_ele_value_vector.empty ())
		return;
	//----------------------------------------------------------------------
	// File handling
	//......................................................................
	string tec_file_name = file_base_name + "_domain" + "_ele";
	if(getProcessType () != FiniteElement::INVALID_PROCESS) // PCS
		// 09/2010 TF msh_type_name;
		tec_file_name += "_" + convertProcessTypeToString (getProcessType());
	if(msh_type_name.size() > 1)          // MSH
		tec_file_name += "_" + msh_type_name;

#if defined(USE_PETSC)  //  8/2015 JOD 
	tec_file_name += "_" + mrank_str;
#endif

	tec_file_name += TEC_FILE_EXTENSION;
	//WW
	if(!_new_file_opened)
		remove(tec_file_name.c_str());
	//......................................................................
	fstream tec_file (tec_file_name.data(),ios::app | ios::out);
	tec_file.setf(ios::scientific,ios::floatfield);
	tec_file.precision(12);
	if (!tec_file.good())
		return;
	tec_file.seekg(0L,ios::beg);
#ifdef SUPERCOMPUTER
	// kg44 buffer the output
	char mybuffer [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
	tec_file.rdbuf()->pubsetbuf(mybuffer,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
	//
#endif

	//--------------------------------------------------------------------
	WriteELEValuesTECHeader(tec_file);
	WriteELEValuesTECData(tec_file);
	//--------------------------------------------------------------------
	_new_file_opened = true;
	tec_file.close();                     // kg44 close file
}

void COutput::WriteELEValuesTECHeader(fstream &tec_file)
{
	// Write Header I: variables
	/*tec_file << "VARIABLES = \"X\",\"Y\",\"Z\",\"VX\",\"VY\",\"VZ\",\"Rho\"";
	for (size_t i = 0; i < _ele_value_vector.size(); i++)
		//WW
		if (_ele_value_vector[i].find("VELOCITY") == string::npos || _ele_value_vector[i].find("DENSITY1") == string::npos)
			tec_file << "," << _ele_value_vector[i];
	*/
	tec_file << "VARIABLES = \"X\",\"Y\",\"Z\"";
	for (size_t i = 0; i < _ele_value_vector.size(); i++)
	{
		//WW
		if (_ele_value_vector[i].find("VELOCITY") != string::npos)
		{
			tec_file << ",\"VELOCITY1_X\",\"VELOCITY1_Y\",\"VELOCITY1_Z\"";
			break;
		}
	}

	for (size_t i = 0; i < _ele_value_vector.size(); i++)
		if (_ele_value_vector[i].find("DENSITY1") != string::npos)
		{
			tec_file << ",\"DENSITY1\"";
			break;
		}

	tec_file << "\n";
	// Write Header II: zone
	tec_file << "ZONE T=\"";
	tec_file << _time << "s\", ";
	tec_file << "I=" << (long) m_msh->ele_vector.size() << ", ";
	tec_file << "F=POINT" << ", ";
	tec_file << "C=BLACK";
	tec_file << "\n";
	//--------------------------------------------------------------------
	// Write Header III: solution time			// 8/2015 JOD
	//tec_file << "STRANDID=1, SOLUTIONTIME=";  
	//tec_file << _time;			// << "s\"";
	//tec_file << "\n";
  // removed by cb
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   09/2004 OK Implementation
   11/2005 OK MSH
   01/2006 OK
**************************************************************************/
void COutput::WriteELEValuesTECData(fstream &tec_file)
{
	CRFProcess* m_pcs = NULL;
	CRFProcess* m_pcs_2 = NULL;
	if (!_ele_value_vector.empty())
		m_pcs = GetPCS_ELE(_ele_value_vector[0]);
	else
		return;



	if (getProcessType() == FiniteElement::HEAT_TRANSPORT || getProcessType() == FiniteElement::MASS_TRANSPORT) // 8/2015 JOD
	{
		ele_gp_flux.clear();
		const size_t mesh_ele_vector_size(m_pcs->m_msh->ele_vector.size());
		for (size_t i = 0; i < mesh_ele_vector_size; i++)
			ele_gp_flux.push_back(new ElementValue(m_pcs, m_pcs->m_msh->ele_vector[i]));

		m_pcs = PCSGet(getProcessType());

		m_pcs->CalIntegrationPointValue();    //  calculate FICK / FOURRIER flux
		m_pcs->Extropolation_GaussValue();    //  and extrapolate to node
	}

	vector <bool> skip; // CB
	size_t no_ele_values = _ele_value_vector.size();
	bool out_element_vel = false;
	bool out_element_density = false; // JOD 2018-8-15

	for (size_t j = 0; j < no_ele_values; j++) //WW
	{
		if(_ele_value_vector[j].find("DENSITY") != string::npos)
		{
			out_element_density = true;
			skip.push_back(false);
		}
		else if(_ele_value_vector[j].find("VELOCITY") != string::npos)
		{
			out_element_vel = true;
			// break;  // CB: allow output of velocity AND other ele values
			skip.push_back(false);
		}
		else
		{
			m_pcs_2 = GetPCS_ELE(_ele_value_vector[j]);
			skip.push_back(true);
		}
   }
	vector<int> ele_value_index_vector(no_ele_values);
	GetELEValuesIndexVector(ele_value_index_vector);

	MeshLib::CElem* m_ele = NULL;
	FiniteElement::ElementValue* gp_ele = NULL;
	for (size_t i = 0; i < m_msh->ele_vector.size(); i++)
	{
		m_ele = m_msh->ele_vector[i];
		double const* xyz(m_ele->GetGravityCenter());
		tec_file << xyz[0] << " " << xyz[1] << " " << xyz[2] << " ";
		if (out_element_vel)      //WW
		{
			if (PCSGet(FiniteElement::FLUID_MOMENTUM)) // PCH 16.11 2009
			{
				CRFProcess* pch_pcs = PCSGet(FiniteElement::FLUID_MOMENTUM);

				tec_file << pch_pcs->GetElementValue(i,
				                                     pch_pcs->GetElementValueIndex(
				                                             "VELOCITY1_X") + 1)
				<< " ";
				tec_file << pch_pcs->GetElementValue(i,
				                                     pch_pcs->GetElementValueIndex(
				                                             "VELOCITY1_Y") + 1)
				<< " ";
				tec_file << pch_pcs->GetElementValue(i,
				                                     pch_pcs->GetElementValueIndex(
				                                             "VELOCITY1_Z") + 1)
				<< " ";
			}
			else
			{
				if (getProcessType() == FiniteElement::HEAT_TRANSPORT || getProcessType() == FiniteElement::MASS_TRANSPORT)
				  gp_ele = ele_gp_flux[i];
				else
				  gp_ele = ele_gp_value[i];

				tec_file << gp_ele->Velocity(0, 0) << " ";
				tec_file << gp_ele->Velocity(1, 0) << " ";
				tec_file << gp_ele->Velocity(2, 0) << " ";
			}
		}  // end out_element_vel

		if(out_element_density)
		{
			tec_file << gp_ele->density << " ";
		}
		
		for (size_t j = 0; j < ele_value_index_vector.size(); j++)
		{
			if (skip[j]) // CB: allow output of velocity AND other ele values
			{
				tec_file
					<< m_pcs_2->GetElementValue(i, ele_value_index_vector[j])
					<< " ";
			}
		}
		/*
		   int j;
		   int eidx;
		   char ele_value_char[20];
		   int no_ele_values = (int)ele_value_vector.size();
		   for(j=0;j<no_ele_values;j++){
		   sprintf(ele_value_char,"%s",ele_value_vector[j].data());
		   eidx = PCSGetELEValueIndex(ele_value_char);
		   tec_file << ElGetElementVal(i,eidx) << " ";
		   }
		 */
		tec_file << "\n";
	}

	if (ele_gp_flux.size() > 0)  // release memory 8/2015 JOD
	{
		for (int i = 0; i < (long)ele_gp_flux.size(); i++)
		{
			gp_ele = ele_gp_flux[i];
			delete gp_ele;
			gp_ele = NULL;
		}
		ele_gp_flux.clear();
	}

	ele_value_index_vector.clear();
    skip.clear();
}
/**************************************************************************
FEMLib-Method:
Task:
Programing:
10/2014 BW Inherited from NODWriteDOMDataTEC()
**************************************************************************/
void COutput::BLOCKWriteDOMDataTEC()
{
	int te = 0;
	string eleType;
	string tec_file_name;
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL)
	char tf_name[10];
	std::cout << "Process " << myrank << " in WriteDOMDataTEC" << "\n";
#endif
	//----------------------------------------------------------------------
	// Tests
	//OK4704
	if ((_nod_value_vector.size() == 0) && (mfp_value_vector.size() == 0)
		&& (_ele_value_vector.size() == 0))
		return;
	//......................................................................
	// MSH
	//m_msh = FEMGet(pcs_type_name);
	//  m_msh = GetMSH();
	if (!m_msh)
	{
		cout << "Warning in COutput::NODWriteDOMDataTEC() - no MSH data" << "\n";
		return;
	}


	// Output files for each mesh type
	//----------------------------------------------------------------------
	// File name handling
	tec_file_name = file_base_name + "_" + "domain" + "_block";
	if (getProcessType() != FiniteElement::INVALID_PROCESS) // PCS
		tec_file_name += "_" + convertProcessTypeToString(getProcessType());
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL)
	sprintf(tf_name, "%d", myrank);
	tec_file_name += "_" + string(tf_name);
	std::cout << "Tecplot filename: " << tec_file_name << "\n";
#endif
#if defined(USE_PETSC)  //|| defined(other parallel libs)//03.3012. WW
	tec_file_name += "_" + mrank_str;
	std::cout << "Tecplot filename: " << tec_file_name << "\n";
#endif
	tec_file_name += TEC_FILE_EXTENSION;
	//WW
	if (!_new_file_opened)
		remove(tec_file_name.c_str());
	fstream tec_file(tec_file_name.data(), ios::app | ios::out);
	tec_file.setf(ios::scientific, ios::floatfield);
	tec_file.precision(12);
	if (!tec_file.good())
		return;
#ifdef SUPERCOMPUTER
	// kg44 buffer the output
	char mybuf1[MY_IO_BUFSIZE * MY_IO_BUFSIZE];
	tec_file.rdbuf()->pubsetbuf(mybuf1, MY_IO_BUFSIZE * MY_IO_BUFSIZE);
#endif
	//
	WriteBLOCKValuesTECHeader(tec_file);
	WriteBLOCKValuesTECData(tec_file);

	// Output files for each mesh type
	//---------------------------------------
	if (tecplot_zone_share)
	{
		if (!_new_file_opened)
			WriteTECBLOCKData(tec_file);
	}
	else
	{
		WriteTECBLOCKData(tec_file);
	}

	_new_file_opened = true;
	tec_file.close();         // kg44 close file
}


/**************************************************************************
FEMLib-Method:
Task: Offer Connectity Information for Block Datapacking in TECPLOT
Programing:
10/2014 BW Implemented and inheritated from WriteTECElementData();
**************************************************************************/
void COutput::WriteTECBLOCKData(fstream &tec_file)
{
	for (size_t i = 0; i < m_msh->ele_vector.size(); i++)
	{
		if (!m_msh->ele_vector[i]->GetMark())
			continue;
		m_msh->ele_vector[i]->WriteIndex_TECBLOCK(tec_file);
	}
}
/**************************************************************************
FEMLib-Method:
Programing:
08/2004 OK Implementation
08/2004 WW Header by the names gives in .out file
03/2005 OK MultiMSH
04/2005 WW Output active elements only
08/2005 WW Output by MSH
12/2005 OK GetMSH
10/2014 BW Implemented for BLOCK DataPacking of TECPLOT
Inheritated from WriteTECHeader();
**************************************************************************/
void COutput::WriteBLOCKValuesTECHeader(fstream &tec_file) //BW: 23.03.2020 please update changes
{
	// MSH
	//	m_msh = GetMSH();

	//OK411
	size_t no_elements = 0;
	const size_t mesh_ele_vector_size(m_msh->ele_vector.size());
	no_elements = mesh_ele_vector_size;
	//--------------------------------------------------------------------
	// Write Header I: variables
	CRFProcess* pcs = NULL;               //WW
	const size_t nName(_nod_value_vector.size());
	tec_file << "VARIABLES = \"X\",\"Y\",\"Z\"";
	for (size_t k = 0; k < nName; k++)
	{
		tec_file << ", \"" << _nod_value_vector[k] << "\"";
		//-------------------------------------WW
		//pcs = GetPCS(_nod_value_vector[k]);                                     // WTP: This is a problem if we want to support 3 Phases
		//if (pcs != NULL)
		//	if ((pcs->type == 1212 ||
		//	     pcs->type == 42 ) && _nod_value_vector[k].find("SATURATION")
		//	    != string::npos)
		//		tec_file << ", SATURATION2";
		//-------------------------------------WW
	}
	const size_t mfp_value_vector_size(mfp_value_vector.size());
	for (size_t k = 0; k < mfp_value_vector_size; k++)
		//NB
		tec_file << ", \"" << mfp_value_vector[k] << "\"";

	// PCON
	//MX
	const size_t nPconName(_pcon_value_vector.size());
	for (size_t k = 0; k < nPconName; k++)
		//MX
		tec_file << ", " << _pcon_value_vector[k] << "";



	//MMP Value
	const size_t mmp_value_vector_size(mmp_value_vector.size());
	for (size_t k = 0; k < mmp_value_vector_size; k++)
		tec_file << ", \"" << mmp_value_vector[k] << "\"";


	//ELement value
	const size_t ele_value_vector_size(_ele_value_vector.size());
	for (size_t k = 0; k < ele_value_vector_size; k++)
		tec_file << ", \"" << _ele_value_vector[k] << "\"";
	tec_file << "\n";

	//--------------------------------------------------------------------
	// Write Header II: zone
	tec_file << "ZONE T=\"";
	tec_file << _time << "s\", ";
	tec_file << "Nodes=" << m_msh->GetNodesNumber(false) << ", ";
	tec_file << "Elements=" << no_elements << ", ";
	tec_file << "ZONETYPE=" << "FEBrick" << ", ";
	tec_file << "DATAPACKING=" << "BLOCK" << ", ";

	if (ele_value_vector_size > 0 || mmp_value_vector_size > 0)
	{
		if ((ele_value_vector_size == 1 && mmp_value_vector_size == 0))
		{
			tec_file << "VARLOCATION=" << "(["
				<< nName + mfp_value_vector_size + nPconName + 4
				<< "]=CELLCENTERED)";
		}

		else if ((ele_value_vector_size == 0 && mmp_value_vector_size == 1))
		{
			tec_file << "VARLOCATION=" << "(["
				<< nName + mfp_value_vector_size + nPconName + 4
				<< "]=CELLCENTERED)";
		}

		else
		tec_file << "VARLOCATION=" << "(["
			<< nName + mfp_value_vector_size + nPconName + 4
			<< "-"
			<< nName + mfp_value_vector_size + nPconName + ele_value_vector_size + mmp_value_vector_size + 3
			<< "]=CELLCENTERED)";
	}
	
	//data accuracy for each variable
	tec_file << "DT=(DOUBLE,DOUBLE,DOUBLE)"; // BW, for the accuracy of the coordinates

	tec_file << "\n";
	//--------------------------------------------------------------------
	// Write Header III: solution time			; BG 05/2011
	tec_file << "STRANDID=1, SOLUTIONTIME=";
	tec_file << _time;			// << "s\"";
	tec_file << "\n";

	if (_new_file_opened && tecplot_zone_share)  // 08.2012. WW
	{
		tec_file << "VARSHARELIST=([1-3]=1)" << "\n";
		tec_file << "CONNECTIVITYSHAREZONE=1" << "\n";
	}
}

/**************************************************************************
FEMLib-Method:
Programing:
08/2004 OK Implementation
08/2004 WW Output node variables by their names given in .out file
03/2005 OK MultiMSH
08/2005 WW Correction of node index
12/2005 OK Mass transport specifics
OK ??? too many specifics
10/2014 BW Write Block Datapacking Format of TECPLOT with Nodal/Ele Data
together, inheritated from WriteTECNodeData()
**************************************************************************/
void COutput::WriteBLOCKValuesTECData(fstream &tec_file) //BW: 23.03.2020 please update changes
{
	const size_t nName(_nod_value_vector.size());
	double val_n = 0.;                    //WW
	int nidx;
	vector<int> NodeIndex(nName);
	string nod_value_name;                //OK
	CNode *node = NULL;
	CRFProcess *deform_pcs = NULL; // 23.01.2012. WW. nulltpr

	int timelevel;
	//	m_msh = GetMSH();
	CRFProcess* m_pcs_out = NULL;
	// MSH
	for (size_t k = 0; k < nName; k++)
	{
		m_pcs = PCSGet(_nod_value_vector[k], true);
		if (m_pcs != NULL)
		{
			NodeIndex[k] = m_pcs->GetNodeValueIndex(_nod_value_vector[k], true); // JT Latest.
			if ((m_pcs->getProcessType() == FiniteElement::DEFORMATION)
				|| (m_pcs->getProcessType() == FiniteElement::DEFORMATION_DYNAMIC)
				|| (m_pcs->getProcessType() == FiniteElement::DEFORMATION_FLOW)
				|| (m_pcs->getProcessType() == FiniteElement::DEFORMATION_H2)
				)
			{
				deform_pcs = m_pcs;
			}
		}
	}

	// Coordinates
	bool out_coord = true;
	if (tecplot_zone_share && _new_file_opened)
		out_coord = false;

	if (out_coord) // List all coordinates
		for (size_t i = 0; i < 3; i++){
			for (size_t j = 0; j < m_msh->GetNodesNumber(false); j++){
				node = m_msh->nod_vector[j];
				// XYZ
				const double *x = node->getData();
				tec_file << x[i] << " ";
				if ((j + 1) % 5 == 0)
					tec_file << '\n';
			}
			tec_file << '\n';
		}


	// NOD values
	// Mass transport
	//     if(pcs_type_name.compare("MASS_TRANSPORT")==0){
	if (getProcessType() == FiniteElement::MASS_TRANSPORT)
		for (size_t i = 0; i < _nod_value_vector.size(); i++)
		{
			std::string nod_value_name = _nod_value_vector[i];
			for (size_t l = 0; l < pcs_vector.size(); l++)
			{
				m_pcs = pcs_vector[l];
				//					if (m_pcs->pcs_type_name.compare("MASS_TRANSPORT") == 0) {
				if (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
				{
					timelevel = 0;
					for (size_t m = 0;
						m < m_pcs->nod_val_name_vector.size(); m++)
						if (m_pcs->nod_val_name_vector[m].compare(
							nod_value_name) == 0)
						{
							m_pcs_out = PCSGet(FiniteElement::MASS_TRANSPORT,
								nod_value_name);
							if (!m_pcs_out)
								continue;
							if (timelevel == 1)
							{
								for (size_t j = 0; j < m_msh->GetNodesNumber(false); j++) // list all nodal vaule
								{
									node = m_msh->nod_vector[j];
									const size_t n_id = node->GetIndex();
									nidx =
										m_pcs_out->
										GetNodeValueIndex(
										nod_value_name)
										+
										timelevel;
									tec_file <<
										m_pcs_out->GetNodeValue(n_id, nidx) << " ";
									if ((j + 1) % 5 == 0)
										tec_file << '\n';
								}
							}
							timelevel++;
						}
				}
				tec_file << "\n";
			}
		}
	else
	{
		for (size_t k = 0; k < nName; k++)
		{
			m_pcs = GetPCS(_nod_value_vector[k]);
			if (m_pcs != NULL) { //WW
				for (size_t j = 0; j < m_msh->GetNodesNumber(false); j++)
				{
					node = m_msh->nod_vector[j];
					const size_t n_id = node->GetIndex();
					if (NodeIndex[k] > -1) {
						if ((m_pcs->type == 14) && (_nod_value_vector[k] == "HEAD"))
						{    // HEAD output for RICHARDS_FLOW (unconfined GW) 5.3.07 JOD 
							double rhow;
							double dens_arg[3];
							CRFProcess* m_pcs_transport = NULL;
							m_pcs_transport = PCSGet("MASS_TRANSPORT");

							if (m_pcs_transport) {
								dens_arg[2] = m_pcs_transport->GetNodeValue(n_id, 1);  // first component!!! 
								rhow = mfp_vector[0]->Density(dens_arg);   // first phase!!!  dens_arg
							}
							else
								rhow = mfp_vector[0]->Density();
							val_n = m_pcs->GetNodeValue(n_id, m_pcs->GetNodeValueIndex("PRESSURE1") + 1) / (rhow * 9.81) + m_msh->nod_vector[m_msh->nod_vector[j]->GetIndex()]->getData()[2];
							tec_file << val_n << " ";
							if ((j + 1) % 5 == 0)
								tec_file << '\n';
						}
						else {
							if (_nod_value_vector[k].find("DELTA") == 0) // BW 24/11/2015 
								val_n = m_pcs->GetNodeValue(n_id, 1) - m_pcs->GetNodeValue(n_id, NodeIndex[k]);
							else
								val_n = m_pcs->GetNodeValue(n_id, NodeIndex[k]); //WW
							tec_file << val_n << " ";
							//if ((m_pcs->type == 1212 || m_pcs->type == 42)
							//	&& _nod_value_vector[k].find("SATURATION") != string::npos) //WW
							//	tec_file << 1. - val_n << " ";
							if ((j + 1) % 5 == 0)
								tec_file << '\n';
						}
					}
				}
			}
			tec_file << '\n';
		}

		//OK4704
		for (size_t k = 0; k < mfp_value_vector.size(); k++) {
			for (size_t j = 0; j < m_msh->GetNodesNumber(false); j++)
			{
				node = m_msh->nod_vector[j];
				if ((j + 1) % 5 == 0)
					tec_file << '\n';
				const size_t n_id = node->GetIndex();
				//tec_file << MFPGetNodeValue(m_msh->nod_vector[j]->GetIndex(),mfp_value_vector[k]) << " "; //NB
				tec_file << MFPGetNodeValue(n_id,
					mfp_value_vector[k], atoi(&mfp_value_vector[k][mfp_value_vector[k].size() - 1]) - 1) << " ";  //NB: MFP output for all phases
			}
			tec_file << "\n";
		}
	}


	MeshLib::CElem* m_ele = NULL;
	FiniteElement::ElementValue* gp_ele = NULL;
	//MMP values
	size_t no_mmp_values = mmp_value_vector.size();
	CMediumProperties* m_mmp;
	double *tensor = NULL;
	int group;
	for (size_t j = 0; j < no_mmp_values; j++)
	{
		for (size_t i = 0; i < m_msh->ele_vector.size(); i++)
		{
			m_ele = m_msh->ele_vector[i];
			group = m_ele->GetPatchIndex();
			m_mmp = mmp_vector[group];
			if (mmp_value_vector[j].compare("PERMEABILITY") == 0)
			{
				tensor = m_mmp->PermeabilityTensor(i);
				tec_file << tensor[0] << " ";
			}
			if (mmp_value_vector[j].compare("PERMEABILITY_Y") == 0)
			{
				tensor = m_mmp->PermeabilityTensor(i);
				tec_file << tensor[4] << " ";
			}
			if (mmp_value_vector[j].compare("PERMEABILITY_Z") == 0)
			{
				tensor = m_mmp->PermeabilityTensor(i);
				tec_file << tensor[8] << " ";
			}
			if (mmp_value_vector[j].compare("POROSITY") == 0)
				tec_file << m_mmp->Porosity(i, 1) << " ";
			if (mmp_value_vector[j].compare("MATERIAL_GROUP") == 0)
				tec_file << group << " ";

			if ((i + 1) % 5 == 0)
				tec_file << '\n';
		}
	}


	//Element values
	//CRFProcess* m_pcs_2 = NULL;
	//if (!_ele_value_vector.empty())
	//	m_pcs = GetPCS_ELE(_ele_value_vector[0]);
	//else
	//	return;
	bool out_element_vel = false;
	if (_ele_value_vector.empty())
		return;
	else
		out_element_vel = true;
	vector <bool> skip; // CB
	size_t no_ele_values = _ele_value_vector.size();

	for (size_t j = 0; j < no_ele_values; j++) //WW
	{
		m_pcs = GetPCS_ELE(_ele_value_vector[j]);
		//
		//if (_ele_value_vector[j].find("VELOCITY") != string::npos)
		//{
		//	out_element_vel = true;
		//	//break;  // CB: allow output of velocity AND other ele values
		//	skip.push_back(false);
		//}
		//else
		//{
		//	m_pcs_2 = GetPCS_ELE(_ele_value_vector[j]);
		//	skip.push_back(true);
		//}
	}
	vector<int> ele_value_index_vector(no_ele_values);
	GetELEValuesIndexVector(ele_value_index_vector);




	if (out_element_vel)
	{
		//for (size_t j = 0; j < 3; j++){
		//	for (size_t i = 0; i < m_msh->ele_vector.size(); i++)
		//	{
		//		m_ele = m_msh->ele_vector[i];
		//		//if (PCSGet(FiniteElement::FLUID_MOMENTUM)){
		//		//	CRFProcess* pch_pcs = PCSGet(FiniteElement::FLUID_MOMENTUM);
		//		switch (j){
		//		case 0:
		//			tec_file << m_pcs->GetElementValue(i,
		//				m_pcs->GetElementValueIndex(
		//				"VELOCITY1_X") + 1)
		//				<< " ";
		//			break;
		//		case 1:
		//			tec_file << m_pcs->GetElementValue(i,
		//				m_pcs->GetElementValueIndex(
		//				"VELOCITY1_Y") + 1)
		//				<< " ";
		//			break;
		//		case 2:
		//			tec_file << m_pcs->GetElementValue(i,
		//				m_pcs->GetElementValueIndex(
		//				"VELOCITY1_Z") + 1)
		//				<< " ";
		//			break;
		//		}
		//		/*}
		//		else {
		//		gp_ele = ele_gp_value[i];
		//		tec_file << gp_ele->Velocity(j, 0) << " ";
		//		}*/
		//		if ((i + 1) % 5 == 0)
		//			tec_file << '\n';
		//	}
		//}
		//tec_file << '\n';


		for (size_t j = 0; j < no_ele_values; j++)
		{
			m_pcs = GetPCS_ELE(_ele_value_vector[j]);
			//if (skip[j]) // CB: allow output of velocity AND other ele values
			//{
			for (size_t i = 0; i < m_msh->ele_vector.size(); i++)
			{
				m_ele = m_msh->ele_vector[i];
				tec_file
					<< m_pcs->GetElementValue(i, ele_value_index_vector[j])
					<< " ";
				if ((i + 1) % 5 == 0)
					tec_file << '\n';
			}
			//}
			//}
		}
	}
	ele_value_index_vector.clear();
	skip.clear();
}
/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   08/2004 OK Implementation
   08/2004 WW Output by the order of distance to one end of the polyline
        OK this should be done in the PLY method
   08/2005 WW Changes due to the Geometry element class applied
           Remove existing files
   12/2005 OK Mass transport specifics
   12/2005 OK VAR,MSH,PCS concept
   12/2005 WW Output stress invariants
   08/2006 OK FLUX
   10/2008 OK MFP values
   07/2010 TF substituted GEOGetPLYByName
**************************************************************************/
double COutput::NODWritePLYDataTEC(int time_step_number)
{
	//WW  int nidx;
	long gnode;
	bool bdummy = false;
	int stress_i[6], strain_i[6];
	double ss[6];
	double val_n = 0.;                    //WW
	                                      //OK4704
	if ((_nod_value_vector.size() == 0) && (mfp_value_vector.size() == 0))
		return 0.0;

	// TF
	GEOLIB::Polyline const* const ply (
	        dynamic_cast<GEOLIB::Polyline const* const>(this->getGeoObj()));
	if (this->getGeoType() != GEOLIB::POLYLINE || ply == NULL)
	{
		std::cerr << "COutput::NODWritePLYDataTEC geometric object is not a polyline" <<
		"\n";
		return 0.0;
	}

	// File handling
	std::string tec_file_name = file_base_name + "_ply_" + geo_name + "_t"
	                            + number2str<size_t> (_id); //OK4709
	if (getProcessType() != FiniteElement::INVALID_PROCESS)
		tec_file_name += "_" + convertProcessTypeToString(getProcessType());
	if (msh_type_name.size() > 0)
		tec_file_name += "_" + msh_type_name;

#if defined(USE_PETSC)  // JOD 8/2015
	tec_file_name += "_" + mrank_str;
	std::cout << "Tecplot filename: " << tec_file_name << "\n";
#endif

	tec_file_name += TEC_FILE_EXTENSION;
	if (!_new_file_opened)
		remove(tec_file_name.c_str());
	else if (time_step_number == 0)  // JOD 2020-05-08
		return 0.;

	//WW
	fstream tec_file(tec_file_name.data(), ios::app | ios::out);

	tec_file.setf(ios::scientific, ios::floatfield);
	tec_file.precision(12);
	if (!tec_file.good())
		return 0.0;
	tec_file.seekg(0L, ios::beg);
#ifdef SUPERCOMPUTER
	// kg44 buffer the output
	char mybuffer [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
	tec_file.rdbuf()->pubsetbuf(mybuffer,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
	//
#endif
	//----------------------------------------------------------------------
	// Tests
	//......................................................................
	// GEO
//   CGLPolyline* m_ply = GEOGetPLYByName(geo_name);//CC
//   if (!m_ply)
//   {
//      cout << "Warning in COutput::NODWritePLYDataTEC - no GEO data" << "\n";
//      tec_file << "Warning in COutput::NODWritePLYDataTEC - no GEO data: "
//         << geo_name << "\n";
//      tec_file.close();
//      return 0.0;
//   }

	// MSH
	//	CFEMesh* m_msh = GetMSH();
	//	m_msh = GetMSH();
	if (!m_msh)
		cout << "Warning in COutput::NODWritePLYDataTEC - no MSH data" << "\n";
	//OKtec_file << "Warning in COutput::NODWritePLYDataTEC - no MSH data: " << geo_name << "\n";
	//OKtec_file.close();
	//OKToDo return;
	else
		m_msh->SwitchOnQuadraticNodes(false);  //WW

	// PCS
	if (getProcessType() == FiniteElement::INVALID_PROCESS)
		m_pcs = NULL;
	else
		m_pcs = PCSGet(getProcessType());

	CRFProcess* dm_pcs = NULL;            //WW
	for (size_t i = 0; i < pcs_vector.size(); i++)
		if (isDeformationProcess (pcs_vector[i]->getProcessType ()))
		{
			dm_pcs = pcs_vector[i];
			break;
		}

	/* //WW
	   // VEL
	   int v_eidx[3];
	   CRFProcess* m_pcs_flow (PCSGetFlow());
	   //m_pcs_flow = PCSGet("GROUNDWATER_FLOW"); //OKToDo
	   if (!m_pcs_flow)
	   {
	   //WW cout << "Warning in COutput::NODWritePLYDataTEC() - no PCS flow data" << "\n";
	   //tec_file << "Warning in COutput::NODWritePLYDataTEC() - no PCS flow data " << "\n";
	   //tec_file.close();
	   //return 0.0;
	   }
	   else
	   {
	   v_eidx[0] = m_pcs_flow->GetElementValueIndex("VELOCITY1_X");
	   v_eidx[1] = m_pcs_flow->GetElementValueIndex("VELOCITY1_Y");
	   v_eidx[2] = m_pcs_flow->GetElementValueIndex("VELOCITY1_Z");
	   }
	 */

//   for (size_t i = 0; i < 3; i++)
//   {
//      if (v_eidx[i] < 0)
//      {
//         //WW cout << "Warning in COutput::NODWritePLYDataTEC() - no PCS flow data" << "\n";
//         //tec_file << "Warning in COutput::NODWritePLYDataTEC() - no PCS flow data " << "\n";
//         //tec_file.close();
//      }
//   }
	//--------------------------------------------------------------------
	// NIDX for output variables
	size_t no_variables (_nod_value_vector.size());
	std::vector<int> NodeIndex(no_variables);
	GetNodeIndexVector(NodeIndex);
	//--------------------------------------------------------------------
	// Write header
	
  bool header = false; 

  if (aktueller_zeitschritt == 0)
    header = true;
  if (aktueller_zeitschritt > 0){
    if (aktueller_zeitschritt == nSteps)
      header = true;
    else if (time_vector.size()>0)
      if(fabs(aktuelle_zeit - time_vector[0]) < MKleinsteZahl) //WW MKleinsteZahl
        header = true;
  }
  
  if (header)       //WW if(number==1)  
  //if (number == 0 || number == 1)
 {
		//project_title;
		std::string project_title_string = "Profiles along polylines";
		tec_file << " TITLE = \"" << project_title_string << "\"" << "\n";
		tec_file << " VARIABLES = \"DIST\" ";

    for (size_t k = 0; k < no_variables; k++)
		{
			tec_file << "\"" << _nod_value_vector[k] << "\" ";
			//-------------------------------------WW
			// WTP
			//m_pcs = GetPCS(_nod_value_vector[k]);
			//if (m_pcs && m_pcs->type == 1212 && _nod_value_vector[k].find(
			//            "SATURATION") != string::npos)
			//	tec_file << "SATURATION2 ";
			//-------------------------------------WW
			if (_nod_value_vector[k].compare("FLUX") == 0)
				tec_file << "FLUX_INNER" << " ";
		}
		//....................................................................
		//OK4709
		for (size_t k = 0; k < mfp_value_vector.size(); k++)
			tec_file << "\"" << mfp_value_vector[k] << "\" ";
		//....................................................................
		// WW: M specific data
		if (dm_pcs)               //WW

			tec_file << " p_(1st_Invariant) " << " q_(2nd_Invariant)  "
			         << " Effective_Strain";
		tec_file << "\n";
	}
	//....................................................................
	// WW: M specific data
	size_t ns = 4;
	if (dm_pcs)                           //WW
	{
		stress_i[0] = dm_pcs->GetNodeValueIndex("STRESS_XX");
		stress_i[1] = dm_pcs->GetNodeValueIndex("STRESS_YY");
		stress_i[2] = dm_pcs->GetNodeValueIndex("STRESS_ZZ");
		stress_i[3] = dm_pcs->GetNodeValueIndex("STRESS_XY");
		strain_i[0] = dm_pcs->GetNodeValueIndex("STRAIN_XX");
		strain_i[1] = dm_pcs->GetNodeValueIndex("STRAIN_YY");
		strain_i[2] = dm_pcs->GetNodeValueIndex("STRAIN_ZZ");
		strain_i[3] = dm_pcs->GetNodeValueIndex("STRAIN_XY");
		if (max_dim == 2)         // 3D
		{
			ns = 6;
			stress_i[4] = dm_pcs->GetNodeValueIndex("STRESS_XZ");
			stress_i[5] = dm_pcs->GetNodeValueIndex("STRESS_YZ");
			strain_i[4] = dm_pcs->GetNodeValueIndex("STRAIN_XZ");
			strain_i[5] = dm_pcs->GetNodeValueIndex("STRAIN_YZ");
		}
	}
	//......................................................................
	// , I=" << NodeListLength << ", J=1, K=1, F=POINT" << "\n";
	tec_file << " ZONE T=\"TIME=" << _time << "\"";// << "\n";
	//----------------------------------------------------------------------
	//tec_file << " SOLUTIONTIME="; // JOD 8/2015 to link frames
	//tec_file << _time;			// << "s\"";
	tec_file << "\n";
	// Write data
	//======================================================================
	double flux_sum = 0.0;                //OK
	double flux_nod;

	m_msh->SwitchOnQuadraticNodes(false); //WW
	// NOD at PLY
	std::vector<long> nodes_vector;

	CGLPolyline* m_ply = GEOGetPLYByName(geo_name);
//   m_msh->GetNODOnPLY(m_ply, old_nodes_vector); // TF

	double tmp_min_edge_length (m_msh->getMinEdgeLength());
	m_msh->setMinEdgeLength (m_ply->epsilon);
	m_msh->GetNODOnPLY(ply, nodes_vector); // TF
	std::vector<double> interpolation_points;
	m_msh->getPointsForInterpolationAlongPolyline (ply, interpolation_points);
	m_msh->setMinEdgeLength(tmp_min_edge_length);

//   std::cout << "size of nodes_vector: " << nodes_vector.size() << ", size of old_nodes_vector: " << old_nodes_vector.size() << "\n";
	//bool b_specified_pcs = (m_pcs != NULL); //NW m_pcs = PCSGet(pcs_type_name);
	for (size_t j(0); j < nodes_vector.size(); j++)
	{
//		tec_file << m_ply->getSBuffer()[j] << " ";
		tec_file << interpolation_points[j] << " ";
		//WW
//		long old_gnode = nodes_vector[m_ply->getOrderedPoints()[j]];
		gnode = nodes_vector[j];
		for (size_t k = 0; k < no_variables; k++)
		{
			//if(!(_nod_value_vector[k].compare("FLUX")==0))  // removed JOD, does not work for multiple flow processes
			//if (!b_specified_pcs) //NW
			if (msh_type_name != "COMPARTMENT") // JOD 4.10.01
				m_pcs = PCSGet(_nod_value_vector[k], bdummy); // BW, here define which process for what secondary variable

			if (!m_pcs)
			{
				cout << "Warning in COutput::NODWritePLYDataTEC - no PCS data"
				     << "\n";
				tec_file
				<< "Warning in COutput::NODWritePLYDataTEC - no PCS data"
				<< "\n";
				return 0.0;
			}
			// WW
//			double old_val_n = m_pcs->GetNodeValue(old_gnode, NodeIndex[k]);
			if (_nod_value_vector[k].find("DELTA") == 0) // JOD 2014-11-10
				val_n = m_pcs->GetNodeValue(gnode, 1) - m_pcs->GetNodeValue(gnode, NodeIndex[k]);
			else
			    val_n = m_pcs->GetNodeValue(gnode, NodeIndex[k]);
//			tec_file << old_val_n << " ";
			tec_file << val_n << " ";
			// WTP
			//if (m_pcs->type == 1212 && (_nod_value_vector[k].find("SATURATION")
			//                            != string::npos))
			//	tec_file << 1. - val_n << " ";

			if (_nod_value_vector[k].compare("FLUX") == 0)
			{
				if (aktueller_zeitschritt == 0) //OK
					flux_nod = 0.0;
				else
					flux_nod = NODFlux(gnode);
				tec_file << flux_nod << " ";
				//flux_sum += abs(m_pcs->eqs->b[gnode]);
				flux_sum += abs(flux_nod);
				//OK cout << gnode << " " << flux_nod << " " << flux_sum << "\n";
			}
		}
		if (dm_pcs) //WW
		{
			for (size_t i = 0; i < ns; i++)
				ss[i] = dm_pcs->GetNodeValue(gnode, stress_i[i]);
			tec_file << -DeviatoricStress(ss) / 3.0 << " ";
			tec_file << sqrt(3.0 * TensorMutiplication2(ss, ss,
			                                            m_msh->GetCoordinateFlag() /
			                                            10) / 2.0) << "  ";
			for (size_t i = 0; i < ns; i++)
				ss[i] = dm_pcs->GetNodeValue(gnode, strain_i[i]);
			DeviatoricStress(ss);
			tec_file << sqrt(3.0 * TensorMutiplication2(ss, ss,
			                                            m_msh->GetCoordinateFlag() /
			                                            10) / 2.0);
		}

		// MFP //OK4704
		//OK4704
		for (size_t k = 0; k < mfp_value_vector.size(); k++)
			//     tec_file << MFPGetNodeValue(gnode,mfp_value_vector[k],0) << " "; //NB
			tec_file << MFPGetNodeValue(gnode, mfp_value_vector[k], atoi(
			                                    &mfp_value_vector[k][mfp_value_vector[k
			                                                         ].size() - 1]) - 1)
			<< " ";  //NB: MFP output for all phases

		tec_file << "\n";
	}
	_new_file_opened = true;
	tec_file.close();                     // kg44 close file
	return flux_sum;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   08/2004 OK Implementation
   08/2005 WW MultiMesh
   12/2005 OK VAR,MSH,PCS concept
   12/2005 WW Output stress invariants
   10/2010 TF changed access to process type
**************************************************************************/
void COutput::NODWritePNTDataTEC(int time_step_number)
{

//#if defined(USE_PETSC)  // JOD 2015-11-17
//	std::cout << "point_" + mrank_str;
//#endif

	long msh_node_number(m_msh->GetNODOnPNT(
	                             static_cast<const GEOLIB::Point*> (getGeoObj())));
        if(msh_node_number < 0)  //11.06.2012. WW
	  return;

	CRFProcess* dm_pcs = NULL;
	for (size_t i = 0; i < pcs_vector.size(); i++)
		//		if (pcs_vector[i]->pcs_type_name.find("DEFORMATION") != string::npos) { TF
		if (isDeformationProcess (pcs_vector[i]->getProcessType()))
		{
			dm_pcs = pcs_vector[i];
			break;
		}

	// File handling
	std::string tec_file_name(file_base_name + "_time_");
	addInfoToFileName(tec_file_name, true, true, true);

	if (!_new_file_opened)
		remove(tec_file_name.c_str());
	else if (time_step_number == 0)  // JOD 2020-05-08
		return;
	//......................................................................
	fstream tec_file(tec_file_name.data(), ios::app | ios::out);

	tec_file.setf(ios::scientific, ios::floatfield);
	tec_file.precision(12);
	if (!tec_file.good())
		return;
	tec_file.seekg(0L, ios::beg);
#ifdef SUPERCOMPUTER
	// kg44 buffer the output
	char mybuffer [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
	tec_file.rdbuf()->pubsetbuf(mybuffer,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
	//
#endif
	//--------------------------------------------------------------------
	// Tests
	//......................................................................
	//	CFEMesh* m_msh = GetMSH();
	//	m_msh = GetMSH();
	if (!m_msh)
	{
		cout << "Warning in COutput::NODWritePNTDataTEC - no MSH data: "
		     << "\n";
		tec_file << "Warning in COutput::NODWritePNTDataTEC - no MSH data: "
		         << "\n";
		tec_file.close();
		return;
	}

	//----------------------------------------------------------------------
	// NIDX for output variables
	size_t no_variables(_nod_value_vector.size());
	vector<int> NodeIndex(no_variables);
	GetNodeIndexVector(NodeIndex);
	//--------------------------------------------------------------------
	// Write header
	if (time_step_number == 0)            //WW  Old: if(time_step_number==1)
	{
		//project_title;
		
			const std::string project_title_string ("Time curves in points");
			tec_file << " TITLE = \"" << project_title_string << "\"" << "\n";
		

		tec_file << " VARIABLES = \"TIME \" ";

		//    if(pcs_type_name.compare("RANDOM_WALK")==0)
		if (getProcessType() == FiniteElement::RANDOM_WALK)
			tec_file << "leavingParticles ";
		for (size_t k = 0; k < no_variables; k++) //WW
		{
			tec_file << " \"" << _nod_value_vector[k] << "\" ";
			//-------------------------------------WW
			// WTP
			//m_pcs = GetPCS(_nod_value_vector[k]);
			//if (m_pcs && m_pcs->type == 1212 && _nod_value_vector[k].find(
			//            "SATURATION") != string::npos)
			//	tec_file << "SATURATION2 ";
			//-------------------------------------WW
		}
		//OK411
		for (size_t k = 0; k < mfp_value_vector.size(); k++)
			//NB MFP data names for multiple phases
			tec_file << " \"" << mfp_value_vector[k] << "\" ";
		//
#ifdef RFW_FRACTURE
		for(i = 0; i < (int)mmp_vector.size(); ++i)
			if( mmp_vector[i]->frac_num > 0)
				for(int j = 0; j < mmp_vector[i]->frac_num; ++j)
					tec_file << mmp_vector[i]->frac_names[j] << "_k " <<
					mmp_vector[i]->frac_names[j] << "_aper "
					         << mmp_vector[i]->frac_names[j] << "_closed ";

#endif

		if (dm_pcs)               //WW
			tec_file << " p_(1st_Invariant) " << " q_(2nd_Invariant)  "
			         << " Effective_Strain";
		tec_file << "\n";

		if (dat_type_name.compare("GNUPLOT") == 0) // 5.3.07 JOD
		  tec_file << "# "; // comment

		if (geo_name.find("POINT") == std::string::npos)
			tec_file << " ZONE T=\"POINT="
			//, I=" << anz_zeitschritte << ", J=1, K=1, F=POINT" << "\n";
			         << getGeoTypeAsString() << geo_name << "\"" << "\n";
		else
			tec_file << " ZONE T=\"POINT=" << geo_name << "\""
			         << "\n";  //, I=" << anz_zeitschritte << ", J=1, K=1, F=POINT" << "\n";
	}

	// For deformation
	size_t ns = 4;
	int stress_i[6], strain_i[6];
	double ss[6];
	if (dm_pcs)                           //WW
	{
		stress_i[0] = dm_pcs->GetNodeValueIndex("STRESS_XX");
		stress_i[1] = dm_pcs->GetNodeValueIndex("STRESS_YY");
		stress_i[2] = dm_pcs->GetNodeValueIndex("STRESS_ZZ");
		stress_i[3] = dm_pcs->GetNodeValueIndex("STRESS_XY");
		strain_i[0] = dm_pcs->GetNodeValueIndex("STRAIN_XX");
		strain_i[1] = dm_pcs->GetNodeValueIndex("STRAIN_YY");
		strain_i[2] = dm_pcs->GetNodeValueIndex("STRAIN_ZZ");
		strain_i[3] = dm_pcs->GetNodeValueIndex("STRAIN_XY");
		if (m_msh->GetCoordinateFlag() / 10 == 3) // 3D
		{
			ns = 6;
			stress_i[4] = dm_pcs->GetNodeValueIndex("STRESS_XZ");
			stress_i[5] = dm_pcs->GetNodeValueIndex("STRESS_YZ");
			strain_i[4] = dm_pcs->GetNodeValueIndex("STRAIN_XZ");
			strain_i[5] = dm_pcs->GetNodeValueIndex("STRAIN_YZ");
		}
	}
	//--------------------------------------------------------------------
	// Write data
	//......................................................................
	tec_file << aktuelle_zeit << " ";
	//......................................................................
	// NOD values
	if (getProcessType() == FiniteElement::RANDOM_WALK)
		tec_file << m_msh->PT->leavingParticles << " ";
	int timelevel;
	CRFProcess* m_pcs_out = NULL;

	// fetch geometric entities, especial the associated GEOLIB::Point vector
	if (pcs_vector[0] == NULL)
		return;

	//11.06.2012. WW// long msh_node_number(m_msh->GetNODOnPNT(
	//                             static_cast<const GEOLIB::Point*> (getGeoObj())));

	// Mass transport
	if (getProcessType() == FiniteElement::MASS_TRANSPORT)
		for (size_t i = 0; i < _nod_value_vector.size(); i++)
		{
			std::string nod_value_name = _nod_value_vector[i];
			for (size_t l = 0; l < pcs_vector.size(); l++)
			{
				m_pcs = pcs_vector[l];
				//				if (m_pcs->pcs_type_name.compare("MASS_TRANSPORT") == 0) { TF
				if (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
				{
					timelevel = 0;
					for (size_t m = 0; m < m_pcs->nod_val_name_vector.size();
					     m++)
						if (m_pcs->nod_val_name_vector[m].compare(
						            nod_value_name) == 0)
						{
							//							m_pcs_out = PCSGet(pcs_type_name, nod_value_name);
							m_pcs_out
							        = PCSGet(FiniteElement::MASS_TRANSPORT,
							                 nod_value_name);
							if (timelevel == 1)
							{
								int nidx =
								        m_pcs_out->
								        GetNodeValueIndex(
								                nod_value_name) +
								        timelevel;
								tec_file << m_pcs_out->GetNodeValue(
								        msh_node_number,
								        nidx) << " ";
							}
							timelevel++;
						}
				}
			}
		}
	else
	{
		double flux_nod, flux_sum = 0.0;
		for (size_t i = 0; i < _nod_value_vector.size(); i++)
		{
			// PCS
			if (!(_nod_value_vector[i].compare("FLUX") == 0)
				|| getProcessType() == FiniteElement::OVERLAND_FLOW) //JOD separate infiltration flux output in overland flow

				m_pcs = GetPCS(_nod_value_vector[i]);
			else
				m_pcs = GetPCS();
			if (!m_pcs)
			{
				cout << "Warning in COutput::NODWritePLYDataTEC - no PCS data"
				     << "\n";
				tec_file
				<< "Warning in COutput::NODWritePLYDataTEC - no PCS data"
				<< "\n";
				return;
			}
			//..................................................................
			// PCS
			if (!(_nod_value_vector[i].compare("FLUX") == 0)
				|| getProcessType() == FiniteElement::OVERLAND_FLOW) // JOD separate infiltration flux output in overland flow
			{
				//-----------------------------------------WW
				double val_n;
				
				if (_nod_value_vector[i].find("DELTA") == 0) //JOD 2014-11-10
					val_n = m_pcs->GetNodeValue(msh_node_number, 1) - m_pcs->GetNodeValue(msh_node_number, NodeIndex[i]);
				else
				    val_n = m_pcs->GetNodeValue(msh_node_number, NodeIndex[i]);
				tec_file << val_n << " ";
				// WTP
				//m_pcs = GetPCS(_nod_value_vector[i]);
				//if (m_pcs->type == 1212 && (_nod_value_vector[i].find(
				//                                    "SATURATION") != string::npos))
				//	tec_file << 1. - val_n << " ";
				//-----------------------------------------WW
			}
			else
			{
				flux_nod = NODFlux(msh_node_number);
				tec_file << flux_nod << " ";
				//flux_sum += abs(m_pcs->eqs->b[gnode]);
				flux_sum += abs(flux_nod);
				//OK cout << gnode << " " << flux_nod << " " << flux_sum << "\n";
			}
		}
		//....................................................................
#ifdef RFW_FRACTURE
		for(i = 0; i < (int)mmp_vector.size(); ++i)
			if( mmp_vector[i]->frac_num > 0)
				for(int j = 0; j < mmp_vector[i]->frac_num; ++j)
					tec_file << mmp_vector[i]->frac_perm[j] << " " <<
					mmp_vector[i]->avg_aperture[j] << " "
					         << mmp_vector[i]->closed_fraction[j] << " ";

#endif
		//....................................................................
		if (dm_pcs)               //WW
		{
			for (size_t i = 0; i < ns; i++)
				ss[i] = dm_pcs->GetNodeValue(msh_node_number, stress_i[i]);
			tec_file << -DeviatoricStress(ss) / 3.0 << " ";
			tec_file << sqrt(3.0 * TensorMutiplication2(ss, ss,
			                                            m_msh->GetCoordinateFlag() /
			                                            10) / 2.0) << "  ";
			for (size_t i = 0; i < ns; i++)
				ss[i] = dm_pcs->GetNodeValue(msh_node_number, strain_i[i]);
			DeviatoricStress(ss);
			tec_file << sqrt(3.0 * TensorMutiplication2(ss, ss,
			                                            m_msh->GetCoordinateFlag() /
			                                            10) / 2.0);
		}
		//OK411
		for (size_t k = 0; k < mfp_value_vector.size(); k++)
			tec_file << MFPGetNodeValue(msh_node_number, mfp_value_vector[k],
			                            atoi(&mfp_value_vector[k][mfp_value_vector[k].
			                                                      size() - 1])
			                            - 1) << " ";  //NB
	}
	tec_file << "\n";
	_new_file_opened = true;
	//----------------------------------------------------------------------
	tec_file.close();                     // kg44 close file
}

void COutput::WriteRFOHeader(fstream &rfo_file)
{
	//#0#0#0#1#0.00000000000000e+000#0#3915###########################################
	rfo_file << "#0#0#0#1#";
	rfo_file << _time;
	rfo_file << "#0#";
	rfo_file << OGS_VERSION;
	rfo_file << "###########################################";
	rfo_file << "\n";
}

void COutput::WriteRFONodes(fstream &rfo_file)
{
	//0 101 100
	rfo_file << 0 << " " << m_msh->nod_vector.size() << " "
	         << m_msh->ele_vector.size() << "\n";
	//0 0.00000000000000 0.00000000000000 0.00000000000000
	for (size_t i = 0; i < m_msh->nod_vector.size(); i++)
	{
		double const* const pnt_i (m_msh->nod_vector[i]->getData());
		rfo_file << i << " " << pnt_i[0] << " " << pnt_i[1] << " "
		         << pnt_i[2] << " " << "\n";
	}
}

void COutput::WriteRFOElements(fstream &rfo_file)
{
	size_t j;
	MeshLib::CElem* m_ele = NULL;
	//0 0 -1 line 0 1
	for(long i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		m_ele = m_msh->ele_vector[i];
		rfo_file << i << " " << m_ele->GetPatchIndex() \
		         << " -1 " \
		         << m_ele->GetName() << " ";
		for(j = 0; j < m_ele->GetNodesNumber(false); j++)
			rfo_file << m_ele->getNodeIndices()[j] << " ";
		rfo_file << "\n";
	}
}

void COutput::WriteRFOValues(fstream &rfo_file)
{
	int p,nidx;
	CRFProcess* m_pcs = NULL;
	//1 2 4
	rfo_file << "1 1 4" << "\n";
	//2 1 1
	int no_processes = (int)pcs_vector.size();
	rfo_file << no_processes;
	for(p = 0; p < no_processes; p++)
		rfo_file << " 1";
	rfo_file << "\n";
	//PRESSURE1, Pa
	// Names and units
	for(p = 0; p < no_processes; p++)
	{
		m_pcs = pcs_vector[p];
		rfo_file << m_pcs->pcs_primary_function_name[0];
		rfo_file << ", ";
		rfo_file << m_pcs->pcs_primary_function_unit[0];
		rfo_file << "\n";
	}
	//0 0.00000000000000 0.00000000000000
	// Node values
	for(long i = 0; i < (long)m_msh->nod_vector.size(); i++)
	{
		rfo_file << i;
		for(p = 0; p < no_processes; p++)
		{
			m_pcs = pcs_vector[p];
			nidx = m_pcs->GetNodeValueIndex(m_pcs->pcs_primary_function_name[0]) + 1;
			rfo_file << " " << m_pcs->GetNodeValue(i,nidx);
		}
		rfo_file << " " << "\n";
	}
}

void COutput::WriteRFO()
{
	m_msh = FEMGet(convertProcessTypeToString (getProcessType()));
	if(!m_msh)
	{
		cout << "Warning in COutput::WriteRFONodes - no MSH data" << "\n";
		return;
	}
	//--------------------------------------------------------------------
	// File handling
	string rfo_file_name;
	rfo_file_name = file_base_name + "." + "rfo";
	//WW
	if(!_new_file_opened)
		remove(rfo_file_name.c_str());
	fstream rfo_file (rfo_file_name.data(),ios::app | ios::out);

	rfo_file.setf(ios::scientific,ios::floatfield);
	rfo_file.precision(12);
	if (!rfo_file.good())
		return;
	rfo_file.seekg(0L,ios::beg);
#ifdef SUPERCOMPUTER
	// kg44 buffer the output
	char mybuffer [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
	rfo_file.rdbuf()->pubsetbuf(mybuffer,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
	//
#endif
	//--------------------------------------------------------------------
	WriteRFOHeader(rfo_file);
	WriteRFONodes(rfo_file);
	WriteRFOElements(rfo_file);
	WriteRFOValues(rfo_file);
	//  RFOWriteELEValues();
	rfo_file.close();                     // kg44 close file
}

void COutput::NODWriteSFCDataTEC(int time_step_number)
{

	/*   CB:   Extended for 2D-Element projection along a regular surface   */
	if (_nod_value_vector.size() == 0)
	{
		std::cout << "Warning - No nodes for surface " << geo_name << "\n";
		return;
	}

	m_msh = FEMGet(convertProcessTypeToString(getProcessType()));
	m_pcs = PCSGet(getProcessType());

	// File handling
	char number_char[6];
	sprintf(number_char, "%i", time_step_number);
	string number_string(number_char);

	//	string tec_file_name = pcs_type_name + "_sfc_" + geo_name + "_t"
	//				+ number_string + TEC_FILE_EXTENSION;
   //std::string tec_file_name = convertProcessTypeToString (getProcessType()) + "_sfc_" + geo_name + "_t"
   //   + number_string + TEC_FILE_EXTENSION;
   // AB SB Use Model name for output file name
   // std::string tec_file_name = convertProcessTypeToString (getProcessType()) 
   //std::string tec_file_name = file_base_name 
		 //                       + "_sfc_" + geo_name + "_t"
	  //                          + number_string + TEC_FILE_EXTENSION;   
   std::string tec_file_name = file_base_name
     + "_sfc_" + geo_name + TEC_FILE_EXTENSION;

   std::fstream tec_file;
	if (!_new_file_opened)
		remove(tec_file_name.c_str());
	else if (time_step_number == 0)  // JOD 2020-05-08
		return;

	tec_file.open(tec_file_name.data(), ios::app);

	tec_file.setf(ios::scientific, ios::floatfield);
	tec_file.precision(12);
	if (!tec_file.good())
	{
		std::cout << "Warning - Could not open file for writing surface data " << geo_name << "\n";
		return;
	}
	tec_file.seekg(0L, ios::beg);
#ifdef SUPERCOMPUTER
	// kg44 buffer the output
	char mybuffer [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
	tec_file.rdbuf()->pubsetbuf(mybuffer,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
	//
#endif
	//--------------------------------------------------------------------
	// Write header
	//project_title;
	string project_title_string = "Profile at surface";
	tec_file << " TITLE = \"" << project_title_string << "\""
	         << "\n";
	tec_file << " VARIABLES = \"X\",\"Y\",\"Z\",";
	for (size_t k = 0; k < _nod_value_vector.size(); k++)
		tec_file << _nod_value_vector[k] << ",";
	tec_file << "\n";
	// , I=" << NodeListLength << ", J=1, K=1, F=POINT" << "\n";
    tec_file << " ZONE T=\"TIME=" << _time << "\""; // << "\n";
	//--------------------------------------------------------------------
	// Write data
	std::vector<long> nodes_vector;
	Surface* m_sfc = NULL;
	m_sfc = GEOGetSFCByName(geo_name);    //CC
	if (m_sfc)
	{
		m_msh->GetNODOnSFC(m_sfc, nodes_vector);
		//for (size_t i = 0; i < m_msh->nod_vector.size(); i++)

		// CB set up data for Element section
		//std::vector < std::vector <long> > eledefvec;
		std::vector <long> eledef;
		bool* elecheck;
		if (aktueller_zeitschritt == 0)
		{

		  MeshLib::CElem* m_ele = NULL;
		  MeshLib::CNode* m_node = NULL;

		  elecheck = new bool[m_msh->ele_vector.size()];
		  for (size_t i = 0; i < m_msh->ele_vector.size(); i++)
			elecheck[i] = false;


		  for (size_t i = 0; i < nodes_vector.size(); i++) // AB SB
		  {
			int nd = nodes_vector[i];
			m_node = m_msh->nod_vector[nd];

			//test all connected elements of a node
			for (size_t j = 0; j < m_node->getNumConnectedElements(); j++) // AB SB
			{
			  m_ele = m_msh->ele_vector[m_node->getConnectedElementIDs()[j]];
			  if (elecheck[m_ele->GetIndex()] == true)
				continue;
			  //test all faces of the element
			  for (size_t k = 0; k < m_ele->GetFacesNumber(); k++) // AB SB
			  {

				int faceNodeIndex_loc[10];
				int faceNodeIndex_glo[10];
				int faceNodeIndex_lis[10];
				// get # face nodes and their local indices
				int n_face_node = m_ele->GetElementFaceNodes(k, faceNodeIndex_loc);
				int n = 0; // reset counter

				// check if all nodes of a face are in the list of surface nodes
				for (size_t l = 0; l < n_face_node; l++) // AB SB
				{
				  //get global Node index for comparison
				  faceNodeIndex_glo[l] = m_ele->GetNodeIndex(faceNodeIndex_loc[l]);
				  // compare with all nodes of surface node list
				  for (size_t m = 0; m < nodes_vector.size(); m++) // AB SB
				  {
					if (faceNodeIndex_glo[l] == nodes_vector[m])
					{
					  n++; // a match --> save the list index
					  faceNodeIndex_lis[l] = m;
					}
					// check if all nodes of face match
					if (n == n_face_node)
					  break;  // face identified, skip rest of loop now
				  }  // end nodes_vector
				  // check if face in surface
				  if (n == n_face_node)
				  {
					for (size_t l = 0; l < n_face_node; l++)
					{
					  eledef.push_back(1 + faceNodeIndex_lis[l]);
					}
					eledefvec.push_back(eledef);
					eledef.clear();
				  }  // end n == n_face_node
				}  // end n_face_node

			  }  // end faces
			  elecheck[m_ele->GetIndex()] = true;
			}  // end elements
		  }  // end nodes_vector

		  delete[] elecheck;
		  // Element section finished
		}  // end aktueller_zeitschritt == 0

		// CB extend header
		tec_file << ", N = " << nodes_vector.size() << ", E = " << eledefvec.size() << ", F = FEPOINT, ET = ";
		if (eledefvec[0].size() == 3)
		  tec_file << "TRIANGLE" << "\n";
		if (eledefvec[0].size() == 4)
		  tec_file << "QUADRILATERAL" << "\n";

		// node values
		for (size_t i = 0; i < nodes_vector.size(); i++) // AB SB
		{
			double const* const pnt_i (m_msh->nod_vector[nodes_vector[i]]->getData());
			tec_file << pnt_i[0] << " ";
			tec_file << pnt_i[1] << " ";
			tec_file << pnt_i[2] << " ";
			for (size_t k = 0; k < _nod_value_vector.size(); k++)
			{
				m_pcs = PCSGet(_nod_value_vector[k], true); // AB SB
				int nidx = m_pcs->GetNodeValueIndex(_nod_value_vector[k]) + 1;

				if (_nod_value_vector[k].find("DELTA") == 0) // JOD 2014-11-10
					tec_file << m_pcs->GetNodeValue(nodes_vector[i], 1) - m_pcs->GetNodeValue(nodes_vector[i], nidx - 1) << " ";
				else
				tec_file << m_pcs->GetNodeValue(nodes_vector[i], nidx) << " ";
			}
			tec_file << "\n";
		}  // end nodes_vector

		// CB print element section
		for (size_t i = 0; i < eledefvec.size(); i++)
		{
		  for (size_t j = 0; j < eledefvec[j].size(); j++)
			tec_file << eledefvec[i][j] << " ";
		  tec_file << "\n";
		}
		//clean up
		//eledefvec.clear();
		eledef.clear();

	}  // end m_sfc
	else
		tec_file << "Error in NODWriteSFCDataTEC: Surface " << geo_name
		         << " not found" << "\n";

	_new_file_opened = true;
	tec_file.close();                     // kg44 close file
}


/**************************************************************************
   FEMLib-Method:
   12/2005 OK Implementation
   04/2006 CMCD no mesh option & flux weighting
   last modification:
   03/2018 JOD rewritten to use FaceIntegration
   	   	   works only for primary variables since index is incremented
**************************************************************************/
void COutput::NODWriteSFCAverageDataTEC(int time_step_number)
{

	/*   CB:   Extended for 2D-Element projection along a regular surface   */
		if (_nod_value_vector.size() == 0)
		{
			std::cout << "Warning - No nodes for surface " << geo_name << "\n";
			return;
		}

		m_msh = FEMGet(convertProcessTypeToString(getProcessType()));
		m_pcs = PCSGet(getProcessType());
		if(m_pcs == NULL)
		{
			std::cout << "Warning - PCS not known for surface-averaged output" << "\n";
			return;
		}

		// File handling
		int number=1;
		char number_char[6];
		sprintf(number_char, "%i", number);
		string number_string(number_char);


		std::string tec_file_name = file_base_name
	     + "_sfc_" + geo_name + "_" + std::string(convertProcessTypeToString(getProcessType()))
	    		 + "_averaged" + TEC_FILE_EXTENSION;

	   std::fstream tec_file;
		if (!_new_file_opened)
			remove(tec_file_name.c_str());
		else if (time_step_number == 0)  // JOD 2020-05-08
			return;

	    tec_file.open(tec_file_name.data(), ios::app|ios::out);

		tec_file.setf(ios::scientific, ios::floatfield);
		tec_file.precision(12);
		if (!tec_file.good())
		{
			std::cout << "Warning - Could not open file for writing surface data " << geo_name << "\n";
			return;
		}
		tec_file.seekg(0L, ios::beg);
	#ifdef SUPERCOMPUTER
		// kg44 buffer the output
		char mybuffer [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
		tec_file.rdbuf()->pubsetbuf(mybuffer,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
		//
	#endif
		//--------------------------------------------------------------------
		// Write header
		if (aktueller_zeitschritt == 0)
		{
			tec_file << "Time "; // << "\n";
			for(int n = 0; n < _nod_value_vector.size(); n++)
				tec_file << _nod_value_vector[n] << " ";
			tec_file << "\n";
		}

		//--------------------------------------------------------------------
		// Write data
		std::vector<long> nodes_vector;
		Surface* m_sfc = NULL;
		m_sfc = GEOGetSFCByName(geo_name);    //CC

		double total_area, average_value;
		m_msh->GetNODOnSFC(m_sfc, nodes_vector);
		std::vector<double> nodes_area_vector(nodes_vector.size(), 1.);

		FaceIntegration(m_msh, nodes_vector, nodes_area_vector, m_sfc, FiniteElement::CONSTANT_NEUMANN,
				m_pcs->m_num->ele_gauss_points);

		total_area = 0.;
		for(int i = 0; i<nodes_vector.size(); i++)
		{
			total_area += nodes_area_vector[i];
		}

		tec_file << aktuelle_zeit;

		for(int n = 0; n < _nod_value_vector.size(); n++)
		{
			average_value = 0;
			for(int i = 0; i<nodes_vector.size(); i++)
			{
				average_value += m_pcs->GetNodeValue(nodes_vector[i],
						m_pcs->GetNodeValueIndex(_nod_value_vector[n]) + 1)   // index increased by one - new value taken
								* nodes_area_vector[i];  // only first variable
			}
			average_value /= total_area;
			tec_file << " " << average_value;
		}

		tec_file << "\n";
		_new_file_opened = true;
		tec_file.close();                     // kg44 close file
}

/**************************************************************************
   FEMLib-Method:
   04/2020 JOD Implementation

restrictions:
	to primary variables (index increased by one)
	axisymmety ignored
**************************************************************************/
void COutput::NODWritePLYAverageDataTEC(int time_step_number)
{
		if (_nod_value_vector.size() == 0)
		{
			std::cout << "Warning - No nodes for polyline " << geo_name << "\n";
			return;
		}

		m_pcs = PCSGet(getProcessType());
		if(m_pcs == NULL)
		{
			std::cout << "Warning - PCS not known for polyline-averaged output" << "\n";
			return;
		}

		// File handling
		int number=1;
		char number_char[6];
		sprintf(number_char, "%i", number);
		string number_string(number_char);


		std::string tec_file_name = file_base_name
	     + "_ply_" + geo_name + "_" + std::string(convertProcessTypeToString(getProcessType()))
	    		 + "_averaged" + TEC_FILE_EXTENSION;

		if (!_new_file_opened)
			remove(tec_file_name.c_str());
		else if (time_step_number == 0)  // JOD 2020-05-08
			return;

	   std::fstream tec_file;
	   tec_file.open(tec_file_name.data(), ios::app|ios::out);

		tec_file.setf(ios::scientific, ios::floatfield);
		tec_file.precision(12);
		if (!tec_file.good())
		{
			std::cout << "Warning - Could not open file for writing polyline data " << geo_name << "\n";
			return;
		}
		tec_file.seekg(0L, ios::beg);
	#ifdef SUPERCOMPUTER
		// kg44 buffer the output
		char mybuffer [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
		tec_file.rdbuf()->pubsetbuf(mybuffer,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
		//
	#endif
		//--------------------------------------------------------------------
		// Write header
		if (aktueller_zeitschritt == 0)
		{
			tec_file << "Time "; // << "\n";
			for(int n = 0; n < _nod_value_vector.size(); n++)
				tec_file << _nod_value_vector[n] << " ";
			tec_file << "\n";
		}

		//--------------------------------------------------------------------
		// Write data
		std::vector<long> nodes_vector;
		
		GEOLIB::Polyline const* const ply (
                dynamic_cast<GEOLIB::Polyline const* const>(this->getGeoObj()));
		m_msh->GetNODOnPLY(ply, nodes_vector);

		double total_area, average_value;
		std::vector<double> nodes_area_vector(nodes_vector.size(), 1.);

		if (m_msh->GetMaxElementDim() == 1)
			DomainIntegration(m_pcs, nodes_vector, nodes_area_vector);
		else
			EdgeIntegration(m_msh, nodes_vector, nodes_area_vector,
					FiniteElement::CONSTANT_NEUMANN,
					FiniteElement::convertPrimaryVariable(m_pcs->GetPrimaryVName(0)),  // not used (for BC)
					true,  // ignore axisymmetry
					false  // means no pressure BC
			);


		total_area = 0.;
		for(int i = 0; i<nodes_vector.size(); i++)
		{
			total_area += nodes_area_vector[i];
		}

		tec_file << aktuelle_zeit;

		for(int n = 0; n < _nod_value_vector.size(); n++)
		{
			average_value = 0;
			for(int i = 0; i<nodes_vector.size(); i++)
			{
				average_value += m_pcs->GetNodeValue(nodes_vector[i],
						m_pcs->GetNodeValueIndex(_nod_value_vector[n]) + 1)   // index increased by one - new value taken
								* nodes_area_vector[i];
			}
			average_value /= total_area;
			tec_file << " " << average_value;
		}

		tec_file << "\n";
		_new_file_opened = true;
		tec_file.close();
}



void COutput::GetNodeIndexVector(vector<int>&NodeIndex)
{
	CRFProcess* pcs = NULL;
	const size_t nName = _nod_value_vector.size();
	if (getProcessType() != FiniteElement::INVALID_PROCESS)
	{
		pcs = PCSGet(getProcessType());
		for (size_t k = 0; k < nName; k++)
		{
			if (getProcessType () == FiniteElement::MASS_TRANSPORT)
				pcs = PCSGet(getProcessType(), _nod_value_vector[k]);
			if (!pcs)
			{
				cout
				<< "Warning in COutput::GetNodeIndexVector - no PCS data: "
				<< _nod_value_vector[k] << "\n";
				return;
			}
			NodeIndex[k] = pcs->GetNodeValueIndex(_nod_value_vector[k],true);  // JT latest
		}
	}
	else if (msh_type_name.size() > 0)
	{
		pcs = PCSGet(msh_type_name);
		if (!pcs)
		{
			cout << "Warning in COutput::GetNodeIndexVector - no PCS data"
			     << "\n";
			return;
		}
		for (size_t k = 0; k < nName; k++)
		{
			NodeIndex[k] = pcs->GetNodeValueIndex(_nod_value_vector[k],true); // JT latest
		}
	}
	else if (fem_msh_vector.size() == 1)
	{
		bool bdummy = true;
		for (size_t k = 0; k < nName; k++)
		{
			pcs = PCSGet(_nod_value_vector[k], bdummy);
			if (!pcs)
			{
				cout
				<< "Warning in COutput::GetNodeIndexVector - no PCS data: "
				<< _nod_value_vector[k] << "\n";
				return;
			}
			NodeIndex[k] = pcs->GetNodeValueIndex(_nod_value_vector[k],true); // JT latest
		}
	}
}

CRFProcess* COutput::GetPCS(const string &var_name)
{
	CRFProcess* m_pcs = NULL;
	if(getProcessType () != FiniteElement::INVALID_PROCESS)
		m_pcs = PCSGet(getProcessType());
	else if(msh_type_name.size() > 0)
		m_pcs = PCSGet(msh_type_name);
	if(!m_pcs)
		m_pcs = PCSGet(var_name,true);
	return m_pcs;
}

// 09/2010 TF
CRFProcess* COutput::GetPCS()
{
	if(getProcessType () != FiniteElement::INVALID_PROCESS)
	{
		if (getProcess() != NULL)
			return getProcess();
		else
			return PCSGet(getProcessType());
	}
	else
	{
		CRFProcess* pcs (NULL);
		if(msh_type_name.size() > 0)
			pcs = PCSGet(msh_type_name);
		//	  if(!pcs)
		//		pcs = PCSGet(var_name,true);
		return pcs;
	}
}

CRFProcess* COutput::GetPCS_ELE(const string &var_name)
{
	string pcs_var_name;
	CRFProcess* m_pcs = NULL;
	//----------------------------------------------------------------------
	//  if(pcs_type_name.size()>0)
	//    m_pcs = PCSGet(pcs_type_name);
	if (getProcessType() != FiniteElement::INVALID_PROCESS)
		m_pcs = PCSGet(getProcessType());
	else if (msh_type_name.size() > 0)
		m_pcs = PCSGet(msh_type_name);

	if (!m_pcs)
		for (size_t i = 0; i < pcs_vector.size(); i++)
		{
			m_pcs = pcs_vector[i];
			for (int j = 0; j < m_pcs->pcs_number_of_evals; j++)
			{
				pcs_var_name = m_pcs->pcs_eval_name[j];
				if (pcs_var_name.compare(var_name) == 0)
					return m_pcs;
			}
		}
	return m_pcs;
}

void COutput::GetELEValuesIndexVector(vector<int>&ele_value_index_vector)
{
	if (_ele_value_vector[0].size() == 0)
		return;
   CRFProcess * m_pcs = NULL;

   // CB THMBM
   //m_pcs = GetPCS_ELE(_ele_value_vector[0]);   // CB this is buggy: not all ele vals are defined with the same (or any) process
   for (size_t i = 0; i < _ele_value_vector.size(); i++)
   {
     m_pcs = GetPCS_ELE(_ele_value_vector[i]);   // CB 
     ele_value_index_vector[i] = m_pcs->GetElementValueIndex(_ele_value_vector[i]);
   }
}

/**************************************************************************
   FEMLib-Method:
   Programing:
   01/2006 OK Implementation
   07/2010 TF substituted GEOGetPLYByName, renamed local variables for CGLPolyline,
         CFEMesh and CRFProcess
         (conflicting with class attributes that have the same name)
**************************************************************************/
void COutput::SetNODFluxAtPLY()
{
	// Tests
	//  CGLPolyline* ply = GEOGetPLYByName(geo_name);
	//	CGLPolyline* ply = polyline_vector[getGeoObjIdx()];
	//	if (!ply) {
	//		cout << "Warning in COutput::SetNODFluxAtPLY() - no PLY data" << "\n";
	//		return;
	//	}

	//	CFEMesh* msh = GetMSH();
	if (!m_msh)
	{
		cout << "Warning in COutput::SetNODFluxAtPLY() - no MSH data" << "\n";
		return;
	}

	CRFProcess* pcs = GetPCS("FLUX");
	if (!pcs)
	{
		cout << "Warning in COutput::SetNODFluxAtPLY() - no PCS data" << "\n";
		return;
	}

	std::vector<long> nodes_vector;
	//	msh->GetNODOnPLY(ply, nodes_vector);
	m_msh->GetNODOnPLY(static_cast<const GEOLIB::Polyline*>(getGeoObj()), nodes_vector);
	std::vector<double> node_value_vector;
	node_value_vector.resize(nodes_vector.size());
	int nidx = pcs->GetNodeValueIndex("FLUX");
	for (size_t i = 0; i < nodes_vector.size(); i++)
		node_value_vector[i] = pcs->GetNodeValue(nodes_vector[i], nidx);
	//----------------------------------------------------------------------
	//m_st->FaceIntegration(m_pcs,nodes_vector,node_value_vector);
}

void COutput::ELEWriteSFC_TEC()
{
	//----------------------------------------------------------------------
	if (_ele_value_vector.size() == 0)
		return;
	//----------------------------------------------------------------------
	// File handling
	//......................................................................
	std::string tec_file_name = file_base_name + "_surface" + "_ele";
	addInfoToFileName(tec_file_name, false, true, true);
	//  if(pcs_type_name.size()>1) // PCS
	//    tec_file_name += "_" + msh_type_name;
	//  if(msh_type_name.size()>1) // MSH
	//    tec_file_name += "_" + msh_type_name;
	//  tec_file_name += TEC_FILE_EXTENSION;

	if (!_new_file_opened)
		remove(tec_file_name.c_str());  //WW
	//......................................................................
	fstream tec_file(tec_file_name.data(), ios::app | ios::out);
	tec_file.setf(ios::scientific, ios::floatfield);
	tec_file.precision(12);
	if (!tec_file.good())
		return;
	tec_file.seekg(0L, ios::beg);
#ifdef SUPERCOMPUTER
	// kg44 buffer the output
	char mybuffer [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
	tec_file.rdbuf()->pubsetbuf(mybuffer,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
	//
#endif
	//--------------------------------------------------------------------
	vector<long> tmp_ele_sfc_vector;
	tmp_ele_sfc_vector.clear();
	//--------------------------------------------------------------------
	ELEWriteSFC_TECHeader(tec_file);
	ELEWriteSFC_TECData(tec_file);
	//--------------------------------------------------------------------
	tec_file.close();                     // kg44 close file
	_new_file_opened = true;
}

void COutput::ELEWriteSFC_TECHeader(fstream &tec_file)
{
	// Write Header I: variables
	tec_file << "VARIABLES = \"X\",\"Y\",\"Z\"";
	for (size_t i = 0; i < _ele_value_vector.size(); i++)
		tec_file << "," << _ele_value_vector[i];
	tec_file << "\n";
	//--------------------------------------------------------------------
	// Write Header II: zone
	tec_file << "ZONE T=\"";
	tec_file << _time << "s\", ";
	tec_file << "I=" << (long) m_msh->ele_vector.size() << ", ";
	tec_file << "F=POINT" << ", ";
	tec_file << "C=BLACK";
	tec_file << "\n";
}

void COutput::ELEWriteSFC_TECData(fstream &tec_file)
{
	tec_file << "COutput::ELEWriteSFC_TECData - implementation not finished" << "\n";

	/* // Make it as comment to avoid compilation warnings. 18.082011 WW
	   long i;
	   int j;
	   MeshLib::CElem* m_ele = NULL;
	   MeshLib::CElem* m_ele_neighbor = NULL;
	   double v[3];
	   CRFProcess* m_pcs = NULL;
	   double v_n;
	   //--------------------------------------------------------------------
	   m_pcs = pcs_vector[0];                         //GetPCS_ELE(ele_value_vector[0]);
	   int nidx[3];
	   nidx[0] = m_pcs->GetElementValueIndex("VELOCITY1_X");
	   nidx[1] = m_pcs->GetElementValueIndex("VELOCITY1_Y");
	   nidx[2] = m_pcs->GetElementValueIndex("VELOCITY1_Z");
	   //--------------------------------------------------------------------
	   for(i=0l;i<(long)m_msh->ele_vector.size();i++)
	   {
	   m_ele = m_msh->ele_vector[i];
	   for(j=0;j<m_ele->GetFacesNumber();j++)
	   {
	      m_ele_neighbor = m_ele->GetNeighbor(j);
	      if((m_ele->GetDimension() - m_ele_neighbor->GetDimension())==1)
	      {
	         v[0] = m_pcs->GetElementValue(m_ele->GetIndex(),nidx[0]);
	         v[1] = m_pcs->GetElementValue(m_ele->GetIndex(),nidx[1]);
	         v[2] = m_pcs->GetElementValue(m_ele->GetIndex(),nidx[2]);
	         m_ele_neighbor->SetNormalVector();

	         v_n = v[0]*m_ele_neighbor->normal_vector[0] \
	 + v[1]*m_ele_neighbor->normal_vector[1] \
	 + v[2]*m_ele_neighbor->normal_vector[2];

	      }
	   }
	   }
	 */
	//--------------------------------------------------------------------
}

/**************************************************************************
   FEMLib-Method:
   08/2006 OK Implementation
   modification:
   05/2010 TF if geo_type_name[3] == 'N' (case point) only a pointer to the
   CGLPoint is fetched, but the pointer is never used
   07/2010 TF substituted GEOGetPLYByName
   09/2010 TF substituted pcs_type_name
**************************************************************************/
void COutput::CalcELEFluxes()
{
	double Test[5];

	const FiniteElement::ProcessType pcs_type (getProcessType());
	if (pcs_type == FiniteElement::INVALID_PROCESS)      // WW moved it here.

		//WW cout << "Warning in COutput::CalcELEFluxes(): no PCS data" << "\n";
		return;

	CRFProcess* pcs = PCSGet(getProcessType());
    //BG 04/2011: MASS_TRANSPORT added to get MASS FLUX for Polylines
    //cout << pcs->Tim->step_current << "\n";
    if (isDeformationProcess(pcs_type) || (!isFlowProcess (pcs_type) && (pcs_type != FiniteElement::MASS_TRANSPORT))
	//if (isDeformationProcess(pcs_type) || !isFlowProcess (pcs_type)
	    //WW
	    || pcs->m_msh->geo_name.find("REGIONAL") != string::npos)
		return;

	//----------------------------------------------------------------------
	switch (getGeoType())
	{
	case GEOLIB::POLYLINE:
	{
		//		CGLPolyline* ply = GEOGetPLYByName(geo_name);
		//		if (!ply)
		//			std::cout << "Warning in COutput::CalcELEFluxes - no GEO data" << "\n";

		 //BG 04/2011: ELEWritePLY_TEC does not work for MASS_TRANSPORT because there is no flux considered
		if (pcs_type != FiniteElement::MASS_TRANSPORT)
		{
			double f_n_sum = 0.0;
			double *PhaseFlux(new double [2]);
			std::string Header[2];
			int dimension = 2;
			Header[0] = "q_Phase1";
			Header[1] = "q_Phase2";

			pcs->CalcELEFluxes(static_cast<const GEOLIB::Polyline*> (getGeoObj()), PhaseFlux);
			if ((pcs_type == FiniteElement::GROUNDWATER_FLOW) || (pcs_type == FiniteElement::FLUID_FLOW))
			{
				ELEWritePLY_TEC();
				f_n_sum = PhaseFlux[0];
				TIMValue_TEC(f_n_sum);
			}
			if (pcs_type == FiniteElement::MULTI_PHASE_FLOW)
			{
				Test[0] = PhaseFlux[0];
				Test[1] = PhaseFlux[1];
				TIMValues_TEC(Test, Header, dimension);
			}
			delete [] PhaseFlux;
		}
		// BG, Output for Massflux added
		else
		{
			double *MassFlux (new double[5]);
			std::string Header[5];
			int dimension = 5;
			Header[0] = "AdvectiveMassFlux";
			Header[1] = "DispersiveMassFlux";
			Header[2] = "DiffusiveMassFlux";
			Header[3] = "TotalMassFlux";
			Header[4] = "TotalMass_sum";

			pcs->CalcELEMassFluxes(static_cast<const GEOLIB::Polyline*> (getGeoObj()), geo_name, MassFlux);
			Test[0] = MassFlux[0];
			Test[1] = MassFlux[1];
			Test[2] = MassFlux[2];
			Test[3] = MassFlux[3];
			Test[4] = MassFlux[4];
			TIMValues_TEC(Test, Header, dimension);
			delete [] MassFlux;
		}

		//double f_n_sum = 0.0;
		//		f_n_sum = pcs->CalcELEFluxes(ply); // TF
		//f_n_sum = pcs->CalcELEFluxes(static_cast<const GEOLIB::Polyline*> (getGeoObj()));

		//ELEWritePLY_TEC();
		//BUGFIX_4402_OK_1
		//TIMValue_TEC(f_n_sum);
		break;
	}
	case GEOLIB::SURFACE:
	{
		//		Surface* m_sfc = GEOGetSFCByName(geo_name);
		//		pcs->CalcELEFluxes(m_sfc);
		break;
	}
	case GEOLIB::VOLUME:
	{
		//		CGLVolume* m_vol = GEOGetVOL(geo_name);
		//		pcs->CalcELEFluxes(m_vol);
		break;
	}
	case GEOLIB::GEODOMAIN:               //domAin
		//pcs->CalcELEFluxes(m_dom);
		break;
	default:
		cout << "Warning in COutput::CalcELEFluxes(): no GEO type data" << "\n";
	}

	// WW   pcs->CalcELEFluxes(ply) changed 'mark' of elements
	for (size_t i = 0; i < fem_msh_vector.size(); i++)
		for (size_t j = 0; j < fem_msh_vector[i]->ele_vector.size(); j++)
			fem_msh_vector[i]->ele_vector[j]->MarkingAll(true);
}

void COutput::ELEWritePLY_TEC()
{
	//----------------------------------------------------------------------
	if(_ele_value_vector.size() == 0)
		return;
	//----------------------------------------------------------------------
	// File handling
	//......................................................................
	string tec_file_name = file_base_name; // + "_ply" + "_ele";
	tec_file_name += "_" + getGeoTypeAsString();
	tec_file_name += "_" + geo_name;
	tec_file_name += "_ELE";
	//  if(pcs_type_name.size()>1) // PCS
	//    tec_file_name += "_" + pcs_type_name;
	if(getProcessType () != FiniteElement::INVALID_PROCESS) // PCS
		tec_file_name += "_" + convertProcessTypeToString (getProcessType ());

	if(msh_type_name.size() > 1)          // MSH
		tec_file_name += "_" + msh_type_name;
	tec_file_name += TEC_FILE_EXTENSION;
	if(!_new_file_opened)
		remove(tec_file_name.c_str());  //WW
	//......................................................................
	fstream tec_file(tec_file_name.data(),ios::app | ios::out);
	tec_file.setf(ios::scientific,ios::floatfield);
	tec_file.precision(12);
	if (!tec_file.good())
		return;
	tec_file.seekg(0L,ios::beg);
#ifdef SUPERCOMPUTER
	// kg44 buffer the output
	char mybuffer [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
	tec_file.rdbuf()->pubsetbuf(mybuffer,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
	//
#endif
	//--------------------------------------------------------------------
	vector<long>tmp_ele_ply_vector;
	tmp_ele_ply_vector.clear();
	//--------------------------------------------------------------------
	ELEWritePLY_TECHeader(tec_file);
	ELEWritePLY_TECData(tec_file);
	//--------------------------------------------------------------------
	_new_file_opened = true;
	tec_file.close();                     // kg44 close file
}

void COutput::ELEWritePLY_TECHeader(fstream &tec_file)
{
	// Write Header I: variables
	tec_file << "VARIABLES = \"X\",\"Y\",\"Z\"";
	for (size_t i = 0; i < _ele_value_vector.size(); i++)
		tec_file << "," << _ele_value_vector[i];
	tec_file << "\n";

	// Write Header II: zone
	tec_file << "ZONE T=\"";
	tec_file << _time << "s\", ";
	tec_file << "\n";
}

/**************************************************************************
   FEMLib-Method:
   06/2006 OK Implementation
   07/2010 TF substituted GEOGetPLYByName
   10/2010 TF changed access to process type
**************************************************************************/
void COutput::ELEWritePLY_TECData(fstream &tec_file)
{
	//	CRFProcess* pcs = PCSGet(pcs_type_name);
	CRFProcess* pcs = PCSGet(getProcessType());
	int v_eidx[3];
	CRFProcess* m_pcs_flow = NULL;

	//	if (pcs->pcs_type_name.find("FLOW") != string::npos) {
	if (isFlowProcess(pcs->getProcessType ()))
		m_pcs_flow = pcs;
	else {
		m_pcs_flow = PCSGet(FiniteElement::GROUNDWATER_FLOW);
		if( m_pcs_flow == NULL)
			m_pcs_flow = PCSGet(FiniteElement::LIQUID_FLOW);
	}
	v_eidx[0] = m_pcs_flow->GetElementValueIndex("VELOCITY1_X");
	v_eidx[1] = m_pcs_flow->GetElementValueIndex("VELOCITY1_Y");
	v_eidx[2] = m_pcs_flow->GetElementValueIndex("VELOCITY1_Z");
	for (size_t i = 0; i < 3; i++)
		if (v_eidx[i] < 0)
		{
			cout <<
			"Fatal error in CRFProcess::CalcELEFluxes(CGLPolyline*m_ply) - abort";
			abort();
		}

	//  CGLPolyline* m_ply = GEOGetPLYByName(geo_name);

	// Get elements at GEO
	//	vector<long> ele_vector_at_geo;
	//	m_msh->GetELEOnPLY(m_ply, ele_vector_at_geo);
	std::vector<size_t> ele_vector_at_geo;
	m_msh->GetELEOnPLY(static_cast<const GEOLIB::Polyline*> (getGeoObj()), ele_vector_at_geo, false);

	// helper variables
	Math_Group::vec<MeshLib::CEdge*> ele_edges_vector(15);
	Math_Group::vec<MeshLib::CNode*> edge_nodes(3);
	double edge_mid_vector[3] = {0.0,0.0,0.0};

	for (size_t i = 0; i < ele_vector_at_geo.size(); i++)
	{
		MeshLib::CElem* m_ele = m_msh->ele_vector[ele_vector_at_geo[i]];
		// x,y,z
		m_ele->GetEdges(ele_edges_vector);
		for (size_t j = 0; j < m_ele->GetEdgesNumber(); j++)
		{
			MeshLib::CEdge* m_edg = ele_edges_vector[j];
			if (m_edg->GetMark())
			{
				m_edg->GetNodes(edge_nodes);
				double const* const pnt0(edge_nodes[0]->getData());
				double const* const pnt1(edge_nodes[1]->getData());
				edge_mid_vector[0] = 0.5 * (pnt1[0] + pnt0[0]);
				edge_mid_vector[1] = 0.5 * (pnt1[1] + pnt0[1]);
				edge_mid_vector[2] = 0.5 * (pnt1[2] + pnt0[2]);
			}
		}
		tec_file << edge_mid_vector[0] << " " << edge_mid_vector[1] << " "
		         << edge_mid_vector[2];
		// ele vector values
		tec_file << " " << m_pcs_flow->GetElementValue(m_ele->GetIndex(),
		                                               v_eidx[0]);
		tec_file << " " << m_pcs_flow->GetElementValue(m_ele->GetIndex(),
		                                               v_eidx[1]);
		tec_file << " " << m_pcs_flow->GetElementValue(m_ele->GetIndex(),
		                                               v_eidx[2]);
		//tec_file << " " << pcs->GetElementValue(m_ele->GetIndex(), f_eidx[0]);
		//tec_file << " " << pcs->GetElementValue(m_ele->GetIndex(), f_eidx[1]);
		//tec_file << " " << pcs->GetElementValue(m_ele->GetIndex(), f_eidx[2]);
		tec_file << "\n";
	}
}

void COutput::TIMValue_TEC(double tim_value)
{
	//----------------------------------------------------------------------
	// File handling
	//......................................................................
	fstream tec_file;
	string tec_file_name = file_base_name; // + "_ply" + "_ele";
	tec_file_name += "_" + getGeoTypeAsString();
	tec_file_name += "_" + geo_name;
	tec_file_name += "_TIM";
	//  if(pcs_type_name.size()>1) // PCS
	//    tec_file_name += "_" + pcs_type_name;
	if(getProcessType () != FiniteElement::INVALID_PROCESS) // PCS
		tec_file_name += "_" + convertProcessTypeToString (getProcessType());
	if(msh_type_name.size() > 1)          // MSH
		tec_file_name += "_" + msh_type_name;
	tec_file_name += TEC_FILE_EXTENSION;
	if(!_new_file_opened)
		remove(tec_file_name.c_str());  //WW
	//......................................................................
    if(aktueller_zeitschritt==0)		//BG:04/2011 deletes the content of the file at the start of the simulation
	   tec_file.open(tec_file_name.data(),ios::trunc|ios::out);
    else
		tec_file.open(tec_file_name.data(),ios::app | ios::out);

	tec_file.setf(ios::scientific,ios::floatfield);
	tec_file.precision(12);
	if (!tec_file.good())
		return;
	tec_file.seekg(0L,ios::beg);
#ifdef SUPERCOMPUTER
	// kg44 buffer the output
	char mybuffer [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
	tec_file.rdbuf()->pubsetbuf(mybuffer,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
	//
#endif
	//--------------------------------------------------------------------
	// Write Header I: variables
    if(aktueller_zeitschritt==0)		//BG:04/2011 bevor it was timestep 1
	{
		tec_file << "VARIABLES = \"Time\",\"Value\"";
		tec_file << "\n";
		//--------------------------------------------------------------------
		// Write Header II: zone
		tec_file << "ZONE T=";
		tec_file << geo_name;
		tec_file << "\n";
	}
	//--------------------------------------------------------------------
	tec_file << aktuelle_zeit << " " << tim_value << "\n";
	//--------------------------------------------------------------------
	tec_file.close();                     // kg44 close file
}

/*-------------------------------------------------------------------------
   GeoSys - Function: TIMValues_TEC
   Task: Can write several values over time
   Return: nothing
   Programming: 10/2011 BG
   Modification:
 -------------------------------------------------------------------------*/
void COutput::TIMValues_TEC(double tim_value[5], std::string *header, int dimension)
{
	double j[10];

    for (int i = 0; i < dimension; i++)
		j[i] = tim_value[i];

   //----------------------------------------------------------------------
   // File handling
   //......................................................................
   fstream tec_file;
   string tec_file_name = file_base_name;         // + "_ply" + "_ele";
   tec_file_name += "_" + getGeoTypeAsString();
   tec_file_name += "_" + geo_name;
   tec_file_name += "_TIM";
   //  if(pcs_type_name.size()>1) // PCS
   //    tec_file_name += "_" + pcs_type_name;
   if(getProcessType () != FiniteElement::INVALID_PROCESS)       // PCS
      tec_file_name += "_" + convertProcessTypeToString (getProcessType());
   if(msh_type_name.size()>1)                     // MSH
      tec_file_name += "_" + msh_type_name;
   tec_file_name += TEC_FILE_EXTENSION;

   if(!_new_file_opened)
      remove(tec_file_name.c_str());              //WW
   //......................................................................
   if(aktueller_zeitschritt==0)		//BG:04/2011 deletes the content of the file at the start of the simulation
	   tec_file.open(tec_file_name.data(),ios::trunc|ios::out);
   else
      tec_file.open(tec_file_name.data(),ios::app|ios::out);
   tec_file.setf(ios::scientific,ios::floatfield);
   tec_file.precision(12);
   if (!tec_file.good()) return;
   tec_file.seekg(0L,ios::beg);
#ifdef SUPERCOMPUTER
   // kg44 buffer the output
   char mybuffer [MY_IO_BUFSIZE*MY_IO_BUFSIZE];
   tec_file.rdbuf()->pubsetbuf(mybuffer,MY_IO_BUFSIZE*MY_IO_BUFSIZE);
   //
#endif
   //--------------------------------------------------------------------
   // Write Header I: variables
   if(aktueller_zeitschritt==0)		//BG:04/2011 bevor it was timestep 1
   {
	  tec_file << "VARIABLES = \"Time\"";
	  for (int i = 0; i < dimension; i++)
         tec_file << ",\"" << header[i] << "\"";
      tec_file << "\n";
      //--------------------------------------------------------------------
      // Write Header II: zone
      tec_file << "ZONE T=";
      tec_file << geo_name;
      tec_file << "\n";
   }
   //--------------------------------------------------------------------
   tec_file << aktuelle_zeit;
   for (int i = 0; i < dimension; i++)
      tec_file << " " << j[i];
   tec_file << "\n";
   _new_file_opened = true;
   //--------------------------------------------------------------------
   tec_file.close();                              // kg44 close file

}

double COutput::NODFlux(long nod_number)
{
	nod_number = nod_number;              //OK411
	/*
	   cout << gnode << " " \
	    << m_pcs->GetNodeValue(gnode,NodeIndex[k]) << end
	   flux_sum += m_pcs->GetNodeValue(gnode,NodeIndex[k]);
	 */
	// All elements at node //OK
#if defined (USE_PETSC) // || defined (other parallel solver lib). 04.2012 WW
	return 0;
#elif defined(NEW_EQS)                                 //WW. 07.11.2008
	return 0.;                            //To do: m_pcs->eqs_new->b[nod_number];
#else
	// Element nodal RHS contributions
	m_pcs->getEQSPointer()->b[nod_number] = 0.0;
	MeshLib::CNode* nod = m_msh->nod_vector[nod_number];
	m_pcs->AssembleParabolicEquationRHSVector(nod);
	return m_pcs->getEQSPointer()->b[nod_number];
#endif
}

void COutput::NODWriteLAYDataTEC(int time_step_number)
{
	// Tests
	const size_t nName (_nod_value_vector.size());
	if (nName == 0)
		return;
	std::vector<int> NodeIndex(nName);

	// PCS
	CRFProcess* m_pcs = PCSGet(getProcessType());
	if (!m_pcs)
		return;
	for (size_t k = 0; k < nName; k++)
		NodeIndex[k] = m_pcs->GetNodeValueIndex(_nod_value_vector[k]);

	// MSH
	//	m_msh = GetMSH();
	if (!m_msh)
	{
		cout << "Warning in COutput::NODWriteLAYDataTEC() - no MSH data"
		     << "\n";
		return;
	}

	// File name handling
	char char_time_step_number[10];
	sprintf(char_time_step_number, "%i", time_step_number);
	string tec_file_name(file_base_name);
	tec_file_name += "_layer_";
	tec_file_name += char_time_step_number;
	tec_file_name += TEC_FILE_EXTENSION;
	fstream tec_file(tec_file_name.data(), ios::trunc | ios::out);

	tec_file.setf(ios::scientific, ios::floatfield);
	tec_file.precision(12);
	if (!tec_file.good())
		return;
#ifdef SUPERCOMPUTER
	//
	// kg44 buffer the output
	char mybuffer [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
	tec_file.rdbuf()->pubsetbuf(mybuffer,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
	//
#endif
	//--------------------------------------------------------------------
	// Write Header I: variables
	tec_file << "VARIABLES = X,Y,Z,N";
	for (size_t k = 0; k < nName; k++)
		tec_file << "," << _nod_value_vector[k] << " ";
	tec_file << "\n";

	long j;
	long no_per_layer = m_msh->GetNodesNumber(false)
	                    / (m_msh->getNumberOfMeshLayers() + 1);
	long jl;
	for (size_t l = 0; l < m_msh->getNumberOfMeshLayers() + 1; l++)
	{
		//--------------------------------------------------------------------
		tec_file << "ZONE T=LAYER" << l << "\n";
		//--------------------------------------------------------------------
		for (j = 0l; j < no_per_layer; j++)
		{
			jl = j + j* m_msh->getNumberOfMeshLayers() + l;
			//..................................................................
			// XYZ
			double const* const pnt (m_msh->nod_vector[jl]->getData());
			tec_file << pnt[0] << " ";
			tec_file << pnt[1] << " ";
			tec_file << pnt[2] << " ";
			tec_file << jl << " ";
			//..................................................................
			for (size_t k = 0; k < nName; k++)
				tec_file << m_pcs->GetNodeValue(
				        m_msh->nod_vector[jl]->GetIndex(), NodeIndex[k]) << " ";
			tec_file << "\n";
		}
	}
	tec_file.close();                     // kg44 close file
}

/**************************************************************************
   FEMLib-Method:
   Task: Write PCON data for ChemApp output
   Programing:
   08/2008 MX Implementation
**************************************************************************/
void COutput::PCONWriteDOMDataTEC()
{
	int te = 0;
	string eleType;
	string tec_file_name;
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL)
	char tf_name[10];
	std::cout << "Process " << myrank << " in WriteDOMDataTEC" << "\n";
#endif

	//----------------------------------------------------------------------
	// Tests
	if(_pcon_value_vector.size() == 0)
		return;
	//......................................................................
	// MSH
	//m_msh = FEMGet(pcs_type_name);
	//  m_msh = GetMSH();
	if(!m_msh)
	{
		cout << "Warning in COutput::NODWriteDOMDataTEC() - no MSH data" << "\n";
		return;
	}
	//======================================================================
	vector<int> mesh_type_list;           //NW
	if (m_msh->getNumberOfLines() > 0)
		mesh_type_list.push_back(1);
	if (m_msh->getNumberOfQuads () > 0)
		mesh_type_list.push_back(2);
	if (m_msh->getNumberOfHexs () > 0)
		mesh_type_list.push_back(3);
	if (m_msh->getNumberOfTris () > 0)
		mesh_type_list.push_back(4);
	if (m_msh->getNumberOfTets () > 0)
		mesh_type_list.push_back(5);
	if (m_msh->getNumberOfPrisms () > 0)
		mesh_type_list.push_back(6);

	// Output files for each mesh type
	//NW
	for (int i = 0; i < (int)mesh_type_list.size(); i++)
	{
		te = mesh_type_list[i];
		//----------------------------------------------------------------------
		// File name handling
		tec_file_name = file_base_name + "_" + "domain_PCON";
		if(msh_type_name.size() > 0) // MultiMSH
			tec_file_name += "_" + msh_type_name;
		//  if(pcs_type_name.size()>0) // PCS
		//    tec_file_name += "_" + pcs_type_name;
		if(getProcessType () != FiniteElement::INVALID_PROCESS) // PCS
			tec_file_name += "_" + convertProcessTypeToString (getProcessType());
		//======================================================================
		switch (te)               //NW
		{
		case 1:
			tec_file_name += "_line";
			eleType = "QUADRILATERAL";
			break;
		case 2:
			tec_file_name += "_quad";
			eleType = "QUADRILATERAL";
			break;
		case 3:
			tec_file_name += "_hex";
			eleType = "BRICK";
			break;
		case 4:
			tec_file_name += "_tri";
			eleType = "QUADRILATERAL";
			break;
		case 5:
			tec_file_name += "_tet";
			eleType = "TETRAHEDRON";
			break;
		case 6:
			tec_file_name += "_pris";
			eleType = "BRICK";
			break;
		}
		/*
		   if(m_msh->msh_no_line>0)
		   {
		      tec_file_name += "_line";
		      eleType = "QUADRILATERAL";
		     te=1;
		   }
		   else if (m_msh->msh_no_quad>0)
		   {
		      tec_file_name += "_quad";
		      eleType = "QUADRILATERAL";
		   te=2;
		   }
		   else if (m_msh->msh_no_hexs>0)
		   {
		   tec_file_name += "_hex";
		   eleType = "BRICK";
		   te=3;
		   }
		   else if (m_msh->msh_no_tris>0)
		   {
		   tec_file_name += "_tri";
		   //???Who was this eleType = "TRIANGLE";
		   eleType = "QUADRILATERAL";
		   te=4;
		   }
		   else if (m_msh->msh_no_tets>0)
		   {
		   tec_file_name += "_tet";
		   eleType = "TETRAHEDRON";
		   te=5;
		   }
		   else if (m_msh->msh_no_pris>0)
		   {
		   tec_file_name += "_pris";
		   eleType = "BRICK";
		   te=6;
		   }
		 */
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL)
		sprintf(tf_name, "%d", myrank);
		tec_file_name += "_" + string(tf_name);
		std::cout << "Tecplot filename: " << tec_file_name << "\n";
#endif
		tec_file_name += TEC_FILE_EXTENSION;
		//WW
		if (!_new_file_opened)
			remove(tec_file_name.c_str());

		fstream tec_file (tec_file_name.data(),ios::app | ios::out);
		tec_file.setf(ios::scientific,ios::floatfield);
		tec_file.precision(12);
		if (!tec_file.good())
			return;
#ifdef SUPERCOMPUTER
		// kg44 buffer the output
		char mybuf1 [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
		tec_file.rdbuf()->pubsetbuf(mybuf1,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
#endif
		//
		WriteTECHeader(tec_file,te,eleType);
		WriteTECNodePCONData(tec_file);
		WriteTECElementData(tec_file,te);
		tec_file.close();         // kg44 close file
		//--------------------------------------------------------------------
		// tri elements
		// ***** 07/2010 TF commented out block since the global variable is always zero
		//    if(msh_no_tris>0){
		//    //string tec_file_name = pcs_type_name + "_" + "domain" + "_tri" + TEC_FILE_EXTENSION;
		//#ifdef SUPERCOMPUTER
		//// buffer the output
		//      char sxbuf1[MY_IO_BUFSIZE*MY_IO_BUFSIZE];
		//#endif
		//      string tec_file_name = file_base_name + "_" + "domain" + "_tri" + TEC_FILE_EXTENSION;
		//      fstream tec_file1 (tec_file_name.data(),ios::app|ios::out);
		//      tec_file1.setf(ios::scientific,ios::floatfield);
		//      tec_file1.precision(12);
		//      if (!tec_file1.good()) return;
		//#ifdef SUPERCOMPUTER
		//      tec_file1.rdbuf()->pubsetbuf(sxbuf1,MY_IO_BUFSIZE*MY_IO_BUFSIZE);
		//#endif
		//      //OK  tec_file1.clear();
		//      //OK  tec_file1.seekg(0L,ios::beg);
		//      WriteTECHeader(tec_file1,4,"TRIANGLE");
		//      WriteTECNodeData(tec_file1);
		//      WriteTECElementData(tec_file1,4);
		//      tec_file1.close(); // kg44 close file
		//    }
		//--------------------------------------------------------------------
		// ***** 07/2010 TF commented out block since the global variable is always zero
		// quad elements
		//    if(msh_no_quad>0){
		//      //string tec_file_name = pcs_type_name + "_" + "domain" + "_quad" + TEC_FILE_EXTENSION;
		//#ifdef SUPERCOMPUTER
		//      char sxbuf2[MY_IO_BUFSIZE*MY_IO_BUFSIZE];
		//#endif
		//      string tec_file_name = file_base_name + "_" + "domain" + "_quad" + TEC_FILE_EXTENSION;
		//      fstream tec_file (tec_file_name.data(),ios::app|ios::out);
		//
		//      tec_file.setf(ios::scientific,ios::floatfield);
		//      tec_file.precision(12);
		//      if (!tec_file.good()) return;
		//#ifdef SUPERCOMPUTER
		//      tec_file.rdbuf()->pubsetbuf(sxbuf2,MY_IO_BUFSIZE*MY_IO_BUFSIZE);
		//#endif
		//      WriteTECHeader(tec_file,2,"QUADRILATERAL");
		//      WriteTECNodeData(tec_file);
		//      WriteTECElementData(tec_file,2);
		//      tec_file.close(); // kg44 close file
		//    }
		//--------------------------------------------------------------------
		// ***** 07/2010 TF commented out block since the global variable is always zero
		// tet elements
		//    if(msh_no_tets>0){
		//      //string tec_file_name = pcs_type_name + "_" + "domain" + "_tet" + TEC_FILE_EXTENSION;
		//#ifdef SUPERCOMPUTER
		//      char sxbuf3[MY_IO_BUFSIZE*MY_IO_BUFSIZE];
		//#endif
		//
		//      string tec_file_name = file_base_name + "_" + "domain" + "_tet";
		//
		//#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL)
		//      sprintf(tf_name, "%d", myrank);
		//      tec_file_name += "_" + string(tf_name);
		//#endif
		//
		//      tec_file_name += TEC_FILE_EXTENSION;
		//
		//      fstream tec_file (tec_file_name.data(),ios::app|ios::out);
		//
		//      tec_file.setf(ios::scientific,ios::floatfield);
		//      tec_file.precision(12);
		//      if (!tec_file.good()) return;
		//#ifdef SUPERCOMPUTER
		//      tec_file.rdbuf()->pubsetbuf(sxbuf3,MY_IO_BUFSIZE*MY_IO_BUFSIZE);
		//#endif
		//
		//      WriteTECHeader(tec_file,5,"TETRAHEDRON");
		//      WriteTECNodeData(tec_file);
		//      WriteTECElementData(tec_file,5);
		//      tec_file.close(); // kg44 close file
		//    }
		//--------------------------------------------------------------------
		// ***** 07/2010 TF commented out block since the global variable is always zero
		// pris elements
		//    if(msh_no_pris>0){
		//      //string tec_file_name = pcs_type_name + "_" + "domain" + "_pris" + TEC_FILE_EXTENSION;
		//#ifdef SUPERCOMPUTER
		//        char sxbuf4[MY_IO_BUFSIZE*MY_IO_BUFSIZE];
		//#endif
		//      string tec_file_name = file_base_name + "_" + "domain" + "_pris" + TEC_FILE_EXTENSION;
		//      fstream tec_file (tec_file_name.data(),ios::app|ios::out);
		//
		//      tec_file.setf(ios::scientific,ios::floatfield);
		//      tec_file.precision(12);
		//      if (!tec_file.good()) return;
		//#ifdef SUPERCOMPUTER
		//      tec_file.rdbuf()->pubsetbuf(sxbuf4,MY_IO_BUFSIZE*MY_IO_BUFSIZE);
		//#endif
		//
		//      WriteTECHeader(tec_file,6,"BRICK");
		//      WriteTECNodeData(tec_file);
		//      WriteTECElementData(tec_file,6);
		//      tec_file.close(); // kg44 close file
		//    }
		//--------------------------------------------------------------------
		// ***** 07/2010 TF commented out block since the global variable is always zero
		// hex elements
		//    if(msh_no_hexs>0){
		//      //string tec_file_name = pcs_type_name + "_" + "domain" + "_hex" + TEC_FILE_EXTENSION;
		//#ifdef SUPERCOMPUTER
		//        char sxbuf5[MY_IO_BUFSIZE*MY_IO_BUFSIZE];
		//#endif
		//
		//      string tec_file_name = file_base_name + "_" + "domain" + "_hex" + TEC_FILE_EXTENSION;
		//      fstream tec_file (tec_file_name.data(),ios::app|ios::out);
		//
		//
		//      tec_file.setf(ios::scientific,ios::floatfield);
		//      tec_file.precision(12);
		//      if (!tec_file.good()) return;
		//#ifdef SUPERCOMPUTER
		//      tec_file.rdbuf()->pubsetbuf(sxbuf5,MY_IO_BUFSIZE*MY_IO_BUFSIZE);
		//#endif
		//      WriteTECHeader(tec_file,3,"BRICK");
		//      WriteTECNodeData(tec_file);
		//      WriteTECElementData(tec_file,3);
		//      tec_file.close(); // kg44 close file
		//    }
	}
	_new_file_opened = true;
}

/**************************************************************************
   FEMLib-Method:
   Task: Node value output of PCON in aquous
   Programing:
   08/2008 MX Implementation
**************************************************************************/
void COutput::WriteTECNodePCONData(fstream &tec_file)
{
	const size_t nName (_pcon_value_vector.size());
	int nidx_dm[3];
	std::vector<int> PconIndex(nName);

	//  m_msh = GetMSH();

#ifdef CHEMAPP
	CEqlink* eq = NULL;

	eq = eq->GetREACTION();
	if (!eq)
		return;
	const int nPCON_aq = eq->NPCON[1];    //GetNPCON(1);
	eq->GetPconNameAq();

	for(i = 0; i < nName; i++)
	{
		for(k = 0; k < nPCON_aq; k++)
			//		 pcon_value_name = PconName_Aq[i];
			if(pcon_value_vector[i].compare(PconName_Aq[k]) == 0)
			{
				PconIndex[i] = k;
				break;
			}
	}
#endif
	// MSH
	//--------------------------------------------------------------------
	for (size_t j = 0l; j < m_msh->GetNodesNumber(false); j++)
	{
		// XYZ
		double x[3] =
		{m_msh->nod_vector[j]->getData()[0], m_msh->nod_vector[j]->getData()[1],
		 m_msh->nod_vector[j]->getData()[2]};
//      x[0] = m_msh->nod_vector[j]->X();
//      x[1] = m_msh->nod_vector[j]->Y();
//      x[2] = m_msh->nod_vector[j]->Z();
		// Amplifying DISPLACEMENTs
		if (M_Process || MH_Process) //WW

			for (size_t k = 0; k < max_dim + 1; k++)
				x[k] += out_amplifier * m_pcs->GetNodeValue(
				        m_msh->nod_vector[j]->GetIndex(), nidx_dm[k]);
		for (size_t i = 0; i < 3; i++)
			tec_file << x[i] << " ";
		// NOD values
#ifdef CHEMAPP
		for(size_t k = 0; k < nName; k++)
			tec_file << eq->GetPconAq_mol_amount(j,PconIndex[k]) << " ";

#endif
		tec_file << "\n";
	}
}

void COutput::checkConsistency ()
{
	if (!_nod_value_vector.empty())
	{
		std::vector<std::string> del_index, alias_del_lindex;
		bool found = false;
		CRFProcess* pcs = NULL;
		const size_t pcs_vector_size(pcs_vector.size());
		const size_t nod_value_vector_size (_nod_value_vector.size());
		for (size_t j = 0; j < nod_value_vector_size; j++)
		{
			found = false; // initialize variable found
			for (size_t l = 0; l < pcs_vector_size; l++)
			{
				pcs = pcs_vector[l];
				for (size_t m = 0; m < pcs->nod_val_name_vector.size(); m++)
				{
					if (pcs->nod_val_name_vector[m].compare(
						_nod_value_vector[j]) == 0)
					{
						found = true;
						del_index.push_back(_nod_value_vector[j]);
						alias_del_lindex.push_back(_alias_nod_value_vector[j]);
						break;
					}
				}
				// end for(m...)
			}             // end for(l...)
			if (!found)
			{
				std::cout << "Warning - no PCS data for output variable "
				          << _nod_value_vector[j] << " in ";
				switch (getGeoType())
				{
				case GEOLIB::POINT:
					std::cout << "POINT " << getGeoName() << "\n";
					break;
				case GEOLIB::POLYLINE:
					std::cout << "POLYLINE " << getGeoName() << "\n";
					break;
				case GEOLIB::SURFACE:
					std::cout << "SURFACE " << getGeoName() << "\n";
					break;
				case GEOLIB::VOLUME:
					std::cout << "VOLUME " << getGeoName() << "\n";
					break;
				case GEOLIB::GEODOMAIN:
					std::cout << "DOMAIN " << getGeoName() << "\n";
					break;
				case GEOLIB::INVALID:
					std::cout <<
					"WARNING: COutput::checkConsistency - invalid geo type" <<
					"\n";
					break;
				}
			}
		}                         // end for(j...)

		// Reduce vector out->_nod_value_vector by elements which have no PCS
		if (del_index.size() < _nod_value_vector.size())
		{
			std::cout << " Reducing output to variables with existing PCS-data " <<
			"\n";
			_nod_value_vector.clear();
			for (size_t j = 0; j < del_index.size(); j++)
				_nod_value_vector.push_back(del_index[j]);
            _alias_nod_value_vector.clear();
            for (size_t j = 0; j < del_index.size(); j++)
                _alias_nod_value_vector.push_back(alias_del_lindex[j]);
		}
		if (!pcs)
			pcs = this->GetPCS();
		if (!pcs)
			cout << "Warning in OUTData - no PCS data" << "\n";
	}                                     // end if(_nod_value_vector.size()>0)
}

/**************************************************************************
   FEMLib-Method:
   Task: Set output variable names for internal use
   Programing:
   11/2011 NW Implementation
**************************************************************************/
void COutput::setInternalVarialbeNames(CFEMesh *msh)
{
#if 0
    if (_alias_nod_value_vector.empty())
        return;
    bool isXZplane = (msh->GetCoordinateFlag()==22);
    bool isPVD = (dat_type_name.compare("PVD") == 0); //currently only for PVD

    if (isXZplane && isPVD) {
        std::cout << "-> recognized XZ plane for PVD output." << "\n";
        map<string,string> map_output_variable_name;
        map_output_variable_name.insert(pair<string, string>("DISPLACEMENT_Y1", "DISPLACEMENT_Z1" ));
        map_output_variable_name.insert(pair<string, string>("DISPLACEMENT_Z1", "DISPLACEMENT_Y1" ));
        map_output_variable_name.insert(pair<string, string>("STRESS_XY", "STRESS_XZ" ));
        map_output_variable_name.insert(pair<string, string>("STRESS_YY", "STRESS_ZZ" ));
        map_output_variable_name.insert(pair<string, string>("STRESS_ZZ", "STRESS_YY" ));
        map_output_variable_name.insert(pair<string, string>("STRESS_XZ", "STRESS_XY" ));
        map_output_variable_name.insert(pair<string, string>("STRAIN_XY", "STRAIN_XZ" ));
        map_output_variable_name.insert(pair<string, string>("STRAIN_YY", "STRAIN_ZZ" ));
        map_output_variable_name.insert(pair<string, string>("STRAIN_ZZ", "STRAIN_YY" ));
        map_output_variable_name.insert(pair<string, string>("STRAIN_XZ", "STRAIN_XY"  ));
        map_output_variable_name.insert(pair<string, string>("VELOCITY_Y1", "VELOCITY_Z1"));
        map_output_variable_name.insert(pair<string, string>("VELOCITY_Z1", "VELOCITY_Y1"));
        map_output_variable_name.insert(pair<string, string>("VELOCITY_Y2", "VELOCITY_Z2"));
        map_output_variable_name.insert(pair<string, string>("VELOCITY_Z2", "VELOCITY_Y2"));

        for (size_t j = 0; j < _alias_nod_value_vector.size(); j++) {
            if (map_output_variable_name.count(_alias_nod_value_vector[j])>0) {
                _nod_value_vector.push_back(map_output_variable_name[_alias_nod_value_vector[j]]);
            } else {
                _nod_value_vector.push_back(_alias_nod_value_vector[j]);
            }
        }
    } else {
        for (size_t j = 0; j < _alias_nod_value_vector.size(); j++) {
            _nod_value_vector.push_back(_alias_nod_value_vector[j]);
        }
    }
#else
    if (_nod_value_vector.empty())
        return;
    bool isXZplane = (msh->GetCoordinateFlag()==22);
    bool isPVD = (dat_type_name.compare("PVD") == 0); //currently only for PVD

    if (isXZplane && isPVD) {
        std::cout << "-> recognized XZ plane for PVD output." << "\n";
        map<string,string> map_output_variable_name;
        map_output_variable_name.insert(pair<string, string>("DISPLACEMENT_Y1", "DISPLACEMENT_Z1" ));
        map_output_variable_name.insert(pair<string, string>("DISPLACEMENT_Z1", "DISPLACEMENT_Y1" ));
        map_output_variable_name.insert(pair<string, string>("STRESS_XY", "STRESS_XZ" ));
        map_output_variable_name.insert(pair<string, string>("STRESS_YY", "STRESS_ZZ" ));
        map_output_variable_name.insert(pair<string, string>("STRESS_ZZ", "STRESS_YY" ));
        map_output_variable_name.insert(pair<string, string>("STRESS_XZ", "STRESS_XY" ));
        map_output_variable_name.insert(pair<string, string>("STRAIN_XY", "STRAIN_XZ" ));
        map_output_variable_name.insert(pair<string, string>("STRAIN_YY", "STRAIN_ZZ" ));
        map_output_variable_name.insert(pair<string, string>("STRAIN_ZZ", "STRAIN_YY" ));
        map_output_variable_name.insert(pair<string, string>("STRAIN_XZ", "STRAIN_XY"  ));
        map_output_variable_name.insert(pair<string, string>("VELOCITY_Y1", "VELOCITY_Z1"));
        map_output_variable_name.insert(pair<string, string>("VELOCITY_Z1", "VELOCITY_Y1"));
        map_output_variable_name.insert(pair<string, string>("VELOCITY_Y2", "VELOCITY_Z2"));
        map_output_variable_name.insert(pair<string, string>("VELOCITY_Z2", "VELOCITY_Y2"));

        for (size_t j = 0; j < _nod_value_vector.size(); j++) {
            if (map_output_variable_name.count(_nod_value_vector[j])>0) {
                _alias_nod_value_vector.push_back(map_output_variable_name[_nod_value_vector[j]]);
            } else {
                _alias_nod_value_vector.push_back(_nod_value_vector[j]);
            }
        }
    } else {
        for (size_t j = 0; j < _nod_value_vector.size(); j++) {
            _alias_nod_value_vector.push_back(_nod_value_vector[j]);
        }
    }
#endif
}

void COutput::addInfoToFileName (std::string& file_name, bool geo, bool process, bool mesh) const
{
	// add geo type name
	if (geo)
		//		file_name += getGeoTypeAsString();
		file_name += geo_name;

	// add process type name
	if (getProcessType() != FiniteElement::INVALID_PROCESS && process)
		file_name += "_" + FiniteElement::convertProcessTypeToString (getProcessType());

	// add mesh type name
	if (msh_type_name.size() > 0 && mesh)
		file_name += "_" + msh_type_name;

//#if defined(USE_PETSC)  // JOD 2015-11-17
//	file_name += "_" + mrank_str;
//#endif

	// finally add file extension
	file_name += TEC_FILE_EXTENSION;
}


	
/**************************************************************************
FEMLib-Method:
Task:   Write output of multiple points in single file
Use:    Specify  $DAT_TYPE as COMBINE_POINTS
Programing:
06/2012 JOD Implementation
**************************************************************************/
void COutput::NODWritePointsCombined(double time_current, int time_step_number)
{

	CFEMesh* m_msh = NULL;
	m_msh = FEMGet(convertProcessTypeToString(getProcessType()));
	CRFProcess* m_pcs_out = NULL;
	m_pcs_out = PCSGet(getProcessType());
	bool output = false;
	//std::string tec_file_name(file_base_name + "_time_");
	//addInfoToFileName(tec_file_name, true, true, true);

	char number_char[3];
	string number_string = number_char;
	string tec_file_name = file_base_name + "_" + convertProcessTypeToString(getProcessType()) + "_time_" + "POINTS"; 

	if (!_new_file_opened)
		remove(tec_file_name.c_str());
	else if (time_step_number == 0)  // JOD 2020-05-08
		return;

	if (time_vector.size() == 0 && (nSteps > 0) && (time_step_number
		% nSteps == 0))
		output = true;

	for (size_t j = 0; j < time_vector.size(); j++)
	if ((fabs(time_current - time_vector[j])) < MKleinsteZahl) //WW MKleinsteZahl
		output = true;

	if (output == true)
	{
		fstream tec_file(tec_file_name.data(), ios::app | ios::out);
		tec_file.setf(ios::scientific, ios::floatfield);
		tec_file.precision(12);
		if (!tec_file.good()) return;
		tec_file.seekg(0L, ios::beg);

	long msh_node_number(m_msh->GetNODOnPNT(static_cast<const GEOLIB::Point*> (getGeoObj())));

	//----------------------------------------------------------------------
	// NIDX for output variables
	size_t no_variables(_nod_value_vector.size());
	vector<int> NodeIndex(no_variables);
	GetNodeIndexVector(NodeIndex);

		tec_file << geo_name << " ";
		std::string nod_value_name;
		int nidx = m_pcs_out->GetNodeValueIndex(nod_value_name) + 1;
		double val_n;

	for (size_t i = 0; i < _nod_value_vector.size(); i++) {
		nod_value_name = _nod_value_vector[i];
		val_n = m_pcs_out->GetNodeValue(msh_node_number, NodeIndex[i]);
		tec_file << "time " << time_current << " " << nod_value_name << " " << val_n << " "
						<< "\n";
	}

		tec_file.close();
		_new_file_opened = true;
	}
}

/**************************************************************************
FEMLib-Method:
Task:   
Use:
Programing:
10/2014 JOD Implementation
**************************************************************************/

void COutput::NODWritePrimaryVariableList(double time_current, int time_step_number)
{


	CFEMesh* m_msh = NULL;
	m_msh = FEMGet(convertProcessTypeToString(getProcessType()));
	CRFProcess* m_pcs_out = NULL;
	if ((m_pcs_out = PCSGet(getProcessType())) == 0)
	{
		cout << "error in NODWritePrimaryVariableList - Process " << convertProcessTypeToString(getProcessType()) << " not known" << endl;
		return;
	}

	vector<long> nodes_vector;
	//////////////
	bool output = false;
	char number_char[3];
	string number_string = number_char;
	if (geo_name.size() == 0)
		geo_name = "domain";
	string tec_file_name = file_base_name + "_" + convertProcessTypeToString(getProcessType()) + "_" + geo_name + "_primary_variables";


#if defined(USE_PETSC) 
	tec_file_name += "_" + mrank_str;
#endif

	tec_file_name += ".txt";

	//if (time_step_number == 0)
	//{
		//remove(tec_file_name.c_str());
		//return;
	//}
	//else {

		if (time_vector.size() == 0 && (nSteps > 0) && (time_step_number
			% nSteps == 0))
			output = true;

		for (size_t j = 0; j < time_vector.size(); j++)
		if ((fabs(time_current - time_vector[j])) < MKleinsteZahl) //WW MKleinsteZahl
			output = true;

		if (output == true)
		{
			remove(tec_file_name.c_str());
			fstream tec_file(tec_file_name.data(), ios::app | ios::out);
			tec_file.setf(ios::scientific, ios::floatfield);
			tec_file.precision(12);
			if (!tec_file.good()) return;
			tec_file.seekg(0L, ios::beg);
			//--------------------------------------------------------------------
			Surface *m_sfc = NULL;
			CGLPolyline* m_polyline = NULL;
			GEOLIB::Polyline const* const ply(
				dynamic_cast<GEOLIB::Polyline const* const> (this->getGeoObj()));

		//tec_file << "TIME " << time_current << "\n";

		switch (getGeoType()) {

		case GEOLIB::GEODOMAIN:

				for (long i = 0; i < (long)m_msh->nod_vector.size(); i++)
				if (convertProcessTypeToString(getProcessType()) == "MULTI_COMPONENTIAL_FLOW")  // PRESSURE1, TEMPERATURE1, CONCENTRATION1 
						tec_file << m_msh->nod_vector[i]->GetIndex() << "        " << m_pcs_out->GetNodeValue(m_msh->nod_vector[i]->GetIndex(), 1) << "        " << m_pcs_out->GetNodeValue(m_msh->nod_vector[i]->GetIndex(), 3)  << "        " << m_pcs_out->GetNodeValue(m_msh->nod_vector[i]->GetIndex(), 5) << "\n";
					else
						tec_file << m_msh->nod_vector[i]->GetIndex() << "        " << m_pcs_out->GetNodeValue(m_msh->nod_vector[i]->GetIndex(), 1) << "\n";

			cout << "Data output: " << convertProcessTypeToString(getProcessType()) << " primary variables - DOMAIN - " << m_msh->nod_vector.size() << " nodes" << endl;
			break;
		case GEOLIB::SURFACE:

			m_sfc = GEOGetSFCByName(geo_name);
			if (m_sfc)
				m_msh->GetNODOnSFC(m_sfc, nodes_vector);

				for (long i = 0; i < (long)nodes_vector.size(); i++)
				if (convertProcessTypeToString(getProcessType()) == "MULTI_COMPONENTIAL_FLOW")  // PRESSURE1, TEMPERATURE1, CONCENTRATION1 
					tec_file << m_msh->nod_vector[i]->GetIndex() << "        " << m_pcs_out->GetNodeValue(m_msh->nod_vector[i]->GetIndex(), 1) << "        " << m_pcs_out->GetNodeValue(m_msh->nod_vector[i]->GetIndex(), 3) << "        " << m_pcs_out->GetNodeValue(m_msh->nod_vector[i]->GetIndex(), 5) << "\n";
				else
					tec_file << nodes_vector[i] << "        " << m_pcs_out->GetNodeValue(nodes_vector[i], 1) << "\n";

			cout << "Data output: " << convertProcessTypeToString(getProcessType()) << " primary variables - SURFACE " << geo_name << " -  " << nodes_vector.size() << " nodes" << endl;
			break;
		case GEOLIB::POLYLINE:

				m_polyline = GEOGetPLYByName(geo_name);
				if (ply) {
					double min_edge_length(m_msh->getMinEdgeLength());
					m_msh->setMinEdgeLength(m_polyline->epsilon);
					m_msh->GetNODOnPLY(ply, nodes_vector);
				}

				for (long i = 0; i < (long)nodes_vector.size(); i++)
				if (convertProcessTypeToString(getProcessType()) == "MULTI_COMPONENTIAL_FLOW")  // PRESSURE1, TEMPERATURE1, CONCENTRATION1 
					tec_file << m_msh->nod_vector[i]->GetIndex() << "        " << m_pcs_out->GetNodeValue(m_msh->nod_vector[i]->GetIndex(), 1) << "        " << m_pcs_out->GetNodeValue(m_msh->nod_vector[i]->GetIndex(), 3) << "        " << m_pcs_out->GetNodeValue(m_msh->nod_vector[i]->GetIndex(), 5) << "\n";
				else
					tec_file << nodes_vector[i] << "        " << m_pcs_out->GetNodeValue(nodes_vector[i], 1) << "\n";

			cout << "Data output: " << convertProcessTypeToString(getProcessType()) << " primary variables - POLYLINE " << geo_name << " - " << nodes_vector.size() << " nodes" << endl;
			break;
		default:
			break;

			} // end case
			//////////////
			tec_file << "#STOP";
			tec_file.close();
		} // end output true
	//} // end timestep != 0
	
}

/**************************************************************************
FEMLib-Method:
Task:   Write water balance for polyline with leakance
Use:
Programing:
06/2012 JOD Implementation
**************************************************************************/

void COutput::WriteTotalFlux(double time_current, int time_step_number)
{
	CRFProcess* m_pcs = NULL;

	if (_nod_value_vector.size() > 0)
	  m_pcs = PCSGet(_nod_value_vector[0], true);
	else
	  m_pcs = PCSGet(getProcessType());

	double total_normal_flux_diff = 0, total_normal_flux_adv = 0;
	bool output = false;

	// File handling
	char number_char[3];
	string number_string = number_char;
	string tec_file_name = file_base_name + "_" + convertProcessTypeToString(getProcessType()) + "_";
	if (_nod_value_vector.size() > 0)
		tec_file_name += _nod_value_vector[0] + "_"; 

	if (this->getGeoType() == GEOLIB::POLYLINE)
		tec_file_name += "ply_";
	else if (this->getGeoType() == GEOLIB::SURFACE)
		tec_file_name += "surf_";
	else
		cout << "WARNING - geotype not supported in TOTAL_FLUX calculation";

	tec_file_name += geo_name + "_TOTAL_FLUX";

#if defined(USE_PETSC)  // JOD 2015-11-18
	tec_file_name += "_" + mrank_str;
#endif

	tec_file_name += ".txt";

	if(!_new_file_opened)
		remove(tec_file_name.c_str());

	fstream tec_file(tec_file_name.data(), ios::app | ios::out);
	tec_file.setf(ios::scientific, ios::floatfield);
	tec_file.precision(12);
	if (!tec_file.good()) return;
	tec_file.seekg(0L, ios::beg);

	if(!_new_file_opened)
	{
		tec_file << "\"TIME\"                   ";
		if (m_pcs->getProcessType() == FiniteElement::HEAT_TRANSPORT || m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
			tec_file << "\"DIFFUSION / DISPERSION FLUX\"             \"ADVECTION FLUX\"";
		else
			tec_file << "\"DARCY FLUX\"";
		tec_file << "\n";
	}
	else
	{
		if (time_vector.size() == 0 && (nSteps > 0) && (time_step_number
			% nSteps == 0))
			output = true;

		for (size_t j = 0; j < time_vector.size(); j++)
		if ((fabs(time_current - time_vector[j])) < MKleinsteZahl) //WW MKleinsteZahl
			output = true;

		if (output == true)
		{
			AccumulateTotalFlux(m_pcs, &total_normal_flux_diff, &total_normal_flux_adv);
			tec_file << time_current << "    " << total_normal_flux_diff << "         ";
			if ((m_pcs->getProcessType() == FiniteElement::HEAT_TRANSPORT) || (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT))
				tec_file << total_normal_flux_adv;
			tec_file << "\n";

			cout << "Data output: " << convertProcessTypeToString(getProcessType());
			if (_nod_value_vector.size() == 1)
				cout << " " << _nod_value_vector[0];
			else if (_nod_value_vector.size() > 1)
				cout << "WARNING - more than one node value specified to calculate";
			cout  << " TOTAL_FLUX " << geo_name << endl;
		}
	}
	tec_file.close();
	_new_file_opened = true;

}

/**************************************************************************
FEMLib-Method:
Task:
Use:
Programing:
10/2014 JOD Implementation
**************************************************************************/

void COutput::SetTotalFluxNodes(std::vector<long>& nodes_vector)
{

	switch (this->getGeoType()) {
	case GEOLIB::POLYLINE:
		SetTotalFluxNodesPLY(nodes_vector);
		break;
	case GEOLIB::SURFACE:
		SetTotalFluxNodesSURF(nodes_vector);
		break;
	default:
		cout << "Warning: Water Balance does not support this geotype" << endl;
	}
}

/**************************************************************************
FEMLib-Method:
Task:
Use:
Programing:
10/2014 JOD Implementation
**************************************************************************/

void COutput::SetTotalFluxNodesPLY(std::vector<long>& nodes_vector)
{

	GEOLIB::Polyline const* const ply(
		dynamic_cast<GEOLIB::Polyline const* const > (this->getGeoObj()));

	if (ply) {
		CGLPolyline* m_polyline = GEOGetPLYByName(geo_name);
		double min_edge_length(m_msh->getMinEdgeLength()); // ?????
		m_msh->setMinEdgeLength(m_polyline->epsilon);
		m_msh->GetNODOnPLY(ply, nodes_vector, true);
	}

}

/**************************************************************************
FEMLib-Method:
Task:
Use:
Programing:
10/2014 JOD Implementation
**************************************************************************/

void COutput::SetTotalFluxNodesSURF(std::vector<long>& nodes_vector)
{

	Surface* m_sfc = NULL;
	m_sfc = GEOGetSFCByName(geo_name);

	if (m_sfc)
		m_msh->GetNODOnSFC(m_sfc, nodes_vector, true);


}

/**************************************************************************
OpenGeoSys - Funktion
Task:   Writes contents for specific mmp_index in file, should work for all processes
Use:    $DAT_TYPE
CONTENT mmp_index  (-1 for whole domain)
Programing:
12/2014 JOD Implementation
**************************************************************************/

void COutput::WriteContent(double time_current, int time_step_number)
{
	CFEMesh* m_msh = NULL;
	m_msh = FEMGet(convertProcessTypeToString(getProcessType()));
	CRFProcess* m_pcs = NULL;
	double factor = 1.;
	if (_nod_value_vector.size() > 0)
		m_pcs = PCSGet(_nod_value_vector[0], true);
	else
		m_pcs = PCSGet(getProcessType());
	bool output = false;

	//--------------------------------------------------------------------
	/*if (m_msh->isAxisymmetry() && !_ignore_axisymmetry)
		factor = 6.283185307; // 2 Pi
	else
		factor = 1.;*/
	// File handling
	string tec_file_name;
	
	tec_file_name = file_base_name + "_" + convertProcessTypeToString(getProcessType()) ;
	if (_nod_value_vector.size() > 0)
		tec_file_name += "_" + _nod_value_vector[0];

	if (mmp_index == -2)
	{
		stringstream ss; ss << "_VOLUME_" << domainIntegration_lowerThreshold << "_" << domainIntegration_upperThreshold;
		tec_file_name += ss.str();
	}
	else
	{
		tec_file_name += "_CONTENT";
	}

	if (mmp_index >= 0)
	{
		stringstream ss; ss << mmp_index;
		tec_file_name += ss.str();
	}


#if defined(USE_PETSC)  // JOD 2015-11-18
	tec_file_name += "_" + mrank_str;
#endif

	tec_file_name += ".txt";

	 if(!_new_file_opened)
		remove(tec_file_name.c_str());

	fstream tec_file(tec_file_name.data(), ios::app | ios::out);
	tec_file.setf(ios::scientific, ios::floatfield);
	tec_file.precision(12);
	if (!tec_file.good()) return;
	tec_file.seekg(0L, ios::beg);

	if(!_new_file_opened)
		tec_file << "\"TIME\"                   \"CONTENT\"" << "\n";
	else
	{
		if (time_vector.size() == 0 && (nSteps > 0) && (time_step_number % nSteps == 0))
		  output = true;

		for (size_t j = 0; j < time_vector.size(); j++)
		if ((fabs(time_current - time_vector[j])) < MKleinsteZahl)
			output = true;

		if (output == true)
		{
			tec_file << time_current << "    " << factor * m_pcs->AccumulateContent(mmp_index,
					domainIntegration_lowerThreshold,
					domainIntegration_upperThreshold,
					_nod_value_vector)
							<< "\n";
			cout << "Data output: " << convertProcessTypeToString(getProcessType());
			if (_nod_value_vector.size() == 1)
				cout << " " << _nod_value_vector[0]; 
			cout << " TOTAL_CONTENT " << mmp_index << endl;
		}
	}
	_new_file_opened = true;
	tec_file.close();

}


/**************************************************************************
OpenGeoSys - Funktion
Task:  For axisymmetric W�rmesonde example
provisional - needs to be put to GEO together with integration in ST
Programming:
2/2015 JOD Implementation
**************************************************************************/
/*
void COutput::InterpolatePoints2Nodes(std::vector<double>&nod_val_vector)
{

	GEOLIB::Polyline const* polyline(dynamic_cast<GEOLIB::Polyline const*>(getGeoObj()));

	std::vector<double> interpolation_points;
	std::vector<double> interpolation_values;

	if (polyline) {
		std::vector<double> nodes_as_interpol_points;
		m_msh->getPointsForInterpolationAlongPolyline(polyline, nodes_as_interpol_points);
		//st->InterpolatePolylineNodeValueVector(nodes_as_interpol_points, ply_nod_val_vector);
		CGLPolyline* m_polyline = GEOGetPLYByName(geo_name);

		for (long i = 0; i < (long)DistribedBC.size(); i++) {

			interpolation_points.push_back(polyline->getLength(i));
			if (fabs(DistribedBC[i]) < MKleinsteZahl)
				interpolation_values.push_back(1.0e-20);
			else
				interpolation_values.push_back(DistribedBC[i]);
		}

		MathLib::PiecewiseLinearInterpolation(interpolation_points, interpolation_values, nodes_as_interpol_points, nod_val_vector);

	}
}
*/

/**************************************************************************
OpenGeoSys - Funktion
Task:  Gives advection and diffusion fluxes at nodes - Used in COutput::AccumulateTotalFlux
Programming:
2/2015 JOD Implementation
7/2015 JOD axisymmetry
**************************************************************************/

void COutput::NODCalcFlux(CRFProcess* m_pcs, CElem *elem, CElem* face, int* nodesFace, int nfn, double *NodeVal, double *NodeVal_adv)
{

	double flux[3], flux_normal, factor;
	CNode* e_node;
	CRFProcess *m_pcs_flow = NULL;
	for (size_t i = 0; i < pcs_vector.size(); i++)
	if (isFlowProcess(pcs_vector[i]->getProcessType()))
		m_pcs_flow = pcs_vector[i];
	int ndx1 = 1;
	if (_nod_value_vector.size() == 1) // for MASS_TRANSPORT
		ndx1 = m_pcs->GetNodeValueIndex(_nod_value_vector[0]) + 1;

	for (int k = 0; k < nfn; k++) {

		e_node = elem->GetNode(nodesFace[k]);

		if (m_pcs->m_msh->isAxisymmetry())
		{
			factor = elem->GetFluxArea() * m_pcs->m_msh->nod_vector[e_node->GetIndex()]->getData()[0] * 6.283185307; // 2 Pi x (x is radius) 
		}
		else
			factor = elem->GetFluxArea();


		if (m_pcs_flow == NULL)
			flux[0] = flux[1] = flux[2] = 0;
		else
		{
			flux[0] = m_pcs_flow->GetNodeValue(e_node->GetIndex(), m_pcs_flow->GetNodeValueIndex("VELOCITY_X1")); // Darcy
			flux[1] = m_pcs_flow->GetNodeValue(e_node->GetIndex(), m_pcs_flow->GetNodeValueIndex("VELOCITY_Y1"));
			flux[2] = m_pcs_flow->GetNodeValue(e_node->GetIndex(), m_pcs_flow->GetNodeValueIndex("VELOCITY_Z1"));
		}
		flux_normal = PointProduction(flux, face->normal_vector);

		switch (m_pcs->getProcessType())
		{
		case FiniteElement::LIQUID_FLOW:
			NodeVal[k] = flux_normal * factor;
			break;

		case FiniteElement::MASS_TRANSPORT: case FiniteElement::HEAT_TRANSPORT:

			NodeVal_adv[k] = m_pcs->GetNodeValue(e_node->GetIndex(), ndx1) * flux_normal * factor;				// advection

			if (m_pcs->getProcessType() == FiniteElement::HEAT_TRANSPORT)
				NodeVal_adv[k] *= mfp_vector[0]->SpecificHeatCapacity(NULL,true) * mfp_vector[0]->Density(); //BW: 23.03.2020 please double check
			// diffusion
			flux[0] = m_pcs->GetNodeValue(e_node->GetIndex(), m_pcs->GetNodeValueIndex("VELOCITY_X1"));  // Fick / Fourrier
			flux[1] = m_pcs->GetNodeValue(e_node->GetIndex(), m_pcs->GetNodeValueIndex("VELOCITY_Y1"));
			flux[2] = m_pcs->GetNodeValue(e_node->GetIndex(), m_pcs->GetNodeValueIndex("VELOCITY_Z1"));

			NodeVal[k] = PointProduction(flux, face->normal_vector) * factor;  // Fick / Fourrier
			break;

		default:
			cout << "ERROR - " << m_pcs->getProcessType() << " not supportet in TOTAL FLUX calculation" << endl;
		}

	}

}

/**************************************************************************
OpenGeoSys - Funktion
Task:  loops over elements and prepares integration, checks mmp_index, all elements, if mmp_index == -1
flux from element Gauss points extrapolated to nodes (and later on interpolated to face Gauss points)
Programming:
12/2014 JOD Implementation
**************************************************************************/
void COutput::AccumulateTotalFlux(CRFProcess* m_pcs, double* normal_flux_diff, double* normal_flux_adv) //const
{

	int nfaces, nfn, nodesFace[8], count;
	double fac, nodesFVal[8], nodesFVal_adv[8];
	int j, k, Axisymm = 1;                               // ani-axisymmetry
	if (m_pcs->m_msh->isAxisymmetry())
		Axisymm = -1;                               // Axisymmetry is true
	CNode* e_node;
	CElem *elem = NULL, *e_nei = NULL, *face = new CElem(1);
	FiniteElement::CElement* element = new FiniteElement::CElement(Axisymm * m_pcs->m_msh->GetCoordinateFlag());
	CFiniteElementStd* fem = new CFiniteElementStd(m_pcs, m_pcs->m_msh->GetCoordinateFlag());
	vector<long> nodes_on_geo, elements_at_geo;
	set<long> set_nodes_on_geo;

	// ----- initialize --------------------------------------------------------------------
	SetTotalFluxNodes(nodes_on_geo);  //get  nodes on geo object (put in preprocess)
	std::vector<double> nod_val_vector;
	nod_val_vector.resize(nodes_on_geo.size());
	//if (dis_type == "LINEAR")  //  used in radial W�rmesonde example with axisymmetry keyword - INACTIVE RIGHT NOW
	//	InterpolatePoints2Nodes(nod_val_vector);

	face->SetFace();
	for (long i = 0; i < (long)nodes_on_geo.size(); i++)
		set_nodes_on_geo.insert(nodes_on_geo[i]);
	m_pcs->m_msh->GetConnectedElements(nodes_on_geo, elements_at_geo);

	ele_gp_flux.clear();
	const size_t mesh_ele_vector_size(m_pcs->m_msh->ele_vector.size());
	for (size_t i = 0; i < mesh_ele_vector_size; i++)
		ele_gp_flux.push_back(new ElementValue(m_pcs, m_pcs->m_msh->ele_vector[i]));
		

	if ((m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT) || (m_pcs->getProcessType() == FiniteElement::HEAT_TRANSPORT)) {
		m_pcs->CalIntegrationPointValue();    //  calculate FICK / FOURRIER flux
		m_pcs->Extropolation_GaussValue();    //  and extrapolate to node
	}

	// face integration
	for (long i = 0; i < (long)elements_at_geo.size(); i++) {

		elem = m_pcs->m_msh->ele_vector[elements_at_geo[i]];
		if (!elem->GetMark())
			continue;
		nfaces = elem->GetFacesNumber();
		elem->SetOrder(m_pcs->m_msh->getOrder());
		
		for (j = 0; j < nfaces; j++) {

			e_nei = elem->GetNeighbor(j);
			nfn = elem->GetElementFaceNodes(j, nodesFace);
			// is element face on surface? 1st check
			if (elem->selected < nfn)
				continue;
			//2nd check
			count = 0;
			for (k = 0; k < nfn; k++)
			{
				e_node = elem->GetNode(nodesFace[k]);
				if (set_nodes_on_geo.count(e_node->GetIndex()) > 0)
					count++;
			}
			if (count != nfn)
				continue;
			// --------
			fac = 1.0;
			if (elem->GetDimension() == e_nei->GetDimension())
				fac = 0.5;   // Not a surface face
			face->SetFace(elem, j);
			face->SetOrder(m_pcs->m_msh->getOrder());
			face->ComputeVolume();
			face->SetNormalVector();
			face->DirectNormalVector();
			element->ConfigElement(face, m_pcs->m_num->ele_gauss_points, true); // 2D fem	
			//element->setOrder(m_pcs->m_msh->getOrder() + 1);
			//face->ComputeVolume();    
			NODCalcFlux(m_pcs, elem, face, nodesFace, nfn, nodesFVal, nodesFVal_adv);
			element->CalculateFluxThroughFace(elements_at_geo[i], fac, nodesFVal, nodesFVal_adv, normal_flux_diff, normal_flux_adv);
		} // end j, faces
	} // end i, elements at surface

	delete element;
	delete face;
	delete fem;

	ElementValue* gp_ele = NULL;
	if (ele_gp_flux.size() > 0)  // release memory
	{
		for (j = 0; j < (int)ele_gp_flux.size(); j++)
		{
			gp_ele = ele_gp_flux[j];
			delete gp_ele;
			gp_ele = NULL;
		}
		ele_gp_flux.clear();
	}
}

/**************************************************************************
OpenGeoSys - Funktion
Task:  PVD output
Programming:
8/2015 JOD Introduce function
**************************************************************************/

void COutput::WritePVD(double time_current, int time_step_number, bool output_by_steps, size_t no_times)
{
	if (vtk == NULL)
		CreateVTKInstance(); //WW m_out->vtk = new CVTK();
	//CVTK* vtk = m_out->vtk;

	bool vtk_appended = false;
	if (dat_type_name.find("PVD_A") != string::npos)
		vtk_appended = true;

	stringstream stm;
	string pvd_vtk_file_name, pvd_vtk_file_path;

	switch (getGeoType())
	{
	case GEOLIB::GEODOMAIN: // domain data
		if (time_step_number == 0)
		{
			std::string pcs_type("");
			if (getProcessType() != FiniteElement::INVALID_PROCESS)
				pcs_type = FiniteElement::convertProcessTypeToString(
				getProcessType());
			vtk->InitializePVD(file_base_name,
				pcs_type,
				vtk_appended);
		}
		// Set VTU file name and path
		pvd_vtk_file_name = vtk->pvd_vtk_file_name_base;
		stm << time_step_number;
		pvd_vtk_file_name += stm.str() + ".vtu";
#ifdef _WIN32
    pvd_vtk_file_path = vtk->pvd_vtk_file_path_base + "\\" + pvd_vtk_file_name;
#else
    pvd_vtk_file_path = vtk->pvd_vtk_file_path_base + "/" + pvd_vtk_file_name;
#endif

		// Output
		if (output_by_steps)
		{
			vtk->WriteXMLUnstructuredGrid(pvd_vtk_file_path, this,
				time_step_number);
			VTK_Info dat;
			dat.timestep = getTime();
			dat.vtk_file = pvd_vtk_file_name;
			vtk->vec_dataset.push_back(dat);
			vtk->UpdatePVD(vtk->pvd_file_name, vtk->vec_dataset);
		}
		else
		{
			for (size_t j = 0; j < no_times; j++)
			if (time_current >= time_vector[j])
			{
				vtk->WriteXMLUnstructuredGrid(
					pvd_vtk_file_name,
					this,
					time_step_number);
				time_vector.erase(
					time_vector.begin()
					+ j);
				VTK_Info dat;
				dat.timestep = getTime();
				dat.vtk_file = pvd_vtk_file_name;
				vtk->vec_dataset.push_back(dat);
				vtk->UpdatePVD(vtk->pvd_file_name,
					vtk->vec_dataset);
				break;
			}
		}
		break;

	default:
		break;
	}

}

/**************************************************************************
OpenGeoSys - Funktion
Task:  VTK output
Programming:
8/2015 JOD Introduce function
**************************************************************************/

void COutput::WriteVTK(double time_current, int time_step_number, bool output_by_steps, size_t no_times)
{

	CFEMesh* m_msh = NULL;
	m_msh = getMesh();

	switch (getGeoType())
	{
	case GEOLIB::GEODOMAIN: // domain data
		if (output_by_steps)
		{

		
			//OK
			//m_out->WriteDataVTK(time_step_number);
			LegacyVtkInterface vtkOutput(m_msh,
				_nod_value_vector,
				_ele_value_vector,
				mmp_value_vector,
				msh_type_name,
				this);
#if defined(USE_PETSC)						
			vtkOutput.WriteDataVTKPETSC(
				time_step_number,
				_time,
				file_base_name);
#else
			vtkOutput.WriteDataVTK(time_step_number,
				_time,
				file_base_name);
#endif
			if (!_new_file_opened)
				//WW
				_new_file_opened = true;
		}
		else
		{
			for (size_t j = 0; j < no_times; j++)
			if (time_current >= time_vector[j])
			{
				//OK
				//m_out->WriteDataVTK(time_step_number);
				LegacyVtkInterface vtkOutput(
					m_msh,
					_nod_value_vector,
					_ele_value_vector,
					mmp_value_vector,
					msh_type_name,
					this);
#if defined(USE_PETSC)						
				vtkOutput.WriteDataVTKPETSC(
					time_step_number,
					_time,
					file_base_name);
				time_vector.erase(
					time_vector.begin()
					+ j);
#else
				vtkOutput.WriteDataVTK(
					time_step_number,
					_time,
					file_base_name);
				time_vector.erase(
					time_vector.begin()
					+ j);

#endif
				if (!_new_file_opened)
					//WW
					_new_file_opened = true;
				break;


			}
		}
		break;
	default: // not domain
		break;
	}

}

/**************************************************************************
OpenGeoSys - Funktion
Task: TECPLOT / BINARY / MATLAB output
   switch between geotypes
Programming:
8/2015 JOD Introduce function
**************************************************************************/

void COutput::WriteTEC(double time_current, int time_step_number, bool output_by_steps, size_t no_times)
{
	void (COutput::*outputFunction)(int) = NULL;

	switch (getGeoType())
	{
		case GEOLIB::GEODOMAIN: // domain data
			cout << "Data output: Domain" << endl;		
			outputFunction = &COutput::WriteTEC_DOMAIN;		
			break;
			//------------------------------------------------------------------
		case GEOLIB::POLYLINE: // profiles along polylines
			cout << "Data output: Polyline profile - " << getGeoName() << endl;
			if (getProcessDistributionType() == FiniteElement::AVERAGE)
			{
				outputFunction = &COutput::NODWritePLYAverageDataTEC;
			}
			else
				outputFunction = &COutput::WriteTEC_POLYLINE;
			break;
			//------------------------------------------------------------------
		case GEOLIB::POINT: // breakthrough curves in points
			cout << "Data output: Breakthrough curves - " << getGeoName() << endl;
			outputFunction = &COutput::NODWritePNTDataTEC;
			//NODWritePNTDataTEC(time_current, time_step_number);
			break;
			//------------------------------------------------------------------
		case GEOLIB::SURFACE: // profiles at surfaces
			cout << "Data output: Surface profile" << endl;
			//..............................................................
			if (getProcessDistributionType() == FiniteElement::AVERAGE)
			{
				outputFunction = &COutput::NODWriteSFCAverageDataTEC;
			}
			//..............................................................
			else
				outputFunction = &COutput::NODWriteSFCDataTEC;
			//..............................................................
			// ELE data
			if (getElementValueVector().size() > 0)
				ELEWriteSFC_TEC();
			//..............................................................
			break;
			//			case 'Y': // Layer
			//				cout << "Data output: Layer" << "\n";
			//				outputFunction = &COutput::NODWriteLAYDataTEC;	
			//				break;
			//------------------------------------------------------------------

		default:
			break;
	} // end switch
	///
	if (outputFunction != NULL) // not for point, averaged surface since output already written
		WritePotentially(time_current, time_step_number, output_by_steps, no_times, outputFunction );

	if (!_new_file_opened)
		//WW
		_new_file_opened = true;
}

/**************************************************************************
OpenGeoSys - Funktion
Task: TECPLOT domain output
    switch between output of node or element (parameter) values 
Programming:
8/2015 JOD Introduce function
**************************************************************************/

void COutput::WriteTEC_DOMAIN(int time_step_number)
{

#if defined (USE_PETSC) // || defined (other parallel solver lib). 12.2012 WW
	if (dat_type_name.compare("BINARY") == 0) // 08.2012. WW
	{
		NODDomainWriteBinary();
	}
	else
	{
#endif
		if (_pcon_value_vector.size() > 0)
			PCONWriteDOMDataTEC();  //MX
		else
		{
			if (tecplot_datapack_block == 0) // BW
			{
				NODWriteDOMDataTEC();
				ELEWriteDOMDataTEC();
			}
			else if (tecplot_datapack_block == 1)  // BW Only block datapack type will be write into tecplot output
				BLOCKWriteDOMDataTEC();
			else if (tecplot_datapack_block == 2)
			{// Only block datapack type  and nodal data will be write into tecplot output
				NODWriteDOMDataTEC();
				BLOCKWriteDOMDataTEC();
			}
			else if (tecplot_datapack_block == 3)
			{// block datapack type, nodal data and ele data will be write into tecplot output
				NODWriteDOMDataTEC();
				BLOCKWriteDOMDataTEC();
				ELEWriteDOMDataTEC();
			}
		}
#if defined (USE_PETSC) // || defined (other parallel solver lib). 12.2012 WW
	}
#endif


}

/**************************************************************************
OpenGeoSys - Funktion
Task: TECPLOT polyline output
   Individual function because of individual function for file handling (TIMVALUE_TEC())
   although I do not see what this is used for
Programming:
8/2015 JOD Introduce function
**************************************************************************/

void COutput::WriteTEC_POLYLINE(int time_step_number)
{
	double tim_value;
	tim_value = NODWritePLYDataTEC(time_step_number);
	if (tim_value > 0.0)
		//OK
		TIMValue_TEC(tim_value);
}


/**************************************************************************
OpenGeoSys - Funktion
Task: 
Programming:
**************************************************************************/

void COutput::setFileBaseName(const std::string& fn)
{
	file_base_name = pathJoin(defaultOutputPath, pathBasename(fn));
}


/**************************************************************************
OpenGeoSys - Funktion
Task:
    calls output function if output by steps or if current time is one of the selected times (is in output vector)
	initial time in time_vector causes additional output step
	values in time_vector must increase with each entry
Programming:
11/2015 JOD Introduce function
**************************************************************************/

void COutput::WritePotentially(double time_current, int time_step_number, bool output_by_steps, size_t no_times, void(COutput::*outputFunction)(int))
{
	
	if (output_by_steps)
	{
		(this->*(outputFunction))(time_step_number);
		if (time_vector.size() > 0) 
		{       // this block added in case initial time is given in output instance (smaller times not permitted)
			if ( fabs(time_current - time_vector[0])
					< MKleinsteZahl) //WW MKleinsteZahl
					time_vector.erase(time_vector.begin());
		}
	}
    else
    {
		for (size_t j = 0; j < no_times; j++)
		{
			if ((time_current > time_vector[j]) ||
				fabs(time_current - time_vector[j])
				< MKleinsteZahl) //WW MKleinsteZahl
			{
				(this->*(outputFunction))(j + 1);
				time_vector.erase(time_vector.begin() + j);
			}
			break;				
		}
    }		
}

/**************************************************************************
OpenGeoSys - Funktion
Task:
Use:    $DAT_TYPE

Programing:
06/2018 JOD Implementation
**************************************************************************/

void COutput::WriteWellDoubletControl(double time_current, int time_step_number)
{	// a_pcs???
	std::cout << "Data output: WDC\n";
	m_pcs = PCSGet(getProcessType());
	if(m_pcs == NULL)
	{
		std::cout << "Warning - PCS not known for WellDoubletControl output" << "\n";
		return;
	}
	if(m_pcs->ogs_WDC_vector.size() == 0)
		std::cout << "Warning - No WDC instance in output\n";

	for(long long unsigned i=0; i<m_pcs->ogs_WDC_vector.size(); ++i)
	{  			// long long unsigned for std::to_string
		// file name
		std::string tec_file_name = file_base_name + "_"
		 + std::string(convertProcessTypeToString(getProcessType()))
				 + "_WellDoublet_" + std::to_string(i) + TEC_FILE_EXTENSION;
		// open file
		std::fstream tec_file;
		if (aktueller_zeitschritt == 0)
			tec_file.open(tec_file_name.data(), ios::out);
		else
			tec_file.open(tec_file_name.data(), ios::out | ios::app);

		tec_file.setf(ios::scientific, ios::floatfield);
		tec_file.precision(12);
		if (!tec_file.good())
		{
			std::cout << "Warning - Could not open file for writing WDC data " << geo_name << "\n";
			return;
		}
		tec_file.seekg(0L, ios::beg);
		#ifdef SUPERCOMPUTER
			// kg44 buffer the output
			char mybuffer [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
			tec_file.rdbuf()->pubsetbuf(mybuffer,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
			//
		#endif
		//--------------------------------------------------------------------
		if (aktueller_zeitschritt == 0)
		{
			tec_file << "TITLE = \"Well doublet " <<  i << "\"\n";
			tec_file << "VARIABLES = \"Step\" \"Time\" \"Scheme\" \"Power adaption flag\" ";
			tec_file << "\"Power rate Q_H\" \"Flow rate Q_w\" ";
			tec_file << "\"Warm well T_1\" \"Cold well T_2\" \"Heat exchanger T_HE\"\n";
		}
		else
		{
			// write results
			const wdc::WellDoubletControl::result_t& result = m_pcs->ogs_WDC_vector[i]->get_WellDoubletControl()->get_result();
			const OGS_WDC::doublet_mesh_nodes_t& doublet_mesh_nodes = m_pcs->ogs_WDC_vector[i]->get_doublet_mesh_nodes();

			tec_file << aktueller_zeitschritt
				<< '\t' << time_current
				<< '\t' << m_pcs->ogs_WDC_vector[i]->get_WellDoubletControl()->scheme_ID()
				<< '\t' << result.storage_state   // 0: powerrate_to_adapt, 1: on_demand
				<< '\t' << result.Q_H
				<< '\t' << result.Q_W
				<< '\t' << m_pcs->GetNodeValue(doublet_mesh_nodes.well1_aquifer, 1)
				<< '\t' << m_pcs->GetNodeValue(doublet_mesh_nodes.well2_aquifer, 1)
				<< '\t' << m_pcs->ogs_WDC_vector[i]->get_extremum(m_pcs, 1, doublet_mesh_nodes.heatExchanger)
				<< '\n';

		}
		tec_file.close();

	}
}

/***************************************************************************
OpenGeoSys - Funktion
Task:
Use:    $DAT_TYPE

Programing:
08/2019 JOD Implementation
**************************************************************************/

void COutput::WriteContraflow(double time_current, int time_step_number)
{	// a_pcs???
	std::cout << "Data output: Contraflow\n";
	m_pcs = PCSGet(getProcessType());
	if(m_pcs == NULL)
	{
		std::cout << "Warning - PCS not known for Contraflow output" << "\n";
		return;
	}
	if(m_pcs->ogs_contraflow_vector.size() == 0)
		std::cout << "Warning - No Contraflow instance in output\n";

	for(long long unsigned ii=0; ii<m_pcs->ogs_contraflow_vector.size(); ++ii)
	{  			// long long unsigned for std::to_string
		std::string tec_file_name_tf = file_base_name + "_"
				 + std::string(convertProcessTypeToString(getProcessType()))
						 + "_Contraflow_" + std::to_string(ii) + TEC_FILE_EXTENSION;
		// open file
		std::fstream tec_file_tf;
		if (aktueller_zeitschritt == 0)
			tec_file_tf.open(tec_file_name_tf.data(), ios::out);
		else
			tec_file_tf.open(tec_file_name_tf.data(), ios::out | ios::app);

		tec_file_tf.setf(ios::scientific, ios::floatfield);
		tec_file_tf.precision(12);
		if (!tec_file_tf.good())
		{
			std::cout << "Warning - Could not open file for writing Contraflow data " << geo_name << "\n";
			return;
		}
		tec_file_tf.seekg(0L, ios::beg);
		#ifdef SUPERCOMPUTER
			// kg44 buffer the output
			char mybuffer [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
			tec_file_tf.rdbuf()->pubsetbuf(mybuffer,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
			//
		#endif

		//--------------------------------------------------------------------
		if (aktueller_zeitschritt == 0)
		{
			tec_file_tf << "TITLE = \"Contraflow instance " <<  ii << "\"\n";
			tec_file_tf << "VARIABLES = \"Time\" \"T_in\" \"T_out\" \"flux_1\" \"flux_2\" \n";
		}
		else
		{

			if(m_pcs->ogs_contraflow_vector[ii]->get_input_list().front().Q > 10e-10)
			{

				// write results
				std::vector<long> nodes_vec =  m_pcs->ogs_contraflow_vector[ii]->get_nodes_vec();
				stru3::DVec T_s = stru3::DVec(nodes_vec.size());
				for(int j=0; j < nodes_vec.size(); ++j)
				{
					T_s[j] = m_pcs->GetNodeValue(nodes_vec[j], 1);
				}

				const contra::Result& result = m_pcs->ogs_contraflow_vector[ii]->get_Contraflow()->get_result();
				std::vector<contra::SegmentData> segment_data_vec = m_pcs->ogs_contraflow_vector[ii]->get_segment_data_vec();

				double total_flux1 = 0., total_flux2 = 0;
				double z = 0;


				int j = 0, k = 0;
				double L_ele = 0.;
				int N = segment_data_vec[0].N;
				double R_0_Delta = result.resistances_vec[0].R_0_Delta;  // segment allocation not verified
				double R_1_Delta = result.resistances_vec[0].R_1_Delta;

				for(int i=0; i < nodes_vec.size(); ++i)
				{
					if(k == 1)
					{
						L_ele = segment_data_vec[j].L/(N);
						R_0_Delta = result.resistances_vec[j].R_0_Delta;
						R_1_Delta = result.resistances_vec[j].R_1_Delta;
					}
					if(k == N)
					{
						R_0_Delta /= 2;
						R_1_Delta /= 2;
						L_ele /=2;
						k = 0;
						++j;
						if(j < segment_data_vec.size())
							N = segment_data_vec[j].N;
					}

					if(k == 0 && j < segment_data_vec.size())
					{
						R_0_Delta += result.resistances_vec[j].R_0_Delta;
						R_1_Delta += result.resistances_vec[j].R_1_Delta;

						L_ele += segment_data_vec[j].L/(2*N);
					}

					double flux1 = result.T_in[i] / R_0_Delta;
					flux1 -=  T_s[i] / R_0_Delta;
					flux1 *= L_ele;

					double flux2 = result.T_out[i] / R_1_Delta;
					flux2 -=  T_s[i] / R_1_Delta;
					flux2 *= L_ele;

					total_flux1 += flux1;
					total_flux2 += flux2;

					if(j < segment_data_vec.size())
						z += segment_data_vec[j].L/(N);

					 ++k;


					//const OGS_WDC::doublet_mesh_nodes_t& doublet_mesh_nodes = m_pcs->ogs_WDC_vector[i]->get_doublet_mesh_nodes();

					/*tec_file << aktueller_zeitschritt
						<< '\t' << time_current
						<< '\t' << m_pcs->ogs_WDC_vector[i]->get_WellDoubletControl()->scheme_ID()
						<< '\t' << result.storage_state   // 0: powerrate_to_adapt, 1: on_demand
						<< '\t' << result.Q_H
						<< '\t' << result.Q_W
						<< '\t' << m_pcs->GetNodeValue(doublet_mesh_nodes.well1_aquifer, 1)
						<< '\t' << m_pcs->GetNodeValue(doublet_mesh_nodes.well2_aquifer, 1)
						<< '\t' << m_pcs->ogs_WDC_vector[i]->get_extremum(m_pcs, 1, doublet_mesh_nodes.heatExchanger)
						<< '\n';
					 */
				} // end nodes_vec

				tec_file_tf << time_current << "\t" << result.T_in[0] << "\t" << result.T_out[0]
								<< "\t" << total_flux1 << "\t" << total_flux2 << "\n";

			} // end Q > 1.e-10
			else
			{
				tec_file_tf << time_current << "\t0\t0\t0\t0\n";
			}
		} // end aktueller_zeitschritt != 0
		tec_file_tf.close();

	}  // end ogs_contraflow_vector
}

/***************************************************************************
OpenGeoSys - Funktion
Task:
Use:    $DAT_TYPE

Programing:
04/2020 JOD Implementation
**************************************************************************/

void COutput::WriteContraflowPolyline(double time_current, int time_step_number)
{	// a_pcs???
	std::cout << "Data output: Contraflow\n";
	m_pcs = PCSGet(getProcessType());
	if(m_pcs == NULL)
	{
		std::cout << "Warning - PCS not known for Contraflow output" << "\n";
		return;
	}
	if(m_pcs->ogs_contraflow_vector.size() == 0)
		std::cout << "Warning - No Contraflow instance in output\n";

	for(long long unsigned ii=0; ii<m_pcs->ogs_contraflow_vector.size(); ++ii)
	{  			// long long unsigned for std::to_string
		// file name polyname
		std::string tec_file_name = file_base_name + "_"
				 + std::string(convertProcessTypeToString(getProcessType()))
				 + "_Contraflow_ply_" + std::to_string(ii) + TEC_FILE_EXTENSION;
		// open file
		std::fstream tec_file;
		if (aktueller_zeitschritt == 0)
			tec_file.open(tec_file_name.data(), ios::out);
		else
			tec_file.open(tec_file_name.data(), ios::out | ios::app);

		tec_file.setf(ios::scientific, ios::floatfield);
		tec_file.precision(12);
		if (!tec_file.good())
		{
			std::cout << "Warning - Could not open file for writing Contraflow data " << geo_name << "\n";
			return;
		}
		tec_file.seekg(0L, ios::beg);
		#ifdef SUPERCOMPUTER
			// kg44 buffer the output
			char mybuffer [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
			tec_file.rdbuf()->pubsetbuf(mybuffer,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
			//
		#endif


		//--------------------------------------------------------------------
		if (aktueller_zeitschritt != 0)
		{

			if(m_pcs->ogs_contraflow_vector[ii]->get_input_list().front().Q > 10e-10)
			{

				tec_file << "ZONE T=\"TIME=" << time_current << "\"\n";

				// write results
				std::vector<long> nodes_vec =  m_pcs->ogs_contraflow_vector[ii]->get_nodes_vec();
				stru3::DVec T_s = stru3::DVec(nodes_vec.size());
				for(int j=0; j < nodes_vec.size(); ++j)
				{
					T_s[j] = m_pcs->GetNodeValue(nodes_vec[j], 1);
				}

				const contra::Result& result = m_pcs->ogs_contraflow_vector[ii]->get_Contraflow()->get_result();
				std::vector<contra::SegmentData> segment_data_vec = m_pcs->ogs_contraflow_vector[ii]->get_segment_data_vec();

				double total_flux1 = 0., total_flux2 = 0;
				double z = 0;


				int j = 0, k = 0;
				double L_ele = 0.;
				int N = segment_data_vec[0].N;
				double R_0_Delta = result.resistances_vec[0].R_0_Delta;  // segment allocation not verified
				double R_1_Delta = result.resistances_vec[0].R_1_Delta;

				for(int i=0; i < nodes_vec.size(); ++i)
				{

					if(k == 1)
					{
						L_ele = segment_data_vec[j].L/(N);
						R_0_Delta = result.resistances_vec[j].R_0_Delta;
						R_1_Delta = result.resistances_vec[j].R_1_Delta;
					}
					if(k == N)
					{
						R_0_Delta /= 2;
						R_1_Delta /= 2;
						L_ele /=2;
						k = 0;
						++j;
						if(j < segment_data_vec.size())
							N = segment_data_vec[j].N;
					}

					if(k == 0 && j < segment_data_vec.size())
					{
						R_0_Delta += result.resistances_vec[j].R_0_Delta;
						R_1_Delta += result.resistances_vec[j].R_1_Delta;

						L_ele += segment_data_vec[j].L/(2*N);
					}

					double flux1 = result.T_in[i] / R_0_Delta;
					flux1 -=  T_s[i] / R_0_Delta;
					flux1 *= L_ele;

					double flux2 = result.T_out[i] / R_1_Delta;
					flux2 -=  T_s[i] / R_1_Delta;
					flux2 *= L_ele;

					total_flux1 += flux1;
					total_flux2 += flux2;
					tec_file << z << "\t" << T_s[i] << "\t" << result.T_in[i] << "\t" << result.T_out[i]
											<< "\t" << flux1 << "\t" << flux2
							//<< "\t" << R_0_Delta << "\t" << R_1_Delta
							<< '\n';

					if(j < segment_data_vec.size())
						z += segment_data_vec[j].L/(N);

					 ++k;


					//const OGS_WDC::doublet_mesh_nodes_t& doublet_mesh_nodes = m_pcs->ogs_WDC_vector[i]->get_doublet_mesh_nodes();

					/*tec_file << aktueller_zeitschritt
						<< '\t' << time_current
						<< '\t' << m_pcs->ogs_WDC_vector[i]->get_WellDoubletControl()->scheme_ID()
						<< '\t' << result.storage_state   // 0: powerrate_to_adapt, 1: on_demand
						<< '\t' << result.Q_H
						<< '\t' << result.Q_W
						<< '\t' << m_pcs->GetNodeValue(doublet_mesh_nodes.well1_aquifer, 1)
						<< '\t' << m_pcs->GetNodeValue(doublet_mesh_nodes.well2_aquifer, 1)
						<< '\t' << m_pcs->ogs_WDC_vector[i]->get_extremum(m_pcs, 1, doublet_mesh_nodes.heatExchanger)
						<< '\n';
					 */
				} // end nodes_vec


			} // end Q > 1.e-10

		} // end aktueller_zeitschritt != 0
		tec_file.close();

	}  // end ogs_contraflow_vector
}
