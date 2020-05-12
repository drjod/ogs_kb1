/**************************************************************************
   FEMLib - Object: OUT
   Task:
   Programing:
   06/2004 OK Implementation
   last modified:
**************************************************************************/
#include "Configure.h"
#include "makros.h"
// C++ STL
#include <cfloat>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <list>
#include <string>
using namespace std;

// Base
#include "StringTools.h"
// Math
#include "matrix_class.h" //JOD 2014-11-10
// FEM-Makros
#include "files0.h"
#include "makros.h"
// GeoSys-GeoLib
#include "GEOObjects.h"
#include "files0.h"
#include "geo_ply.h"
#include "geo_sfc.h"
// GeoSys-FEMLib
#include "LegacyVtkInterface.h"
#include "Output.h"
#include "fem_ele_std.h"
#include "mathlib.h"
#include "rf_msp_new.h"
#include "rf_out_new.h"
#include "rf_pcs.h"
#include "rf_pcs.h"
#include "rf_random_walk.h"
#include "rf_tim_new.h"
// GeoSys-MSHLib
#include "msh_lib.h"

// FileIO/FEMIO
#include "FEMIO/GeoIO.h"

#include "problem.h"

// Base
#include "StringTools.h"
#include "FileTools.h"

extern size_t max_dim;                            //OK411 todo


#ifdef CHEMAPP
#include "eqlink.h"
#endif
#include "vtk.h"
// MPI Parallel
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL)
#include "par_ddc.h"
#endif

#if defined(USE_PETSC) ||  defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL)//|| defined(other parallel libs)//03.3012. WW
#include <mpi.h>
#endif

#ifdef SUPERCOMPUTER
// kg44 this is usefull for io-buffering as endl flushes the buffer
#define endl '\n'
#define MY_IO_BUFSIZE 4096
#endif
#ifdef GEM_REACT
#include "rf_REACT_GEM.h"
#endif
using MeshLib::CFEMesh;
//==========================================================================
vector<COutput*>out_vector;


std::string defaultOutputPath = ""; // CL


/**************************************************************************
   FEMLib-Method:
   Task: OUT read function
   Programing:
   06/2004 OK Implementation
   08/2004 WW Remove the old files
   01/2005 OK Boolean type
   01/2005 OK Destruct before read
   06/2006 WW Remove the old files by new way
   06/2010 TF reformated, restructured, signature changed, use new GEOLIB data structures
**************************************************************************/
bool OUTRead(const std::string& file_base_name,
             const GEOLIB::GEOObjects& geo_obj,
             const std::string& unique_name)
{
	char line[MAX_ZEILE];
	std::string line_string;
	ios::pos_type position;
	bool output_version = false; // 02.2011. WW

#if defined(USE_PETSC) || defined(USE_MPI) //|| defined(other parallel libs)//03.3012. WW
	int rank , msize;
	string rank_str;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &msize);
	std::ifstream is;
	stringstream ss (stringstream::in | stringstream::out);
	ss.clear(); 
	ss.str("");
	ss << rank;
	rank_str = ss.str();
	ss.clear();
#endif
	// File handling
	std::string out_file_name = file_base_name + OUT_FILE_EXTENSION;
	std::ifstream out_file(out_file_name.data(), ios::in);
	if (!out_file.good())
		return false;
	out_file.seekg(0L, ios::beg);


	// Keyword loop
	cout << "OUTRead" << "\n";
	while (!out_file.eof())
	{
		out_file.getline(line, MAX_ZEILE);
		line_string = line;
		if (line_string.find("#STOP") != string::npos)
			return true;

		// Give version in file name
		//15.01.2008. WW
		if (line_string.find("#VERSION") != string::npos)
			output_version = true;  // 02.2011. WW
		//----------------------------------------------------------------------
		// keyword found
		if (line_string.find("#OUTPUT") != string::npos)
		{
			COutput* out(new COutput(out_vector.size()));

#if defined(USE_PETSC) || defined(USE_MPI) //|| defined(other parallel libs)//03.3012. WW
			out->setMPI_Info(rank, msize, rank_str);
#endif
			out->setFileBaseName(file_base_name);

			position = out->Read(out_file, geo_obj, unique_name);

			if(output_version) //// 02.2011. WW
			{
				std::string VersionStr = OGS_VERSION; //02.2011 WX
				int curPos = 0;
				int pos = 0;
				while((pos = VersionStr.find("/",curPos)) != -1)
				{
					VersionStr.replace(pos, 1, "_");
					curPos = pos + 1;
				}
				out->setFileBaseName(out->getFileBaseName() + "(V" + VersionStr + ")");
			}

			out_vector.push_back(out);

			//			char number_char[3]; //OK4709
			//			sprintf(number_char, "%i", (int) out_vector.size() - 1); //OK4709
			//			out->ID = number_char; //OK4709
			//			out->setID (out_vector.size() - 1);

			out_file.seekg(position, ios::beg);
		}                         // keyword found
	}                                     // eof
	return true;
}

/**************************************************************************
   FEMLib-Method: OUTWrite
   Task: master write function
   Programing:
   06/2004 OK Implementation
   last modification:
**************************************************************************/
void OUTWrite(string base_file_name)
{
	//========================================================================
	// File handling
	string out_file_name = base_file_name + OUT_FILE_EXTENSION;
	fstream out_file (out_file_name.data(),ios::trunc | ios::out);
	out_file.setf(ios::scientific,ios::floatfield);
	out_file.precision(12);
	if (!out_file.good())
		return;
	out_file.seekg(0L,ios::beg);
#ifdef SUPERCOMPUTER
	// kg44 buffer the output
	char mybuffer [MY_IO_BUFSIZE * MY_IO_BUFSIZE];
	out_file.rdbuf()->pubsetbuf(mybuffer,MY_IO_BUFSIZE * MY_IO_BUFSIZE);
	//
#endif
	//========================================================================
	out_file << "GeoSys-OUT: Output ------------------------------------------------\n";
	//========================================================================
	// OUT vector
	size_t out_vector_size (out_vector.size());
	for(size_t i = 0; i < out_vector_size; i++)
		out_vector[i]->Write(&out_file);
	out_file << "#STOP";
	out_file.close();
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   08/2004 OK Implementation
   08/2004 WW Output by steps given in .out file
   03/2005 OK RFO format
   05/2005 OK MSH
   05/2005 OK Profiles at surfaces
   12/2005 OK VAR,MSH,PCS concept
   03/2006 WW Flag to remove existing files
   08/2006 OK FLX calculations
   08/2007 WW Output initial values of variables
**************************************************************************/
void OUTData(double time_current, int time_step_number, bool force_output)
{
#if defined(USE_MPI) // JT2012
	if(myrank != 0) return;
#endif
	//
	COutput* m_out = NULL;
	//CRFProcess* m_pcs = NULL; 
	//CFEMesh* m_msh = NULL;
	bool OutputBySteps = false;
	//double tim_value;

	for (size_t i = 0; i < out_vector.size(); i++)
	{
	    OutputBySteps = false; // reset this flag for each COutput
		m_out = out_vector[i];
		// MSH               8/2016 removed by JOD since WARNING for common case
		//		m_msh = m_out->GetMSH();
		//m_msh = m_out->getMesh();
		//if (!m_msh)
		//	cout << "Warning in OUTData - no MSH data" << "\n";
		//OK continue;
		// PCS
		//if (m_out->_nod_value_vector.size() > 0)
		//	m_pcs = m_out->GetPCS(m_out->_nod_value_vector[0]);
		//if (m_out->getElementValueVector().size() > 0)
		//	m_pcs = m_out->GetPCS_ELE(m_out->getElementValueVector()[0]);
		//if (!m_pcs)
		//	m_pcs = m_out->GetPCS();  //OK
		//if (!m_pcs)
		//	cout << "Warning in OUTData - no PCS data" << "\n";


		//OK4704 continue;
		//--------------------------------------------------------------------'
		m_out->setTime (time_current);
		size_t no_times (m_out->time_vector.size());
		//--------------------------------------------------------------------
		if (no_times == 0 && (m_out->nSteps > 0) && (time_step_number
		                                             % m_out->nSteps == 0))
			OutputBySteps = true;
		if (time_step_number == 0 || force_output) //WW//JT
			OutputBySteps = true;
		//======================================================================
		// TECPLOT
		if (m_out->dat_type_name.compare("TECPLOT") == 0
		    || m_out->dat_type_name.compare("MATLAB") == 0 
		    || m_out->dat_type_name.compare("BINARY") == 0 // 08.2012. WW
           )
		   m_out->WriteTEC(time_current, time_step_number, OutputBySteps, no_times); // 8/2015 JOD 
		
		//--------------------------------------------------------------------
		// vtk
		else if (m_out->dat_type_name.compare("VTK") == 0)
		{
			m_out->WriteVTK(time_current, time_step_number, OutputBySteps, no_times);
		}                           // PVD (ParaView)
		else if (m_out->dat_type_name.find("PVD") != string::npos)
		{
			m_out->WritePVD(time_current, time_step_number, OutputBySteps, no_times);
		}
		else if (m_out->dat_type_name.compare("TOTAL_FLUX") == 0)
			m_out->WriteTotalFlux(time_current, time_step_number); // 6/2012 JOD, MW 
		else if (m_out->dat_type_name.compare("CONTENT") == 0)
			m_out->WriteContent(time_current, time_step_number); //JOD 2/2015
		else if (m_out->dat_type_name.compare("VOLUME") == 0)
			m_out->WriteContent(time_current, time_step_number); //JOD 2020-1-16
		else if (m_out->dat_type_name.compare("COMBINE_POINTS") == 0) 
			m_out->NODWritePointsCombined(time_current, time_step_number);	// 6/2012 for calibration JOD
		else if (m_out->dat_type_name.compare("PRIMARY_VARIABLES") == 0)
			m_out->NODWritePrimaryVariableList(time_current, time_step_number); //JOD 2014-11-10
		else if (m_out->dat_type_name.compare("WELL_DOUBLET_CONTROL") == 0)
					m_out->WriteWellDoubletControl(time_current, time_step_number); // 2018-06-27 JOD
		else if (m_out->dat_type_name.compare("CONTRAFLOW") == 0)
					m_out->WriteContraflow(time_current, time_step_number); // JOD 2019-08-23 JOD
		else if (m_out->dat_type_name.compare("CONTRAFLOW_POLYLINE") == 0)
					m_out->WriteContraflowPolyline(time_current, time_step_number); // JOD 2020-04-30 JOD
		// ELE values, only called if ele values are defined for output, 05/2012 BG
		if (m_out->getElementValueVector().size() > 0)
			m_out->CalcELEFluxes();

	}                                     // OUT loop
    
	//======================================================================
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2005 OK Implementation
   last modified:
**************************************************************************/
void OUTDelete()
{
	const size_t no_out = out_vector.size();
	for(size_t i = 0; i < no_out; i++)
		delete out_vector[i];
	out_vector.clear();
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2005 OK Implementation
   last modification: 03/2010 JT
   09/2010 TF
**************************************************************************/
COutput* OUTGet(const std::string & out_name)
{
	FiniteElement::ProcessType pcs_type (FiniteElement::convertProcessType (out_name));
	for (size_t i = 0; i < out_vector.size(); i++)
		if (out_vector[i]->getProcessType() == pcs_type)
			return out_vector[i];
	return NULL;
}

/**************************************************************************
   FEMLib-Method:
   Task: Return output object for variable name (based on OUTGet)
   Programing:
   03/2010 JT Implementation
   last modification:
**************************************************************************/
COutput* OUTGetRWPT(const std::string & out_name)
{
	for (size_t i = 0; i < out_vector.size(); i++)
	{
		COutput* out(out_vector[i]);
		for (size_t j = 0; j < out->getRandomWalkParticleTracingValueVector().size(); j++)
			if (out->getRandomWalkParticleTracingValueVector()[j].compare(out_name) ==
			    0)
				return out;
	}
	return NULL;
}

/***************************************************************************************
   Function checking the consistency of the output data as specified in the input file
   This means up to now, that data for missing processes is not written
   05/09	SB
 *****************************************************************************************/
void OUTCheck()
{
	std::cout << "Checking output data " << "\n";
	// Go through all out objects (#OUTPUT-section in input file)
	for (size_t i = 0; i < out_vector.size(); i++)
		out_vector[i]->checkConsistency();
}

