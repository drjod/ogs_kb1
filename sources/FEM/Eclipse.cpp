//1. Eclipse starten (die Files liegen alle vorbereitet vor
//2. Daten auslesen -> Tool ecl2gs
//3. Datenumwandeln (Uebertragung auf Knoten, Geschwindigkeit berechnen)
//4. Geosys daten schreiben
//5. Berechnung Geosys
//6. Eclipse restart

// declarations for ecl_interface ECLIPSE - GeoSys

//#include "stdafx.h"
#include "Eclipse.h"
#include <algorithm>
#include <cstdlib>
#include <fstream>  // Datei streams
#include <iostream>  // Bildschirmausgabe
#include <sstream>
#include <sstream>   // string streams (in)
#include <string>
#include <sys/types.h>
#include <vector>
//#include "dirent.h"
#include "display.h"
#include <errno.h>
#include <math.h>
//#include "Windows.h"
#include "mathlib.h"
#include <cfloat>
#include "rf_react.h" //SB redo by wtp
#include "rf_react_int.h"
#include <ctime>

#include "rf_mmp_new.h" // MAT redo by wtp
#include "rfmat_cp.h" // SB redo by wtp
#include <sys/stat.h> //for check if files exist

#include "fem_ele_std.h"
#include "fem_ele_vec.h"

using namespace std;

namespace process
{
	class CRFProcessDeformation;
}
using process::CRFProcessDeformation;
CECLIPSEBlock::CECLIPSEBlock(long Nodelength, long Facelength)
{
	(void)Facelength; // unused
	this->row = 0;
	this->column = 0;
	this->index = 0;
	this->layer = 0;
	this->Nodeindex.resize(Nodelength);//WB
	this->correspondingelenum = 0;//WB
	this->decomposetype = 0;
	this->x_coordinates.resize(Nodelength);
	this->y_coordinates.resize(Nodelength);
	this->z_coordinates.resize(Nodelength);
	this->x_barycentre = 0;
	this->y_barycentre = 0;
	this->z_barycentre = 0;
	this->active = 0;
	this->NeighbourElement.resize(6);
	this->FaceIndex.resize(6);
	this->volume = 0;
}

CECLIPSEBlock::~CECLIPSEBlock(void)
{
}

/*-------------------------------------------------------------------------
   GeoSys - Function: CReadTextfiles
   Task: Read textfiles in a vector of strings
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
CReadTextfiles_ECL::CReadTextfiles_ECL(void)
{
	this->NumberOfRows = 0;
	this->Data.clear();
}

CReadTextfiles_ECL::~CReadTextfiles_ECL(void)
{
	this->Data.clear();
}

bool CReadTextfiles_ECL::Read_Text(std::string Filename)
{
	char Line[MAX_ZEILEN];
	bool error = false;
	bool abort = false;
	std::string tempstring;

	// .data ... provides the filename as it is necessary for C, ios::in ... reads the file
	ifstream in_file(Filename.data(), ios::in);
	if (in_file.fail())
	{
		error = true;
		std::cout << "\n";
		std::cout << " ERROR: The file " << Filename.data() << " can not be opened!" << "\n";
	}
	else
	{
		error = false;
		while ((abort == false) && (in_file.eof() == false))
		{
			tempstring.clear();
			in_file.getline(Line, MAX_ZEILEN);
			tempstring = Line;
			// You could basically use this later to avoid line lenghts of pre-defined length only
			// getline(in_file,tempstring);
			if (tempstring.length() == MAX_ZEILEN - 1)
			{
				std::cout <<
					" ERROR: Increase MAX_ZEILEN in order to read ECLIPSE data file "
					<< "\n";
				std::cout << " or shorten the line in " << Filename.data() << ": " <<
					tempstring << " to " << MAX_ZEILEN << " characters " << "\n";
				exit(0);
			}
			if (tempstring.compare("#STOP") != 0)
			{
				this->Data.push_back(Line);
				this->NumberOfRows = this->NumberOfRows + 1;
			}
			else
				abort = true;
		}
	}
	in_file.close();
	return error;
}

CWriteTextfiles_ECL::CWriteTextfiles_ECL(void)
{
}

CWriteTextfiles_ECL::~CWriteTextfiles_ECL(void)
{
}

/*-------------------------------------------------------------------------
   GeoSys - Function: Write_Text
   Task: Writes a vector of strings to a text file
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
void CWriteTextfiles_ECL::Write_Text(std::string Filename, vector<std::string> Text)
{
	ofstream textfile;
	textfile.open(Filename.c_str(), ios::out);

	//WTP 02/2014 for (int i = 0; i < int(Text.size()); i++)
	for (int i = 0; i < static_cast<int>(Text.size()); i++)
	{
		//WTP 02/2014 if (i == int(Text.size()) - 1)
		if (i == static_cast<int>(Text.size()) - 1)
			textfile << Text[i].data();
		else
			textfile << Text[i].data() << "\n";
		//file.write(Text[i]);
	}
	textfile.close();
}

/*-------------------------------------------------------------------------
   Constructor and Destructor of the class CECLIPSEData
   -------------------------------------------------------------------------*/
CECLIPSEData::CECLIPSEData(void)
{
	rows = 0;
	columns = 0;
	elements = 0;
	layers = 0;
	//WTP 04/2014 this->times = 0;
	timestep_adjust_initial = 0;
	timestep_adjust_iteration_tot = 0;
	numberOutputParameters = 0;
	activeCells = 0;
	eclgridelenum = 0;
	NodeData.clear(); //resize(1); // warum das?
	Data = NULL;
	ProcessIndex_CO2inLiquid = -1;
	ProcessIndex_CO2inGas = -1; //redo wtp
	ProcessIndex_GasinLiquid = -1;
	Radial_I = false;
	Radial_J = false;
	RadialModellIpos = false;
	RadialModellJpos = false;
	phase_shift_flag = false; //redo wtp
	eclipse_ele_active_flag.clear(); // CB redo wtp
	PoroPermIncludeFile = false; //redo wtp
	//WTP 04/2014 this->dissolved_co2_pcs_name_ECL = ""; // SB_KB redo wtp
	//WTP 04/2014 this->dissolved_co2_ingas_pcs_name_ECL = ""; // SB_KB redo wtp
	TempIncludeFile = false;
	// WTP 04/2014 
	no_comps_ECL = 0;
	dm_pcs = NULL;
	option_CO2STORE = false;
	reservoir_conditions = true;
	Water_phase_exists = false;
	Oil_phase_exists = false;
	Gas_phase_exists = false;
	dissolved_comps = false;
	existWells = false;
	UsePrecalculatedFiles = false;
	UseSaveEclipseDataFiles = false;
	UseEclrun = false; // WTP 02/15
	surface_pressure = 0.0;
	surface_temperature = 0.0;
	pathECLFFile = "";
	pathECLFolder = "";
	pathECLProject = "";
	pathECLWellFile = "";
	nameECLProject = "";

    ReadEclipseGrid_called = false;
    ReadCorrespondingList_called = false;
    DetermineNeighbourElements_called = false;
    CompareElementsGeosysEclipse_called = false;
    CreateFaces_called = false;
    ReadWellData_called = false;
    ConnectFacesToElements_called = false;
	verbosity = 0;

	// WTP previously uninitialized valuesBGEIN
	IdentificationNumber = -1;
	FlowFaceNum = -1;
	Radialmodell = false;
	Molweight_CO2 = 44.009E-3; // [g/mol]
	Molweight_H2O = 18.0148E-3; // [g/mol]
	Molweight_NaCl = 58.443E-3; // [g/mol]
	SurfaceDensity_Gas_E100 = 1.9055; // [kg/m3] CO2
	E100 = false;
	actual_time = 0;;
	Windows_System = true;
}

CECLIPSEData::~CECLIPSEData(void)
{
}
/*-------------------------------------------------------------------------
   GeoSys - Function: GetVariableIndex
   Task: Returns the index of the given variable
   Return: Variableindex
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
int CECLIPSEData::GetVariableIndex(std::string Variablename)
{
	int index_Variable = -1;

	for (unsigned int i = 0; i < this->Variables.size(); i++)
	if (this->Variables[i].compare(Variablename) == 0)
		index_Variable = i;
	return index_Variable;
}

/*-------------------------------------------------------------------------
GeoSys - Function: GetECLFlowVariableName
Task: Returns the name of the flow variable
Return: string
Programming: 06/2015 WTP
Modification:
-------------------------------------------------------------------------*/
std::string CECLIPSEData::GetECLFlowVariableName(int phase, int direction)
{
	// phase flags: 1=water, 2=gas, 3=oil
	// direction flags: 1=I+, 2=J+, 3=K+
	switch (phase)
	{
	case 1:
		switch (direction)
		{
		case 1:
			if (reservoir_conditions)
				return "FLRWATI+";
			else
				return "FLOWATI+";

		case 2:
			if (reservoir_conditions)
				return "FLRWATJ+";
			else
				return "FLOWATJ+";

		case 3:
			if (reservoir_conditions)
				return "FLRWATK+";
			else
				return "FLOWATK+";
		}
	case 2:
		switch (direction)
		{
		case 1:
			if (reservoir_conditions)
				return "FLRGASI+";
			else
				return "FLOGASI+";

		case 2:
			if (reservoir_conditions)
				return "FLRGASJ+";
			else
				return "FLOGASJ+";

		case 3:
			if (reservoir_conditions)
				return "FLRGASK+";
			else
				return "FLOGASK+";
		}
	case 3:
		switch (direction)
		{
		case 1:
			if (reservoir_conditions)
				return "FLROILI+";
			else
				return "FLOOILI+";

		case 2:
			if (reservoir_conditions)
				return "FLROILJ+";
			else
				return "FLOOILJ+";

		case 3:
			if (reservoir_conditions)
				return "FLROILK+";
			else
				return "FLOOILK+";
		}
	}

	return "error";
}
/*-------------------------------------------------------------------------
   GeoSys - Function: SplitStrings
   Task: Separate a string with a delimiter
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
void CECLIPSEData::SplitStrings(const std::string str, std::string delimiter)
{
	this->SplittedString.clear();

	// Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiter, 0);
	// Find first "non-delimiter".
	string::size_type pos = str.find_first_of(delimiter, lastPos);

	while (string::npos != pos || string::npos != lastPos)
	{
		// Found a token, add it to the vector.
		this->SplittedString.push_back(str.substr(lastPos, pos - lastPos));
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiter, pos);
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiter, lastPos);
	}
}

/*-------------------------------------------------------------------------
   GeoSys - Function: Round
   Task: Round numbers
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
double CECLIPSEData::Round(double Number, int Decimalplaces)
{
	Number *= pow(10, double(Decimalplaces));
	Number = floor(Number + 0.5);
	Number *= pow(10, double(-Decimalplaces));
	return Number;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: Round
   Task: Round numbers
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
typeExponentialNumber CECLIPSEData::RoundEXP(double Number, int Decimalplaces)
{
	typeExponentialNumber Result;
	int Exponent;
	double tempNumber;
	//WW int sign;

	//WW sign = 1;
	//WW if (Number < 0)
	//WW	sign = -1;

	Exponent = 0;
	tempNumber = fabs(Number);
	if (tempNumber != 0)
	{
		if (tempNumber >= 1)
		{
			do {
				tempNumber /= 10;
				Exponent -= 1;
			} while (tempNumber >= 1);
		}
		else if (tempNumber < 0.1)
		{
			do {
				tempNumber *= 10;
				Exponent += 1;
			} while (tempNumber < 0.1);
		}
		Number *= pow(10, double(Exponent + Decimalplaces));
		Number = floor(Number + 0.5);
		Number *= pow(10, double(-Decimalplaces));
	}
	else
		Exponent = -0;
	Result.Exponent = -Exponent;
	Result.Number = Number;
	return Result;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: AddZero
   Task: Adds zero before or after the given number
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
std::string CECLIPSEData::AddZero(double Number, int Places, bool before)
{
	std::string tempstring;
	std::string Result;
	stringstream NumberString, tempNumberString;
	size_t position, position_e;
	double tempNumber;
	int precision;

	precision = Places + 5;
	NumberString.precision(precision);
	NumberString << Number;
	tempstring = NumberString.str();
	Result = tempstring;

	if (before == true)
		// WTP while (int(Result.length()) < Places)
	while (static_cast<int>(Result.length()) < Places)
	{
		if (Result.substr(0, 1) == "-")
			Result = "-0" + Result.substr(1, Result.length());
		else
			Result = "0" + Result;
	}
	else
	{
		position = Result.find(".");
		// if the string is longer than the string is cutted and rounded up
		if (Result.length() >= Places + position + 1)
		{
			tempNumber = this->Round(Number, Places);
			tempstring = "";
			precision = Places + position;
			tempNumberString.precision(precision);
			tempNumberString << tempNumber;
			tempstring = tempNumberString.str();
			Result = tempstring;
		}
		//else {
		while (size_t(Result.length()) < Places + position + 1)
		{
			position = Result.find(".");
			if (position < Result.length())
			{
				position_e = Result.find("e");
				if (position_e < Result.length())
					Result = Result.substr(0, position_e) + "0" + Result.substr(
					position_e,
					Result.length() - (position_e));
				else
					Result = Result + "0";
			}
			else
			{
				position_e = Result.find("e");
				if (position_e < Result.length())
					Result = Result.substr(0, position_e) + "." + Result.substr(
					position_e,
					Result.length() - (position_e));
				else
					Result = Result + ".";
			}
		}
		//}
	}
	return Result;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: SetFilenamesAndPaths
   Task: Sets all filenames and paths needed for the eclipse-ogs exectution
   Return: true if succeded
   Programming: 01/2014 WTP
   Modification:
   -------------------------------------------------------------------------*/
bool CECLIPSEData::SetFilenamesAndPaths(CRFProcess* m_pcs, long Timestep)
{
	bool success = true;
	int position;
	std::string root_folder;
	std::string geosys_folder;

	//check if eclipse folder is subfolder of the geosys folder
	if (Timestep == 1 && m_pcs->iter_outer_cpl == 0)
	{
		// Get timestep adjustment variable
		if (m_pcs->ecl_time_adjust > 0)
			this->timestep_adjust_initial = m_pcs->ecl_time_adjust;
		if (timestep_adjust_initial > 0)
			std::cout << " Attention: Eclipse simulation is a restart-based simulation" << "\n";

		if (m_pcs->simulator_path.find("eclrun") != std::string::npos)
			this->UseEclrun = true;
		// set the comeplete file name once
		this->pathECLProject = m_pcs->simulator_model_path;

		position = static_cast<int>(pathECLProject.find_last_of("\\"));
		if (position >= 0)
			this->Windows_System = true;
		else
			this->Windows_System = false;

		if (this->Windows_System == true)
			position = static_cast<int>(pathECLProject.find_last_of("\\"));
		else
			position = static_cast<int>(pathECLProject.find_last_of("/"));
		root_folder = pathECLProject.substr(0, position);
		this->pathECLFolder = pathECLProject.substr(0, position + 1);
		nameECLProject = pathECLProject.substr(position + 1, pathECLProject.length() - pathECLFolder.length()-5);

		if (this->Windows_System == true)
			position = static_cast<int>(root_folder.find_last_of("\\"));
		else
			position = static_cast<int>(root_folder.find_last_of("/"));
		root_folder = pathECLProject.substr(0, position);

		if (this->Windows_System == true)
			position = static_cast<int>(m_pcs->file_name_base.find_last_of("\\"));
		else
			position = static_cast<int>(m_pcs->file_name_base.find_last_of("/"));
		geosys_folder = m_pcs->file_name_base.substr(0, position);

		if (root_folder != geosys_folder && verbosity > 1)
		{
			// WTP: Why is this a problem?
			std::cout <<
				" ATTENTION: The Eclipse simulator model path is not part of the GeoSys model path!!!"
				<< "\n";
			std::cout << root_folder << "\n";
			std::cout << geosys_folder << "\n";
		}

		// check if filename is given with or without extension
		std::string tempstring = pathECLProject.substr(pathECLProject.length() - 5, pathECLProject.length());

		if ((tempstring.compare(".data") == 0) || (tempstring.compare(".DATA") == 0))
			pathECLProject = pathECLProject.substr(0, pathECLProject.length() - 5);

		// Now set the well filename if applicable
		pathECLWellFile = m_pcs->simulator_well_path; // set the path to the eclipse input data
		// if a well file name is given, set the corresponding bool flag
		if (pathECLWellFile != "")
			this->existWells = true;
	}

	if (Timestep > 1 || m_pcs->iter_outer_cpl > 0)
	{
		std::string tempstring = pathECLProject.substr(pathECLProject.length() - 7, pathECLProject.length());

		if (tempstring.compare("RESTART") != 0)
		{
			pathECLProject += "_RESTART";
			nameECLProject += "_RESTART";
		}
	}

	//check if file exists
	if (CheckIfFileExists(pathECLProject + ".DATA") == false)
	{
		std::cout << " ERROR: The ECLIPSE data file could not be found! (" << pathECLProject + ".DATA" << ")" << "\n";
		success = false;
	}

	// Update the *.F* file name
	//set the filename to the current output file -> is there a diff between e100 and e300?
	if (Timestep>1 || m_pcs->iter_outer_cpl > 0)
	{
		if (this->UsePrecalculatedFiles != true)//BW:running benchmark the file cannot be read, there is no temporaryresults.F, is that right? WTP: Yes, correct
		{
			if (m_pcs->Iterative_Eclipse_coupling == true)
				this->pathECLFFile = this->pathECLFolder + "TEMPORARYRESULTS.F" + AddZero(this->timestep_adjust_iteration_tot - 1 + timestep_adjust_initial, 4, true);
			else
				this->pathECLFFile = this->pathECLFolder + "TEMPORARYRESULTS.F" + AddZero(m_pcs->Tim->step_current - 1 + timestep_adjust_initial, 4, true);
		}
		else
			this->pathECLFFile = pathECLProject + ".F" + AddZero(m_pcs->Tim->step_current + timestep_adjust_initial, 4, true);
	};
	return success;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: CheckIfFileExists
   Task: Check if the file exists
   Return: true if exist
   Programming: 11/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
bool CECLIPSEData::CheckIfFileExists(string strFilename)
{
	// code source: http://www.techbytes.ca/techbyte103.html (6.11.2009) no restriction for use

	struct stat stFileInfo;
	bool blnReturn;
	int intStat;

	// Attempt to get the file attributes
	intStat = stat(strFilename.c_str(), &stFileInfo);
	if (intStat == 0)
		// We were able to get the file attributes
		// so the file obviously exists.
		blnReturn = true;
	else
		// We were not able to get the file attributes.
		// This may mean that we don't have permission to
		// access the folder which contains this file. If you
		// need to do that level of checking, lookup the
		// return values of stat which will give you
		// more details on why stat failed.
		blnReturn = false;

	return blnReturn;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: ReadEclipseGrid
   Task: Reads the Eclipse model grid
   Return: nothing
   Programming: 09/2009 BG
   Modification: Add collapse cell locator 04/2014 BW
   -------------------------------------------------------------------------*/
void CECLIPSEData::ReadEclipseGrid(std::string Filename)
{
	std::string tempstring;
	CECLIPSEBlock* m_eclipseblock = NULL;
	CReadTextfiles_ECL* TextFile;
	int Corners_row = 0;
	std::stringstream in;
	bool error = false;
	vector<std::string> files = vector<std::string>();
	clock_t start, finish;
	double time;

	start = clock();

	if (verbosity > 2)
	std::cout << "        ReadEclipseGrid()";

    if (this->ReadEclipseGrid_called == true)
    {
		if (verbosity > 2)
			std::cout << " - Function already called once, skip this time" << "\n";    
    }
    else
    {

        TextFile = new CReadTextfiles_ECL;
        error = TextFile->Read_Text(Filename);

        if (error == true)
        {
            std::cout << " ERROR: The program is canceled since the grid file could not be read." << "\n";
            //system("Pause");
            exit(0);
        }

        //Read in the Grid
        for (long i = 0; i < TextFile->NumberOfRows - 1; i++)
        {
   //         //wtp debug
			//std::cout << i << "\n";

            tempstring = TextFile->Data[i].substr(2, 8);

            //Check if the grid is radial or not
            if (tempstring.compare("RADIAL  ") == 0)
            {
                if (TextFile->Data[i + 1].substr(2, 8) == "FALSE   ")
                    Radialmodell = false;
                else if (TextFile->Data[i + 1].substr(2, 8) == "TRUE    ")
                    Radialmodell = true;
                else
                {
                    std::cout <<
						" ERROR: The phrase below the keyword RADIAL does not contain TRUE or FALSE! The program is aborted!"
                        << "\n";
                    //system("Pause");
                    exit(0);
                }
            }

            //Get the Dimesion of this ECLIPSE Grid (I,J,K) -> (Columns,Rows,Layers)
            if (tempstring.compare("DIMENS  ") == 0)
            {
                in.str((string)TextFile->Data[i + 1]);
                in >> this->columns >> this->rows >> this->layers;
                in.clear();

                // Initialize the index array as 3D vector
                this->IndexArray.resize(this->columns);
                for (long i = 0; i < this->columns; i++)
                {
                    this->IndexArray[i].resize(this->rows);
                    for (long j = 0; j < this->rows; ++j)
                        this->IndexArray[i][j].resize(this->layers);
                }
            }

            if (tempstring.compare("COORDS  ") == 0)
            {
                m_eclipseblock = new CECLIPSEBlock(8, 6); // new member of eclblock
                in.str((string)TextFile->Data[i + 1]);
                in >> m_eclipseblock->column >> m_eclipseblock->row >>
                    m_eclipseblock->layer >> m_eclipseblock->index >> m_eclipseblock->active;
                in.clear();

                // SB und KB
                //cout << m_eclipseblock->column << "\n";
                //this->SplitStrings(TextFile->Daten[i+1]," ");
                //if (this->SplittedString.size()<=4) {
                //	std::cout << "\n" << "The grid output format (GRIDFILE) for Eclipse should be set to 2 0 instead of 1 1! The run is terminated now!" << "\n";
                //	//system("Pause");
                //	exit(0);
                //}
                //m_eclipseblock->column = atoi(this->SplittedString[0].data());
                //m_eclipseblock->row = atoi(this->SplittedString[1].data());
                //m_eclipseblock->layer = atoi(this->SplittedString[2].data());
                //m_eclipseblock->index = atoi(this->SplittedString[3].data());
                //m_eclipseblock->active = atoi(this->SplittedString[4].data());
                //index of Eclipse start with 1 -> change to 0
                //m_eclipseblock->index = m_eclipseblock->index - 1;

                // inactive cells are not red anymore
                if (m_eclipseblock->active == 1)
                {
                    /*if (this->columns < m_eclipseblock->column)
                        this->columns = m_eclipseblock->column;
                        if (this->rows < m_eclipseblock->row)
                        this->rows = m_eclipseblock->row;
                        // layer number increases with depth -> layer 1 = top layer
                        if (this->layers < m_eclipseblock->layer)
                        this->layers = m_eclipseblock->layer;
                        */
                    if (m_eclipseblock->active == 1)
                        this->activeCells += 1;

                    m_eclipseblock->index = this->elements;
                    this->elements += 1;
                    //read in alle corner coordinates
                    //Debug.Print(Strings.Mid(Gridfile.Daten(0, i + 2), 3, 8))

                    if (TextFile->Data[i + 2].substr(2, 8).compare("CORNERS ") == 0)
                        Corners_row = 0;
                    else
                    {
                        if (TextFile->Data[i + 3].substr(2,
                            8).compare("CORNERS ") ==
                            0)
                            Corners_row = 1;
                        else
                            std::cout <<
                            " ERROR: Structure of the grid file not compatible with interface." <<
                            "\n";
                    }
                    //Reads the corner points; order of the points (x1<x2, y1<y2, z1<z2):
                    //	0..x1, y1, z2; 1..x2, y1, z2; 2..x1, y2, z2; 3..x2, y2, z2
                    //	4..x1, y1, z1; 5..x2, y1, z1; 6..x1, y2, z1; 7..x2, y2, z1
                    if (Radialmodell == false)  //Non Radial Model, Nodes are saved as the following order
                    {
                        in.str((string)TextFile->Data[i + 3 + Corners_row]);
                        in >> m_eclipseblock->x_coordinates[0] >>
                            m_eclipseblock->y_coordinates[0] >>
                            m_eclipseblock->z_coordinates[0] >>
                            m_eclipseblock->x_coordinates[1];
                        in.clear();

                        in.str((string)TextFile->Data[i + 4 + Corners_row]);
                        in >> m_eclipseblock->y_coordinates[1] >>
                            m_eclipseblock->z_coordinates[1] >>
                            m_eclipseblock->x_coordinates[2] >>
                            m_eclipseblock->y_coordinates[2];
                        in.clear();

                        in.str((string)TextFile->Data[i + 5 + Corners_row]);
                        in >> m_eclipseblock->z_coordinates[2] >>
                            m_eclipseblock->x_coordinates[3] >>
                            m_eclipseblock->y_coordinates[3] >>
                            m_eclipseblock->z_coordinates[3];
                        in.clear();

                        in.str((string)TextFile->Data[i + 6 + Corners_row]);
                        in >> m_eclipseblock->x_coordinates[4] >>
                            m_eclipseblock->y_coordinates[4] >>
                            m_eclipseblock->z_coordinates[4] >>
                            m_eclipseblock->x_coordinates[5];
                        in.clear();

                        in.str((string)TextFile->Data[i + 7 + Corners_row]);
                        in >> m_eclipseblock->y_coordinates[5] >>
                            m_eclipseblock->z_coordinates[5] >>
                            m_eclipseblock->x_coordinates[6] >>
                            m_eclipseblock->y_coordinates[6];
                        in.clear();

                        in.str((string)TextFile->Data[i + 8 + Corners_row]);
                        in >> m_eclipseblock->z_coordinates[6] >>
                            m_eclipseblock->x_coordinates[7] >>
                            m_eclipseblock->y_coordinates[7] >>
                            m_eclipseblock->z_coordinates[7];
                        in.clear();
                    }
                    else
                    {
                        /*----------------------------------------------------------------------------------------------------
                        Radial Model (8-nodes Block).The Numbering of 8-nodes is different with non-radial modal
                        In order to Keep the cosistence of following work, i.e., creat faces, the Nodes are saved as follows:
                        x0 y0 z0 x1;  ->   x2 y2 z2 x3;
                        y1 z1 x2 y2;  ->   y3 z3 x0 y0;
                        z2 x3 y3 z3;  ->   z0 x1 y1 z1;
                        x4 y4 z4 x5;  ->   x6 y6 z6 x7;
                        y5 z5 x6 y6;  ->   y7 z7 x4 y4;
                        z6 x7 y7 z7;  ->   z4 x5 y5 z5;
                        */

                        in.str((string)TextFile->Data[i + 3 + Corners_row]);
                        in >> m_eclipseblock->x_coordinates[2] >>
                            m_eclipseblock->y_coordinates[2] >>
                            m_eclipseblock->z_coordinates[2] >>
                            m_eclipseblock->x_coordinates[3];
                        in.clear();

                        in.str((string)TextFile->Data[i + 4 + Corners_row]);
                        in >> m_eclipseblock->y_coordinates[3] >>
                            m_eclipseblock->z_coordinates[3] >>
                            m_eclipseblock->x_coordinates[0] >>
                            m_eclipseblock->y_coordinates[0];
                        in.clear();

                        in.str((string)TextFile->Data[i + 5 + Corners_row]);
                        in >> m_eclipseblock->z_coordinates[0] >>
                            m_eclipseblock->x_coordinates[1] >>
                            m_eclipseblock->y_coordinates[1] >>
                            m_eclipseblock->z_coordinates[1];
                        in.clear();

                        in.str((string)TextFile->Data[i + 6 + Corners_row]);
                        in >> m_eclipseblock->x_coordinates[6] >>
                            m_eclipseblock->y_coordinates[6] >>
                            m_eclipseblock->z_coordinates[6] >>
                            m_eclipseblock->x_coordinates[7];
                        in.clear();

                        in.str((string)TextFile->Data[i + 7 + Corners_row]);
                        in >> m_eclipseblock->y_coordinates[7] >>
                            m_eclipseblock->z_coordinates[7] >>
                            m_eclipseblock->x_coordinates[4] >>
                            m_eclipseblock->y_coordinates[4];
                        in.clear();

                        in.str((string)TextFile->Data[i + 8 + Corners_row]);
                        in >> m_eclipseblock->z_coordinates[4] >>
                            m_eclipseblock->x_coordinates[5] >>
                            m_eclipseblock->y_coordinates[5] >>
                            m_eclipseblock->z_coordinates[5];
                        in.clear();
                    }

                    //m_eclipseblock->CalcBarycentre(); // For block to node interpolation
                    //m_eclipseblock->CalculateFaceCentres();

                    //std::cout << m_eclipseblock->x_coordinates[0] << " " << m_eclipseblock->x_coordinates[1] << " " << m_eclipseblock->x_coordinates[2] << "\n";
                    // adds 1 dataset of typ eclblocks to vec_eclblocks
                    this->IndexArray[m_eclipseblock->column - 1][m_eclipseblock->row - 1][m_eclipseblock->layer - 1]
                        = m_eclipseblock->index;

                    this->eclgrid.push_back(m_eclipseblock);
                    this->eclipse_ele_active_flag.push_back(true); //SB redo WTP
                }
                else
                {
                    this->IndexArray[m_eclipseblock->column - 1][m_eclipseblock->row - 1][m_eclipseblock->layer - 1]
                        = -1;
                    this->eclipse_ele_active_flag.push_back(false); //SB redo WTP
                }
            }
        }

        //recalculate coordinates if it is a radial grid based on angle and radius
        double alpha = 0;
        double radius = 0;
        if (Radialmodell == true)
        {
            for (long i = 0; i < this->elements; i++)
                //WTP for (int j = 0; j < int(this->eclgrid[i]->x_coordinates.size()); j++)
            for (int j = 0; j < static_cast<int>(this->eclgrid[i]->x_coordinates.size()); j++)
            {
                //coordinates are defined as radius, angle (clockwise in degree), heigth
                //recalculate alpha from degree to radian (Bogenma�)
                //std::cout << "Element: " << i << " Point: " << j << " alpha: " << this->eclgrid[i]->y_coordinates[j] << " radius: " << this->eclgrid[i]->x_coordinates[j];
                alpha = this->eclgrid[i]->y_coordinates[j] * PI / 180;
                radius = this->eclgrid[i]->x_coordinates[j];
                this->eclgrid[i]->x_coordinates[j] = sin(alpha) * radius;
                this->eclgrid[i]->y_coordinates[j] = cos(alpha) * radius;
                //std::cout << " radius: " << radius << " x: " << this->eclgrid[i]->x_coordinates[j] << " y: " << this->eclgrid[i]->y_coordinates[j] << "\n";
            }
            //this->eclgrid[i]->CalcBarycentre();
            double sum = 0;
            // check if radial model is in positiv I (EAST) direction  (KB)
			for (long i = 0; i < this->elements; i++)
            {
                //WTP for (int k = 0; k < int(this->eclgrid[i]->y_coordinates.size()); k++)
                for (int k = 0; k < static_cast<int>(this->eclgrid[i]->y_coordinates.size()); k++)
                    sum = sum + this->eclgrid[i]->y_coordinates[k];
                //for (int j = 0; j < int(this->eclgrid[i]->x_coordinates.size()); j++)
                for (int j = 0; j < static_cast<int>(this->eclgrid[i]->x_coordinates.size()); j++)
                {
                    if ((static_cast<int>(this->eclgrid[i]->x_coordinates[j]) >= 0) &&
                        (sum > -0.00000001) && (sum < 0.0000001)) // SB had them at *10^-1, redo wtp
                        this->RadialModellIpos = true;
                    else
                    {
                        std::cout <<
							" ERROR: The radial model is not perpendicular to the positive I axis (East direction)."
                            << "\n";
                        std::cout.flush();
                        exit(0);
                    }
                }
            }
        }

        //Release memory
        delete (TextFile);
        // set bool flag
        this->ReadEclipseGrid_called = true;

        finish = clock();
        time = (static_cast<double>(finish)-static_cast<double>(start)) / CLOCKS_PER_SEC;
		if (verbosity > 2)
			std::cout << "                     Time: " << time << " seconds." << "\n";
    }
}
/*-------------------------------------------------------------------------
   GeoSys - Function: ReadCorresponding info List
   Task: Reads the Eclipse Corresonding Elements info and Node Index Info
   Return: nothing
   Programming: 05/2014 BW
   -------------------------------------------------------------------------*/
void CECLIPSEData::ReadCorrespondingList(std::string Filename)
{
	std::string tempstring;

	CReadTextfiles_ECL* TextFile;

	std::stringstream in;
	bool error = false;
	vector<std::string> files = vector<std::string>();
	clock_t start, finish;
	double time;

	//TICK
	start = clock();
	
	if (verbosity > 2)
		std::cout << "        ReadCorrespondingList()";

    if (this->ReadCorrespondingList_called == true)
    {
		if (verbosity > 2)
			std::cout << " - Function already called once, skip this time" << "\n";
    }
    else
    {

        TextFile = new CReadTextfiles_ECL;
        error = TextFile->Read_Text(Filename);

        if (error == true)
        {
            std::cout << " ERROR:: Corresponding *.list file does not exist." << "\n";
            //system("Pause");
            exit(0);
        }

        //Read in the list
        int localindex = 0;
        int correspondingelenum = 0;
        //Read in the IndentificationNumber
        in.str((string)TextFile->Data[1]);
        in >> this->IdentificationNumber;
        in.clear();

        for (long i = 0; i < TextFile->NumberOfRows - 2; i++)
        {
            in.str((string)TextFile->Data[i + 2]);
            in >> localindex >> correspondingelenum;

            //Index needs to be identical
            if (localindex != this->eclgrid[i]->index)
            {
                std::cout << " ERROR: The index of ECLIPSE block is not identical." << "\n";
                //system("Pause");
                exit(0);
            }
            else{

                //Read Corresponding List According to the Index of ECLIPSE Block
                this->eclgrid[i]->correspondingelenum = correspondingelenum;
                this->eclgrid[i]->correspondingeleindex.resize(correspondingelenum);

                for (int j = 0; j < correspondingelenum; j++){
                    in >> this->eclgrid[i]->correspondingeleindex[j];
                }

                //Read Node Index into ECLIPSE Block Data
                in >> this->eclgrid[i]->Nodeindex[0]
                    >> this->eclgrid[i]->Nodeindex[1]
                    >> this->eclgrid[i]->Nodeindex[3]
                    >> this->eclgrid[i]->Nodeindex[2]
                    >> this->eclgrid[i]->Nodeindex[4]
                    >> this->eclgrid[i]->Nodeindex[5]
                    >> this->eclgrid[i]->Nodeindex[7]
                    >> this->eclgrid[i]->Nodeindex[6];
                in.clear();
            }
        }

        //Release memory
        delete (TextFile);

        //Kick out Double index in Nodexindex save as MshNodeindex
        for (size_t i = 0; i < this->eclgrid.size(); i++){
            this->eclgrid[i]->MshNodeindex.push_back(this->eclgrid[i]->Nodeindex[0]);
            for (size_t j = 1; j < this->eclgrid[i]->Nodeindex.size(); j++) {
                bool kick = false;
                for (size_t k = 0; k < this->eclgrid[i]->MshNodeindex.size(); k++){
                    if (this->eclgrid[i]->MshNodeindex[k] == this->eclgrid[i]->Nodeindex[j])
                    {
                        kick = true;
                        break;
                    }
                }
                if (kick == false)
                    this->eclgrid[i]->MshNodeindex.push_back(this->eclgrid[i]->Nodeindex[j]);
            }
        }
        this->ReadCorrespondingList_called = true;
        // TOCK
        finish = clock();
        time = (static_cast<double>(finish)-static_cast<double>(start)) / CLOCKS_PER_SEC;
		if (verbosity > 2)
	        std::cout << "               Time: " << time << " seconds." << "\n";
    }
}

/*-------------------------------------------------------------------------
   GeoSys - Function: DetermineNeighbourElements
   Task: Checks the consitency of both grids
   Return: true or false
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
void CECLIPSEData::DetermineNeighbourElements(string Projectname)
{
	string Filename;
	string tempstring;
	CReadTextfiles_ECL* TextFile;
	std::stringstream in;
	//double Element; double X; double Y; double Z; double N1; double N2; double N3; double N4; double N5; double N6;
	double X, Y, Z;
	long Element, N[6];
	bool error = false;
	bool neighbours_exist = true;
	bool generate_neighbours = false;
	clock_t start, finish;
	double time;
	// WTP double fractpart, param, intpart;
	start = clock();

	if (verbosity > 2)
		std::cout << "        DetermineNeighbourElements()";
    if(this->DetermineNeighbourElements_called == true)
    {
		if (verbosity > 2)
			std::cout << " - Function already called once, skip this time" << "\n";
    }
    else
    {
        //check if file *.neighbours exists, which stores the neighbours of each element
        Filename = Projectname + ".neighbours";
		neighbours_exist = CheckIfFileExists(Filename);
		if (neighbours_exist)
		{
			TextFile = new CReadTextfiles_ECL;
			error = TextFile->Read_Text(Filename);
			long tempindentification;
			size_t tempblocknumber;
			//	bool Error_NeighbourFile = false;
			//	double epsilon = 1e-3;
			if (error == true)
				generate_neighbours = true;
			else
			{
				// Check Indentification Number, if false, give warning information
				in.str((string)TextFile->Data[1]);
				in >> tempindentification;
				in.clear();
				if (this->IdentificationNumber != tempindentification && verbosity > 1)
				{
					std::cout << " WARNING: '*.neighbours' file is not generated at the same time as the '*.list' and '*.bfaces' file!" << "\n";
					std::cout << "          This might cause wrong information on the geometry of ECLIPSE grid!" << "\n";
					//system("Pause");
					//exit(0);
				}

				// Check if the total number of blocks are identical, if false, give ERROR information and stop the program
				in.str((string)TextFile->Data[2]);
				in >> tempblocknumber;
				in.clear();
				if (this->eclgrid.size() != tempblocknumber)
				{
					std::cout << " ERROR: The number of ECLIPSE blocks in the '*.neighbours' file is not identical to the grid read!" << "\n";
					std::cout << flush;
					//system("Pause");
					exit(0);
				}

				//Read the neighbours of each element
				for (long i = 3; i < TextFile->NumberOfRows; i++)
				{
					in.str((string)TextFile->Data[i]);
					in >> Element >> X >> Y >> Z;
					for (int j = 0; j < 6; j++)
						in >> N[j];
					in.clear();

					if (Element == i - 3)
					{
						for (int j = 0; j < 6; j++)
							this->eclgrid[i]->NeighbourElement[j] = N[j];
						this->eclgrid[Element]->x_barycentre = X;
						this->eclgrid[Element]->y_barycentre = Y;
						this->eclgrid[Element]->z_barycentre = Z;
					}
					else
					{
						std::cout << " ERROR: Block " << Element << " index is not identical to the index in the *.neighbours file!" << "\n";
						std::cout << "        Please re-run the OGS-ECL mesh converter!" << "\n";
						std::cout << flush;
						//system("Pause");
						exit(0);
					}
				}
				delete (TextFile);
			}
		}
		else
			generate_neighbours = true;


        if (generate_neighbours == true)
        {
            //std::cout << "        Creating new *.neighbours file";

            //Initailize the NeighbourElements grid as -1*6
            for (size_t i = 0; i < this->eclgrid.size(); i++)
            {
                for (size_t j = 0; j < this->eclgrid[i]->NeighbourElement.size(); j++)
                    this->eclgrid[i]->NeighbourElement[j] = -2;
            }

            //Order of the faces: 0..left(x); 1..right(x); 2..front(y); 3..back(y); 4..top(z); 5..bottom(z)
            //Loop over faces, determine neighbours according to the face block index
            for (long i = 0; i < this->FlowFaceNum; i++)
            {
                //determine according to face category
                switch (this->faces[i]->category)
                {
                case 1:	 //I+
                    if (this->faces[i]->connected_blocks.size() == 1)
                        this->eclgrid[this->faces[i]->connected_blocks[0]]->NeighbourElement[1] = -1;
                    if (this->faces[i]->connected_blocks.size() == 2)
                    {
                        this->eclgrid[this->faces[i]->connected_blocks[0]]->NeighbourElement[1] = this->faces[i]->connected_blocks[1];
                        this->eclgrid[this->faces[i]->connected_blocks[1]]->NeighbourElement[0] = this->faces[i]->connected_blocks[0];
                    }
                    break;
                case 2:	//I-
                    if (this->faces[i]->connected_blocks.size() == 1)
                        this->eclgrid[this->faces[i]->connected_blocks[0]]->NeighbourElement[0] = -1;
                    else
                    {
                        cout << " ERROR: I-face " << i << " has more than one block attached to it." << '\n';
                        std::cout << flush;
                        exit(0);
                    }
                    break;
                case 3:	//J+				
                    if (this->faces[i]->connected_blocks.size() == 1)
                        this->eclgrid[this->faces[i]->connected_blocks[0]]->NeighbourElement[3] = -1;
                    if (this->faces[i]->connected_blocks.size() == 2)
                    {
                        this->eclgrid[this->faces[i]->connected_blocks[0]]->NeighbourElement[3] = this->faces[i]->connected_blocks[1];
                        this->eclgrid[this->faces[i]->connected_blocks[1]]->NeighbourElement[2] = this->faces[i]->connected_blocks[0];
                    }
                    break;
                case 4:	 //J-
                    if (this->faces[i]->connected_blocks.size() == 1)
                        this->eclgrid[this->faces[i]->connected_blocks[0]]->NeighbourElement[2] = -1;
                    else
                    {
                        cout << " ERROR: I-face " << i << " has more than one block attached to it." << '\n';
                        std::cout << flush;
                        exit(0);
                    }
                    break;
                case 5:	 //K+				
                    if (this->faces[i]->connected_blocks.size() == 1)
                        this->eclgrid[this->faces[i]->connected_blocks[0]]->NeighbourElement[5] = -1;
                    if (this->faces[i]->connected_blocks.size() == 2)
                    {
                        this->eclgrid[this->faces[i]->connected_blocks[0]]->NeighbourElement[5] = this->faces[i]->connected_blocks[1];
                        this->eclgrid[this->faces[i]->connected_blocks[1]]->NeighbourElement[4] = this->faces[i]->connected_blocks[0];
                    }
                    break;
                case 6:	 //K-
                    if (this->faces[i]->connected_blocks.size() == 1)
                        this->eclgrid[this->faces[i]->connected_blocks[0]]->NeighbourElement[4] = -1;
                    else
                    {
						cout << " ERROR: I-face " << i << " has more than one block attached to it." << '\n';
                        std::cout << flush;
                        exit(0);
                    }
                    break;
                default:
                    cout << " ERROR: Unknown face category found!" << '\n';
                    std::cout << flush;
                    exit(0);
                }
            }


            //Order of the faces: 0..left(x); 1..right(x); 2..front(y); 3..back(y); 4..top(z); 5..bottom(z)
            for (size_t i = 0; i < this->eclgrid.size(); i++)
            {

                //Left Neighbour
				if (this->eclgrid[i]->NeighbourElement[0] == -2)
				{
					if (this->eclgrid[i]->FaceIndex[1] == -1) //I- Face is non-exist, Neighbour is -1
						this->eclgrid[i]->NeighbourElement[0] = -1;
					else if (verbosity > 1)
						std::cout << " WARNING: Neighbour " << 0 << " of block " << i << " is not determined." << '\n';
				}

                //Right Neighbour
				if (this->eclgrid[i]->NeighbourElement[1] == -2)
				{
					if (this->eclgrid[i]->FaceIndex[0] == -1) //I+ Face is non-exist, Neighbour is -1
						this->eclgrid[i]->NeighbourElement[1] = -1;
					else if (verbosity > 1)
						std::cout << " WARNING: Neighbour " << 1 << " of block " << i << " is not determined." << '\n';
				}

                //Front Neighbour	
				if (this->eclgrid[i]->NeighbourElement[2] == -2)
				{
					if (this->eclgrid[i]->FaceIndex[3] == -1) //I- Face is non-exist, Neighbour is -1
						this->eclgrid[i]->NeighbourElement[2] = -1;
					else if (verbosity > 1)
						std::cout << " WARNING: Neighbour " << 2 << " of block " << i << " is not determined." << '\n';
				}

                //Back Neighbour
				if (this->eclgrid[i]->NeighbourElement[3] == -2)
				{
					if (this->eclgrid[i]->FaceIndex[2] == -1) //I- Face is non-exist, Neighbour is -1
						this->eclgrid[i]->NeighbourElement[3] = -1;
					else if (verbosity > 1)
						std::cout << " WARNING: Neighbour " << 3 << " of block " << i << " is not determined." << '\n';
				}

                //Up Neighbour
				if (this->eclgrid[i]->NeighbourElement[4] == -2)
				{
					if (this->eclgrid[i]->FaceIndex[5] == -1) //I- Face is non-exist, Neighbour is -1
						this->eclgrid[i]->NeighbourElement[4] = -1;
					else if (verbosity > 1)
						std::cout << " WARNING: Neighbour " << 4 << " of block " << i << " is not determined." << '\n';
				}

                //Down Neighbour
                if (this->eclgrid[i]->NeighbourElement[5] == -2)
				{
					if (this->eclgrid[i]->FaceIndex[4] == -1) //I- Face is non-exist, Neighbour is -1
	                    this->eclgrid[i]->NeighbourElement[5] = -1;
					else if (verbosity > 1)
	                    std::cout << " WARNING: Neighbour " << 5 << " of block " << i << " is not determined." << '\n';
				}
            }

            ////data output
            //vector <string> vec_string;
            //ostringstream temp;
            //tempstring = "Element  X  Y  Z  N1  N2  N3  N4  N5  N6";
            //vec_string.push_back(tempstring);
            //// Loop over all elements
            //for (long i = 0; i < this->elements; i++)
            //{
            //	temp.precision(12);
            //	temp.str("");
            //	temp.clear();
            //	temp << this->eclgrid[i]->index;
            //	tempstring = temp.str();
            //	temp.str("");
            //	temp.clear();
            //	temp << this->eclgrid[i]->x_barycentre;
            //	tempstring += "  " + temp.str();
            //	temp.str("");
            //	temp.clear();
            //	temp << this->eclgrid[i]->y_barycentre;
            //	tempstring += "  " + temp.str();
            //	temp.str("");
            //	temp.clear();
            //	temp << this->eclgrid[i]->z_barycentre;
            //	tempstring += "  " + temp.str();
            //	for (int j = 0; j < 6; j++)
            //	{
            //		temp.str("");
            //		temp.clear();
            //		temp << this->eclgrid[i]->NeighbourElement[j];
            //		tempstring += "  " + temp.str();
            //	}
            //	vec_string.push_back(tempstring);
            //} // end element loop

            //// Test Output
            //Filename = Projectname + ".neighbours";
            //ofstream aus;
            //aus.open(Filename.data(), ios::out);
            //for (unsigned int i = 0; i < vec_string.size(); i++)
            //	aus << vec_string[i] << "\n";
            //aus.close();
        }

        DetermineNeighbourElements_called = true;

        //TOCK;
        finish = clock();
        time = (static_cast<double>(finish)-static_cast<double>(start)) / CLOCKS_PER_SEC;
		if (verbosity > 2)
			std::cout << "          Time: " << time << " seconds." << "\n";
    }
}

/*-------------------------------------------------------------------------
   GeoSys - Function: AnalyzeDataFromInputFile
   Task: Read the necessary data from .data- File: density of CO2 at surface conditions (for E100),
   necessary for calculating co2 concentration from RS value, Well Data
   Return: true or false
   Programming: 01/2011 BG
   Modification: 02/2011 KB Change from ReadSurfaceDensity to ReadDataFromInputFile
   08/2013 WTP Generalization & Adaptation of the interface (added E300 Comp support)
   04/2014 WTP Changed name from ReadDataFromInputFile to AnalyzeDataFromInputFile
   -------------------------------------------------------------------------*/
int CECLIPSEData::AnalyzeDataFromInputFile(CReadTextfiles_ECL* eclDataFile, CRFProcess* m_pcs)
{
	std::string tempstring;
	//CBoundaryConditions* m_BoundaryConditions = NULL; // unused
	//CReadTextfiles_ECL *TextFile;
	std::stringstream in;
	double density = -1.0;
	double pressure_s = -1.0;
	double temperature_s = -1.0;
	//double mw = -1.0;
	const double epsilon = 1e-3;
	//    bool error;
	//bool save_schedule_section = false;

	//WW bool success = false;
	std::string dummy_rate, dummy_zeile, rest_zeile, name, phase, open_flag, control_mode;
	//int jj = 0;
	std::string incl("GRID_PROPS.INC"); //WTP Maybe we should deal with this with a seperate keyword
	size_t found;    //SB redo wtp

	double time;
	clock_t start, finish;
	start = clock();
	//TICK;
	if (verbosity > 2)
		std::cout << "        AnalyzeDataFromInputFile()";

	if (reservoir_conditions == false)
	{
		std::cout << "\n";
		std::cout << " -----------------------------------------------------------------------------" << "\n";
		std::cout << " WARNING: Flow units are set to surface conditions for the ogs-ecl coupling!" << "\n";
		std::cout << "          Be advised that this functionality should only be used with extreme" << "\n";
		std::cout << "          care as the results obtained are not correct for pretty much every" << "\n";
		std::cout << "          subsurface system! In order to get correct results delete the keyword " << "\n";
		std::cout << "          $USE_ECL_SURFACE_UNITS in the *.pcs file." << "\n";
		std::cout << " -----------------------------------------------------------------------------" << "\n";
		if (verbosity > 2)
			std::cout << "                                  ";
	};

	//Read in the input file and search for various keyword
 	for (long i = 0; i < eclDataFile->NumberOfRows; i++)
	{
		this->SplittedString.clear();
		this->SplitStrings(eclDataFile->Data[i], " ");
		if (this->SplittedString.size() > 0)
		{
			long zeilen = 0;
			if (this->SplittedString[0] == "WATER")
			{
				Water_phase_exists = true;
				this->E100 = true;
			}
			if (this->SplittedString[0] == "GAS")
			{
				Gas_phase_exists = true;
				this->E100 = true;
			}
			if (this->SplittedString[0] == "OIL")
			{
				Oil_phase_exists = true;
				this->E100 = true;
			}
			if (this->SplittedString[0] == "CO2STORE")
			{
				this->Water_phase_exists = true;        // in the CO2STORE case water and gas are always included
				this->Gas_phase_exists = true;
				option_CO2STORE = true;
				this->E100 = false;
			}
			// wtp: added functionality to read a variable composition
			if (this->SplittedString[0] == "COMPS")
			{
				this->E100 = false;
				no_comps_ECL = atoi(eclDataFile->Data[i + 1].data());
			}
			if (this->SplittedString[0] == "INLCUDE")
			{
				do
				{
					zeilen = zeilen + 1;
				} while (eclDataFile->Data[i + zeilen] == ""); // look for possible empty rows between the Keyword and the items
				while ((eclDataFile->Data[i + zeilen] != "") && (eclDataFile->Data[i + zeilen] != "/"))
				{
					if (eclDataFile->Data[i + zeilen].substr(0, 2) != "--")  // ignore comments
					{
						found = eclDataFile->Data[i + zeilen].find(incl);
						if (found != string::npos)
						{
							PoroPermIncludeFile = true;
							break;
						}
					}
					zeilen = zeilen + 1;
				}
			}
			if (this->SplittedString[0] == "CNAMES")
			{
				do
				{
					zeilen = zeilen + 1;
				} while (eclDataFile->Data[i + zeilen] == ""); // look for possible empty rows between the Keyword and the items
				while ((eclDataFile->Data[i + zeilen] != "") && (eclDataFile->Data[i + zeilen] != "/"))
				{
					string delimiter = " ";
					SplittedString.clear();
					SplitStrings(eclDataFile->Data[i + zeilen], delimiter);
					if (eclDataFile->Data[i + zeilen].substr(0, 2) != "--") // ignore comments
					{
						for (int j = 0; j < no_comps_ECL; j++)
						{
							std::string dummy_str = "";
							dummy_str = SplittedString[j].data();
							vec_cnames_ecl.push_back(dummy_str);
						}
					}
					zeilen = zeilen + 1;
				}
				// checking with the already read pcs names from the *.pcs file
				for (unsigned int j = 0; j < m_pcs->vec_component_pcs_names.size(); j++)
				{
					std::string test_str_pcs = m_pcs->vec_component_pcs_names[j][0];
					if (vec_components_ECL_OGS[j][1] == 1 || vec_components_ECL_OGS[j][2] == 1 || vec_components_ECL_OGS[j][3] == 1)
					{
						bool found_comp = false;
						for (int k = 0; k < no_comps_ECL; k++)
						{
							std::string test_str_cnames = vec_cnames_ecl[k];
							if (test_str_pcs == test_str_cnames)
							{
								found_comp = true;
								vec_components_ECL_OGS[j][0] = k + 1;
							}
						}
						if (found_comp == false)
						{
							std::cout << " ERROR: Component \"" << m_pcs->vec_component_pcs_names[j][0] << "\" not found within the *.data file \n";
							exit(0);
						}
					}
				}
			}
			// reading the molweights from Eclipse
			if (this->SplittedString[0] == "MW")
			{
				do
				{
					zeilen = zeilen + 1;
				} while (eclDataFile->Data[i + zeilen] == ""); // look for possible empty rows between the Keyword and the items

				//zeilen = zeilen + 1;//End of the block is characterised with "" or "/"
				while ((eclDataFile->Data[i + zeilen] != "") && (eclDataFile->Data[i + zeilen] != "/"))
				{
					string delimiter = " ";
					SplittedString.clear();
					SplitStrings(eclDataFile->Data[i + zeilen], delimiter);
					if (eclDataFile->Data[i + zeilen].substr(0, 2) != "--") // ignore comments
					{
						for (int j = 0; j < no_comps_ECL; j++)
						{
							double dummy_dbl = 0.;
							//dummy_dbl = atof(SplittedString[j].data());
							dummy_dbl = atof(SplittedString[j].data());
							vec_ecl_mw_comps.push_back(dummy_dbl * 1.0E-3);
						}
					}
					zeilen = zeilen + 1;
				}
				for (unsigned int j = 0; j < vec_cnames_ecl.size(); j++)
				{
					for (unsigned int k = 0; k < this->vec_components_ECL_OGS_pcs_names.size(); k++)
					{
						/*string tempstring_wasser = vec_components_ECL_OGS_pcs_names[k][1];
						string tempstring_gas = vec_components_ECL_OGS_pcs_names[k][2];*/

						if (vec_cnames_ecl[j] == vec_components_ECL_OGS_pcs_names[k][0])
						{
							for (unsigned int l = 0; l<cp_vec.size(); l++)
							{
								std::string dummy_str = cp_vec[l]->compname;
								if (dummy_str == vec_components_ECL_OGS_pcs_names[k][1] || dummy_str == vec_components_ECL_OGS_pcs_names[k][2] || dummy_str == vec_components_ECL_OGS_pcs_names[k][3])
								{
									// If the component is used in both OGS and ECL the molweight specified in each simulator 
									// must be checked for uniformity
									if ((fabs(cp_vec[l]->molar_weight / vec_ecl_mw_comps[j])) > 1. + epsilon)
									{
										std::cout << " ERROR: The Molar weight of the component " << vec_components_ECL_OGS_pcs_names[k][0]
											<< " differs between OGS and Eclipse! " << "\n";
										std::cout << " PCS name: " << cp_vec[l]->compname << "\n";
										std::cout << "     Molweight OGS:     " << cp_vec[l]->molar_weight << "\n";
										std::cout << "     Molweight ECL:     " << vec_ecl_mw_comps[j] << "\n";
										std::cout << "     Difference:     " << fabs(cp_vec[l]->molar_weight - vec_ecl_mw_comps[j]) << "\n";
										std::cout << "     accepted epsilon:     " << epsilon << "\n";
                                        std::cout << flush;
										exit(0);
									}
								}
							}
						}
					}
				}
			}

			// reading the standard conditions specified in the eclipse.data
			if (this->SplittedString[0] == "STCOND")
			{
				do {
					zeilen = zeilen + 1;
				} while (eclDataFile->Data[i + zeilen] == ""); // look for possible empty rows between the Keyword and the items

				//zeilen = zeilen + 1;//End of the block is characterised with "" or "/"
				while ((eclDataFile->Data[i + zeilen] != "") &&
					(eclDataFile->Data[i + zeilen] != "/"))
				{
					string delimiter = " ";
					SplittedString.clear();
					SplitStrings(eclDataFile->Data[i + zeilen], delimiter);
					if (eclDataFile->Data[i + zeilen].substr(0, 2) != "--") // ignore comments
					{
						temperature_s = atof(SplittedString[0].data());
						temperature_s += 273.15;
						pressure_s = atof(SplittedString[1].data());
						pressure_s *= 100000.0;
						break;
					}
					zeilen = zeilen + 1;
				}
			}
			// reading the phase densities from the eclipse.data
			if (this->SplittedString[0] == "DENSITY")
			{
				do {
					zeilen = zeilen + 1;
				} while (eclDataFile->Data[i + zeilen] == ""); // look for possible empty rows between the Keyword and the items

				//zeilen = zeilen + 1;//End of the block is characterised with "" or "/"
				while ((eclDataFile->Data[i + zeilen] != "") &&
					(eclDataFile->Data[i + zeilen] != "/"))
				{
					string delimiter = " ";
					SplittedString.clear();
					SplitStrings(eclDataFile->Data[i + zeilen], delimiter);
					if (eclDataFile->Data[i + zeilen].substr(0, 2) != "--") // ignore comments
					{
						density = atof(SplittedString[2].data());
						// WTP This is dangerous since the (gas)density can, but does not have to be included in an E300.DATA file!
						// Thus a correct usage of Eclipse could lead to unwated program crashes
						this->SurfaceDensity_Gas_E100 = density;
						//check if the number of bc is plausibel
						if (density > 800 || density < 0)
						{
							if(this->E100 == true)    // WTP or just read the DENSITY keyword if it is E100?
							{
								std::cout <<
									" ERROR: The surface density of the gas seems not to be correct: "
									<< density << "\n";
							 std::cout << flush;
								return 0;
							}
						}
						else
							break;
					}
					zeilen = zeilen + 1;
				}
			}

			// if a well file is used read and save the whole SCHEDULE section of the input file
			//if (existWells == true && this->SplittedString[0] == "SCHEDULE")
			//    save_schedule_section = true;    // WTP: currently this is a rather weird approach...
			if (existWells == true)
			{
				if (this->SplittedString[0] == "SCHEDULE")
				{
					std::string dummy_str = eclDataFile->Data[i];
					this->ECL_schedule_section.push_back(dummy_str);
				}
				// save the messages keyword
				if (this->SplittedString[0] == "MESSAGES")
				{
					do
					zeilen = zeilen + 1;
					while (eclDataFile->Data[i + zeilen] == ""); // look for empty rows between the Keyword and the items
					while ((eclDataFile->Data[i + zeilen] != "") && (eclDataFile->Data[i + zeilen] != "/"))
					{
						string delimiter = " ";
						SplittedString.clear();
						SplitStrings(eclDataFile->Data[i + zeilen], delimiter);
						if (eclDataFile->Data[i + zeilen].substr(0, 2) != "--") // ignore comments
						{
							std::string loc_dummy_str = "";
							loc_dummy_str = eclDataFile->Data[i + zeilen];
							ECL_keyword_messages.push_back(loc_dummy_str);
						}
						zeilen = zeilen + 1;
					}
				}

				// save the rptrst keyword
				if (this->SplittedString[0] == "RPTRST")
				{
					do
					zeilen = zeilen + 1;
					while (eclDataFile->Data[i + zeilen] == ""); // look for empty rows between the Keyword and the items
					while ((eclDataFile->Data[i + zeilen] != "") && (eclDataFile->Data[i + zeilen] != "/"))
					{
						string delimiter = " ";
						SplittedString.clear();
						SplitStrings(eclDataFile->Data[i + zeilen], delimiter);
						if (eclDataFile->Data[i + zeilen].substr(0, 2) != "--") // ignore comments
						{
							std::string loc_dummy_str = "";
							loc_dummy_str = eclDataFile->Data[i + zeilen];
							ECL_keyword_rptrst.push_back(loc_dummy_str);
						}
						zeilen = zeilen + 1;
					}
				}
				// save the WELSPECS keyword
				if (this->SplittedString[0] == "WELSPECS")
				{
					do
					zeilen = zeilen + 1;
					while (eclDataFile->Data[i + zeilen] == ""); // look for empty rows between the Keyword and the items

					while ((eclDataFile->Data[i + zeilen] != "") && (eclDataFile->Data[i + zeilen] != "/"))
					{
						string delimiter = " ";
						SplittedString.clear();
						SplitStrings(eclDataFile->Data[i + zeilen], delimiter);
						if (eclDataFile->Data[i + zeilen].substr(0, 2) != "--") // ignore comments
						{
							std::string loc_dummy_str = "";
							loc_dummy_str = eclDataFile->Data[i + zeilen];
							ECL_keyword_welspecs.push_back(loc_dummy_str);
						}
						zeilen = zeilen + 1;
					}
				}
				// save the COMPDAT keyword
				if (this->SplittedString[0] == "COMPDAT")
				{
					do
					zeilen = zeilen + 1;
					while (eclDataFile->Data[i + zeilen] == ""); // look for empty rows between the Keyword and the items

					while ((eclDataFile->Data[i + zeilen] != "") && (eclDataFile->Data[i + zeilen] != "/"))
					{
						string delimiter = " ";
						SplittedString.clear();
						SplitStrings(eclDataFile->Data[i + zeilen], delimiter);
						if (eclDataFile->Data[i + zeilen].substr(0, 2) != "--") // ignore comments
						{
							std::string loc_dummy_str = "";
							loc_dummy_str = eclDataFile->Data[i + zeilen];
							ECL_keyword_compdat.push_back(loc_dummy_str);
						}
						zeilen = zeilen + 1;
					}
				}
				// save the WELLSTRE keyword
				if (this->SplittedString[0] == "WELLSTRE")
				{
					do
					zeilen = zeilen + 1;
					while (eclDataFile->Data[i + zeilen] == ""); // look for empty rows between the Keyword and the items

					while ((eclDataFile->Data[i + zeilen] != "") && (eclDataFile->Data[i + zeilen] != "/"))
					{
						string delimiter = " ";
						SplittedString.clear();
						SplitStrings(eclDataFile->Data[i + zeilen], delimiter);
						if (eclDataFile->Data[i + zeilen].substr(0, 2) != "--") // ignore comments
						{
							std::string loc_dummy_str = "";
							loc_dummy_str = eclDataFile->Data[i + zeilen];
							ECL_keyword_wellstre.push_back(loc_dummy_str);
						}
						zeilen = zeilen + 1;
					}
				}
				// save the winjgas keyword
				if (this->SplittedString[0] == "WINJGAS")
				{
					do
					zeilen = zeilen + 1;
					while (eclDataFile->Data[i + zeilen] == ""); // look for empty rows between the Keyword and the items

					while ((eclDataFile->Data[i + zeilen] != "") && (eclDataFile->Data[i + zeilen] != "/"))
					{
						string delimiter = " ";
						SplittedString.clear();
						SplitStrings(eclDataFile->Data[i + zeilen], delimiter);
						if (eclDataFile->Data[i + zeilen].substr(0, 2) != "--") // ignore comments
						{
							std::string loc_dummy_str = "";
							loc_dummy_str = eclDataFile->Data[i + zeilen];
							ECL_keyword_winjgas.push_back(loc_dummy_str);
						}
						zeilen = zeilen + 1;
					}
				}
				// save the wconinje keyword
				if (this->SplittedString[0] == "WCONINJE")
				{
					do
					zeilen = zeilen + 1;
					while (eclDataFile->Data[i + zeilen] == ""); // look for empty rows between the Keyword and the items

					while ((eclDataFile->Data[i + zeilen] != "") && (eclDataFile->Data[i + zeilen] != "/"))
					{
						string delimiter = " ";
						SplittedString.clear();
						SplitStrings(eclDataFile->Data[i + zeilen], delimiter);
						if (eclDataFile->Data[i + zeilen].substr(0, 2) != "--") // ignore comments
						{
							std::string loc_dummy_str = "";
							loc_dummy_str = eclDataFile->Data[i + zeilen];
							ECL_keyword_wconinje.push_back(loc_dummy_str);
						}
						zeilen = zeilen + 1;
					}
				}

				// save the zippy keyword
				if (this->SplittedString[0] == "ZIPPY2")
				{
					do
					zeilen = zeilen + 1;
					while (eclDataFile->Data[i + zeilen] == ""); // look for empty rows between the Keyword and the items

					while ((eclDataFile->Data[i + zeilen] != "") && (eclDataFile->Data[i + zeilen] != "/"))
					{
						string delimiter = " ";
						SplittedString.clear();
						SplitStrings(eclDataFile->Data[i + zeilen], delimiter);
						if (eclDataFile->Data[i + zeilen].substr(0, 2) != "--") // ignore comments
						{
							std::string loc_dummy_str = "";
							loc_dummy_str = eclDataFile->Data[i + zeilen];
							ECL_keyword_zippy2.push_back(loc_dummy_str);
						}
						zeilen = zeilen + 1;
					}
				}
				// save the tuning keyword
				if (this->SplittedString[0] == "TUNING")
				{
					int current_rec = 1;
					// the tuning keyword is special since it consists of always three data records
					//    for(int current_rec = 0; current_rec<3; current_rec++) 
					//    {
					do
					zeilen = zeilen + 1;
					while (eclDataFile->Data[i + zeilen] == ""); // look for empty rows between the Keyword and the items

					while (eclDataFile->Data[i + zeilen] != "")
					{
						if (eclDataFile->Data[i + zeilen] == "/")
							current_rec++;
						if (current_rec <= 3)
						{
							string delimiter = " ";
							SplittedString.clear();
							SplitStrings(eclDataFile->Data[i + zeilen], delimiter);
							if (eclDataFile->Data[i + zeilen].substr(0, 2) != "--") // ignore comments
							{
								std::string loc_dummy_str = "";
								loc_dummy_str = eclDataFile->Data[i + zeilen];
								ECL_keyword_tuning.push_back(loc_dummy_str);
							}
							zeilen = zeilen + 1;
						}
						else
							break;
					}
				}

				// save the TSTEP keyword
				if (this->SplittedString[0] == "TSTEP")
				{
					do
					zeilen = zeilen + 1;
					while (eclDataFile->Data[i + zeilen] == ""); // look for empty rows between the Keyword and the items

					while ((eclDataFile->Data[i + zeilen] != "") && (eclDataFile->Data[i + zeilen] != "/"))
					{
						string delimiter = " ";
						SplittedString.clear();
						SplitStrings(eclDataFile->Data[i + zeilen], delimiter);
						if (eclDataFile->Data[i + zeilen].substr(0, 2) != "--") // ignore comments
						{
							std::string loc_dummy_str = "";
							loc_dummy_str = eclDataFile->Data[i + zeilen];
							ECL_keyword_tstep.push_back(loc_dummy_str);
						}
						zeilen = zeilen + 1;
					}
				}
				// save the WCONPROD keyword
				if (this->SplittedString[0] == "WCONPROD")
				{
					do
					zeilen = zeilen + 1;
					while (eclDataFile->Data[i + zeilen] == ""); // look for empty rows between the Keyword and the items

					while ((eclDataFile->Data[i + zeilen] != "") && (eclDataFile->Data[i + zeilen] != "/"))
					{
						string delimiter = " ";
						SplittedString.clear();
						SplitStrings(eclDataFile->Data[i + zeilen], delimiter);
						if (eclDataFile->Data[i + zeilen].substr(0, 2) != "--") // ignore comments
						{
							std::string loc_dummy_str = "";
							loc_dummy_str = eclDataFile->Data[i + zeilen];
							ECL_keyword_wconprod.push_back(loc_dummy_str);
						}
						zeilen = zeilen + 1;
					}
				}
			}

			//// Read Well Data
			//if (this->existWells == true)
			//{
			//    if (this->SplittedString[0] == "WCONINJE")
			//    {
			//        long zeilen = 0;
			//        do {
			//            zeilen = zeilen + 1;
			//        }
			//        while (TextFile->Data[i + zeilen] != ""); // look for possible empty rows between the Keyword and the items
			//        jj = 0;

			//        /* WTP: this loop starts at the wrong position and therefore does not read the data
			//        for(long unsigned j = l + zeilen + 1;
			//            j < (l + this->ecl_well.size() + zeilen + 1); j++)*/
			//        for(long unsigned j = i + 1;j < (i + this->ecl_well.size() + 1); j++)
			//        {
			//            if (this->actual_time == 0)
			//            {    
			//                in.str(TextFile->Data[j]);
			//                in >> name;
			//                in >> this->ecl_well[jj]->phase;
			//                in >> this->ecl_well[jj]->open_flag;
			//                in >> this->ecl_well[jj]->control_mode;
			//                in >> dummy_rate;
			//                in >> dummy_zeile;
			//                if (in.gcount() != 0)
			//                    rest_zeile += dummy_zeile;
			//                in.clear();
			//            }
			//            jj++;
			//        }
			//    }
			//}
		}
	}

	// Now, depending on the just read data write all needed variables into the corresponding vector
	//fill the phase vector with the wetting phase first
	if (this->E100 == true)
	{
		if (Water_phase_exists == true)
			this->Phases.push_back("WATER");
		else
		{
			// The software is not checked if water is not the wetting phase
			std::cout << " ERROR: The program is canceled because no water phase is defined in ECLIPSE!" << "\n";
			system("Pause");
			exit(0);
		}

		if (Oil_phase_exists == true)
			this->Phases.push_back("OIL");
		if (Gas_phase_exists == true)
			this->Phases.push_back("GAS");


		// set the correct indicies
		bool found_comp = false;
		for (unsigned int j = 0; j < m_pcs->vec_component_pcs_names.size(); j++)
		if (vec_components_ECL_OGS[j][1] == 1 || vec_components_ECL_OGS[j][2] == 1 || vec_components_ECL_OGS[j][3] == 1)
		{
			found_comp = true;
			vec_components_ECL_OGS[j][0] = 1;
		}
		if (found_comp == false && verbosity > 1)
			std::cout << "  MESSAGE: No component is transported." << "\n";

		//Define necessary Variables
		// WTP: To add functionalty this should be more felxible...espacially for E300
		// e.g. the Variables should be read from the .DATA (under RPTRST) and written in "Variables"
		// of course there should be a consistency check afterwards to check if every needed variable is applicable

		//The total pressure is only defined once
		this->Variables.push_back("PRESSURE");
		this->Variables.push_back("ACAQNUM");

		for (unsigned int i = 0; i < this->Phases.size(); i++)
		{
			if (this->Phases[i].compare("WATER") == 0)
			{
				if (reservoir_conditions)
				{
					this->Variables.push_back("FLRWATI+");
					this->Variables.push_back("FLRWATJ+");
					this->Variables.push_back("FLRWATK+");
				}
				else
				{
					this->Variables.push_back("FLOWATI+");
					this->Variables.push_back("FLOWATJ+");
					this->Variables.push_back("FLOWATK+");
				}
				this->Variables.push_back("PWAT");
				this->Variables.push_back("WAT_DEN");
				this->Variables.push_back("BW");
				this->Variables.push_back("WAT_VISC");
				if (this->Phases.size() > 1)
					this->Variables.push_back("SWAT");
			}
			if (this->Phases[i].compare("OIL") == 0)
			{
				if (reservoir_conditions)
				{
					this->Variables.push_back("FLROILI+");
					this->Variables.push_back("FLROILJ+");
					this->Variables.push_back("FLROILK+");
				}
				else
				{
					this->Variables.push_back("FLOOILI+");
					this->Variables.push_back("FLOOILJ+");
					this->Variables.push_back("FLOOILK+");
				}

				this->Variables.push_back("POIL");
				this->Variables.push_back("OIL_DEN");
				this->Variables.push_back("BO");
				this->Variables.push_back("OIL_VISC");
				if (this->Phases.size() > 1)
					this->Variables.push_back("SOIL");
			}
			if (this->Phases[i].compare("GAS") == 0)
			{
				if (reservoir_conditions)
				{
					this->Variables.push_back("FLRGASI+");
					this->Variables.push_back("FLRGASJ+");
					this->Variables.push_back("FLRGASK+");
				}
				else
				{
					this->Variables.push_back("FLOGASI+");
					this->Variables.push_back("FLOGASJ+");
					this->Variables.push_back("FLOGASK+");
				}

				this->Variables.push_back("PGAS");
				this->Variables.push_back("GAS_DEN");
				this->Variables.push_back("BG");
				this->Variables.push_back("GAS_VISC");
				if (this->Phases.size() > 1)
					this->Variables.push_back("SGAS");
			}
		}

		// variable for reading capillary pressures in the case of Eclipse E100
		if (static_cast<int>(this->Phases.size()) == 2) // wtp change to c++ casts
		{
			if ((this->Phases[0] == "WATER") && (this->Phases[1] == "OIL"))
				this->Variables.push_back("PCOW");
			else
			{
				if ((this->Phases[0] == "OIL") && (this->Phases[1] == "GAS"))
				{
					this->Variables.push_back("PCOG");
					this->Variables.push_back("RS"); //gas saturation in oil in the black oil mode
					this->Variables.push_back("BG");
					this->Variables.push_back("GAS_DEN"); //gas density in the black oil mode
				}
				else
					std::cout <<
					" ERROR: Currently only Oil-Water or Oil-Gas is available in Eclipse E100."
					<< "\n";
                std::cout << flush;
				exit(0);
			}
		}
		if (static_cast<int>(this->Phases.size()) == 3) // wtp change to c++ casts
		{
			this->Variables.push_back("PCOW");
			this->Variables.push_back("PCOG");
			this->Variables.push_back("RS"); //gas saturation in oil in the black oil mode
			//this->Variables.push_back("GAS_DEN"); //gas density in the black oil mode

			this->phase_shift_flag = true;    // WTP moved due to identical if-conditions
		}

		//SB adapt to three phase mode in black oil using oil as water and gas as co2
		//if (static_cast<int>(this->Phases.size()) == 3)  // wtp c++ casts
		//    // three phases present, but use only oil as water and gas as CO2
		//    this->phase_shift_flag = true;
	}
	else
	{
		//E300 mode -> it is assumed that water and CO2 phase are considered --> NOT ANYMORE :)

		// WTP: Not used and obsolete anyway since all available comps are saved under vec_cnames_ecl()
		//    this->Components.push_back("H2O");
		//    this->Components.push_back("CO2");
		//    this->Components.push_back("NaCl");

		//phases
		if (Water_phase_exists == true)
			this->Phases.push_back("WATER");
		if (Gas_phase_exists == true)
			this->Phases.push_back("GAS");
		if (Oil_phase_exists == true)
			this->Phases.push_back("OIL");
		//variables
		this->Variables.push_back("PRESSURE"); // Total pressure
		this->Variables.push_back("PWAT"); // water phase pressure
		this->Variables.push_back("PGAS"); // gas phase pressure
		this->Variables.push_back("POIL"); // oil phase pressure

		this->Variables.push_back("PCOW"); // oil water capillary pressure
		this->Variables.push_back("PCOG"); // oil gas capillary pressure

		//this->Variables.push_back("TEMP");    // WTP: Grid Block Temperature

		this->Variables.push_back("SWAT"); // water phase saturation
		this->Variables.push_back("SGAS"); // gas phase saturation
		this->Variables.push_back("SOIL"); // gas phase saturation

		//this->Variables.push_back("RS"); // aqueous gas concentration WTP: not usable unless CO2STORE is selected or Oil Phase is present

		this->Variables.push_back("RPORV"); // pore volume of the cell

		this->Variables.push_back("DENW"); // water phase density
		this->Variables.push_back("DENG"); // gas phase density
		this->Variables.push_back("DENO"); // oil phase density

		this->Variables.push_back("BW");    // Water formation volume factor - used for density conversion res -> surf    
		this->Variables.push_back("BO");    // Oil formation volume factor - needed for CO2Store optione...
		this->Variables.push_back("BG");    // Gas formation volume factor

		this->Variables.push_back("VWAT");    // water phase viscosity
		this->Variables.push_back("VGAS");    // gas phase viscosity
		this->Variables.push_back("VOIL");    // oil phase viscosity

        this->Variables.push_back("TEMP");      //WTP active grid block temperatures
        
		// Write all variables needed for the gaseous composition (mass fraction) -> YFW
		// always needed if gas phase is present (CO2STORE and GASSOL)
		for (unsigned int k = 0; k < this->vec_components_ECL_OGS.size(); k++)
		{
			if (vec_components_ECL_OGS[k][0] != -1)
			{
				std::stringstream ss;
				ss << vec_components_ECL_OGS[k][0];
				std::string dummy_str = "YFW" + ss.str();
				Variables.push_back(dummy_str);
				ss.clear();
			}
		}

		// this might be needed in writedatabacktoeclipse()
		for (unsigned int k = 0; k < this->vec_components_ECL_OGS.size(); k++)
		{
			if (vec_components_ECL_OGS[k][0] != -1)
			{
				std::stringstream ss;
				ss << vec_components_ECL_OGS[k][0];
				std::string dummy_str = "MLSC" + ss.str();
				Variables.push_back(dummy_str);
				ss.clear();
			}
		}

		//if(this->option_CO2STORE==true) // not needed for CO2Store but for oil phase in general now
		//{
		// Write all variables needed for the liquid composition (mass fraction) -> XFW
		// only needed if CO2STORE  option is selected
		for (unsigned int k = 0; k < this->vec_components_ECL_OGS.size(); k++)
		{
			if (vec_components_ECL_OGS[k][0] != -1)
			{
				std::stringstream ss;
				ss << vec_components_ECL_OGS[k][0];
				std::string dummy_str = "XFW" + ss.str();
				Variables.push_back(dummy_str);
				ss.clear();
			}
		}
		//}

		// Write all variables needed for the inter-block component flow for the components
		for (unsigned int k = 0; k < this->vec_components_ECL_OGS.size(); k++)
		{
			if (vec_components_ECL_OGS[k][0] != -1)
			{
				long dummy_index;
				std::stringstream ss;
				ss << vec_components_ECL_OGS[k][0];
				std::string dummy_str = "FLOC" + ss.str();
				Variables.push_back(dummy_str + "I+");
				Variables.push_back(dummy_str + "J+");
				Variables.push_back(dummy_str + "K+");
				dummy_index = this->GetVariableIndex("FLOC" + ss.str() + "I+");
				vec_floI_indicies.push_back(dummy_index);
				dummy_index = this->GetVariableIndex("FLOC" + ss.str() + "J+");
				vec_floJ_indicies.push_back(dummy_index);
				dummy_index = this->GetVariableIndex("FLOC" + ss.str() + "K+");
				vec_floK_indicies.push_back(dummy_index);
				ss.clear();
			}
		}

		if (reservoir_conditions)
		{
			// ecl keywords for flow output at reservoir conditions for water, gas & oil
			this->Variables.push_back("FLRWATI+");
			this->Variables.push_back("FLRWATJ+");
			this->Variables.push_back("FLRWATK+");

			this->Variables.push_back("FLRGASI+");
			this->Variables.push_back("FLRGASJ+");
			this->Variables.push_back("FLRGASK+");

			this->Variables.push_back("FLROILI+");
			this->Variables.push_back("FLROILJ+");
			this->Variables.push_back("FLROILK+");
		}
		else
		{
			// Variables for outpur of surface values
			this->Variables.push_back("FLOWATI+");
			this->Variables.push_back("FLOWATJ+");
			this->Variables.push_back("FLOWATK+");

			this->Variables.push_back("FLOGASI+");
			this->Variables.push_back("FLOGASJ+");
			this->Variables.push_back("FLOGASK+");

			this->Variables.push_back("FLOOILI+");
			this->Variables.push_back("FLOOILJ+");
			this->Variables.push_back("FLOOILK+");
		}

	}

	this->numberOutputParameters = static_cast<int>(Variables.size());
	//WTP this->times = 1;

	finish = clock();
	time = (static_cast<double>(finish)-static_cast<double>(start)) / CLOCKS_PER_SEC;
	if(verbosity > 2)
		std::cout << "            Time: " << time << " seconds." << "\n";

	return 1;
}
/*-------------------------------------------------------------------------
   GeoSys - Function: ReadWellData
   Task: Read the well information and the time schedule for injection
   Return: true or false
   Programming: 02/2011 KB
   Modification: 12/2013 WTP added functionality for multi comp injection schemes
   -------------------------------------------------------------------------*/
void CECLIPSEData::ReadWellData(std::string Filename_Wells)
{
	structWell* ecl_single_well;
	char line[MAX_ZEILEN];
	std::string tempstring;
	std::string line_string;
	int comps = 0;
	streampos position;
	std::string tempstring_name;
	std::string tempstring_fluidtype;
	std::string tempstring_controlmode;
	std::string tempstring_rate;
	std::string tempstring_limit;
	std::vector <std::string> tempstringvec_mixture;
	std::string dummy_str = "";
	std::vector <vector <double> > inj_gas_composition;
	double tempdouble_time = 0.0;

	clock_t start, finish;
	double time;

	start = clock();

	if (verbosity > 2)
	std::cout << "        ReadWellData()";

    if (this->ReadWellData_called == true)
    {
		if (verbosity > 2)
			std::cout << " - Function already called once, skip this time" << "\n";
    }
    else
    {
        // get the number of components within eclipse if it is not E300
        if (this->E100 == true || this->option_CO2STORE == true)
            comps = 1;
        else
            comps = this->no_comps_ECL;
        // open filename
        ifstream in(Filename_Wells.data(), ios::in);

        if (!in)
            std::cout << " ERROR: File containing well data not found." << "\n";
        else
        while (!in.eof())
        {
            in.getline(line, MAX_ZEILEN);
            position = in.tellg();
            line_string = line;

            // reached the end of the file
            //if (line_string.find("#STOP") != std::string::npos)
            //    std::cout << "ok" << "\n";

            // Look for well names...
            if (line_string.find("$NAME") != string::npos) //subkeyword found
            {
                // create a new instance of a well structure and save the well name
                ecl_single_well = new structWell;
                in >> tempstring_name;
                ecl_single_well->name = tempstring_name;
                in.clear();
            }

            //TIMECURVE
            if (line_string.find("$TIMECURVE") != string::npos) //subkeyword found
            {
                while (tempstring.find("#") != 0)
                {
                    position = in.tellg();
                    in.getline(line, MAX_ZEILEN);
                    line_string = line;
                    if (line_string.find("#") != string::npos)
                        break;
                    if (line_string.find("$") != string::npos)
                        break;
                    in.seekg(position);
                    // first get the time
                    in >> tempdouble_time;
                    ecl_single_well->time.push_back(tempdouble_time);
                    // now get the fluid type
                    in >> tempstring_fluidtype;
                    ecl_single_well->phase.push_back(tempstring_fluidtype);
                    // now the run mode
                    in >> tempstring_controlmode;
                    ecl_single_well->control_mode.push_back(tempstring_controlmode);
                    // now the flow rate
                    in >> tempstring_rate;
                    ecl_single_well->rate.push_back(tempstring_rate);
                    // limiting pressure
                    in >> tempstring_limit;
                    ecl_single_well->limit.push_back(tempstring_limit);
                    // component mix
                    for (int i = 0; i < comps; i++)
                    {
                        std::string dummy_str;
                        in >> dummy_str;
                        tempstringvec_mixture.push_back(dummy_str);
                    }
                    ecl_single_well->comp_mix.push_back(tempstringvec_mixture);
                    tempstringvec_mixture.clear();

                    //if ((double(ecl_single_well->time.begin())) != "0"){
                    //    std::cout << "Warning: first item of time in timecurve is not 0!" << "\n";
                    //    break;
                    //}
                    in.clear();
                    in.getline(line, MAX_ZEILEN);
                    line_string = line;
                }
                if (!tempstring_name.empty())
                    // save the data
                    this->ecl_well.push_back(ecl_single_well);
            }
        }

        ReadWellData_called = true;
        finish = clock();
        time = (static_cast<double>(finish)-static_cast<double>(start)) / CLOCKS_PER_SEC;
		if (verbosity > 2)
			std::cout << "                Time: " << time << " seconds." << "\n";
    }
}

/*-------------------------------------------------------------------------
   GeoSys - Function: ReadPositionBoundaryCondition
   Task: Read the position of boundary conditions
   Return: true or false
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
bool CECLIPSEData::ReadPositionBoundaryCondition(std::string Filename)
{
	std::string tempstring;
	CBoundaryConditions* m_BoundaryConditions = NULL;
	CReadTextfiles_ECL* TextFile;
	std::stringstream in;
	bool error = false;
	int zeile;
	int number;
	int i1;
	int i2;
	int j1;
	int j2;
	int k1;
	int k2;
	int index;
	std::string boundary_position;

	TextFile = new CReadTextfiles_ECL;
	error = TextFile->Read_Text(Filename);

	if (error == true)
	{
		std::cout << " ERROR: The file: " << Filename <<
			"could not been read." << "\n";
		//system("Pause");
		exit(0);
	}

	//Read in the input file and search for the keyword AQUANCON
	for (long l = 0; l < TextFile->NumberOfRows; l++)
	{
		tempstring = TextFile->Data[l].substr(0, 8);
		zeile = 1;
		if (tempstring.compare("AQUANCON") == 0)
		{
			//End of the block is characterised with "" or "/"

			while ((TextFile->Data[l + zeile] != "") &&
				(TextFile->Data[l + zeile] != "/"))
			{
				std::string delimiter = " ";
				SplittedString.clear();
				SplitStrings(TextFile->Data[l + zeile], delimiter);

				if (TextFile->Data[l + zeile].substr(0, 2) != "--") // ignore comments
				{
					number = atoi(SplittedString[0].data());
					//check if the number of bc is plausibel
					if (number > 200)
					{
						std::cout <<
							"There were more than 200 bc found in the Eclipse model. Please make sure that the definitions at the keyword AQUANCON are correct!"
							<< "\n";
						return false;
					}
					i1 = atoi(SplittedString[1].data());
					i2 = atoi(SplittedString[2].data());
					j1 = atoi(SplittedString[3].data());
					j2 = atoi(SplittedString[4].data());
					k1 = atoi(SplittedString[5].data());
					k2 = atoi(SplittedString[6].data());
					boundary_position = SplittedString[7].substr(1, 2);
					index = -1;
					//order of cells in output files: z1, y1, x1..xn; z1, y2, x1..xn; z2, y1, x1..xn
					for (int k = k1; k <= k2; k++)
					{
						for (int j = j1; j <= j2; j++)
						for (int i = i1; i <= i2; i++)
						{
							index = index + 1;
							m_BoundaryConditions =
								new CBoundaryConditions(); // new member of eclblock
							m_BoundaryConditions->index = index;
							m_BoundaryConditions->number =
								number;
							m_BoundaryConditions->
								boundary_position =
								boundary_position;
							for (long m = 0; m < this->elements;
								m++)
							if ((this->eclgrid[m]->
								column == i) &&
								(this->eclgrid[m]->row
								==
								j) &&
								(this->eclgrid[m]->
								layer ==
								k))
								m_BoundaryConditions
								->connected_element
								= this->
								eclgrid[m
								]->
								index;

							this->BC.push_back(
								m_BoundaryConditions);
							long index_boundary =
								long(this->BC.size() - 1);
							this->eclgrid[m_BoundaryConditions
								->connected_element]
								->
								ConnectedBoundaryCondition.
								push_back(
								index_boundary);
						}
					}
					if (m_BoundaryConditions->connected_element == -1)
						return false;
				}
				zeile = zeile + 1;
			}
			l = TextFile->NumberOfRows;
		}
	}
	return true;
}
/*-------------------------------------------------------------------------
   GeoSys - Function: ReadBoundaryData
   Task: Read the flow data for boundary conditions
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
bool CECLIPSEData::ReadBoundaryData(int index_boundary, vector <string> Data)
{
	int counter;
	int number_values = 4;
	int number_cell;
	int number_value;

	//4 numbers are given for each boundary cell (order of the cells at the moment unknown!!!)
	counter = -1;
	for (unsigned long i = 0; i < Data.size(); i++)
	{
		this->SplittedString.clear();
		if (i < (Data.size() - 1))
		{
			std::string SplitBoundaryData[3];

			//this->SplitBoundaryData.push_back(Data[i].substr(3,20));
			//this->SplitBoundaryData.push_back(Data[i].substr(26,20));
			//this->SplitBoundaryData.push_back(Data[i].substr(49,20));

			SplitBoundaryData[0] = Data[i].substr(3, 20);
			SplitBoundaryData[1] = Data[i].substr(26, 20);
			SplitBoundaryData[2] = Data[i].substr(49, 20);

			//std::cout << SplitBoundaryData.size() << "\n";

			//for (unsigned int j = 0; j < this->SplittedString.size(); j++){
			for (unsigned int j = 0; j < 3; j++)
			{
				counter = counter + 1;
				//WTP number_cell = int(counter / number_values);
				number_cell = static_cast<int>(counter / number_values);
				number_value = counter - (number_cell * number_values);
				//search for the corresponding BC (the same boundary number and the index indentical to the counter
				for (unsigned int k = 0; k < this->BC.size(); k++)
				{
					if ((this->BC[k]->number == index_boundary) &&
						(this->BC[k]->index == number_cell))
					{
						//this->BC[k]->value[number_value] = atof(this->SplittedString[j].data());
						this->BC[k]->value[number_value] = atof(
							SplitBoundaryData[j].data());
						//WTP k = int(this->BC.size());
						k = static_cast<int>(this->BC.size());
					}
				}
			}
		}
		else
		{
			this->SplittedString.clear();
			this->SplitStrings(Data[i], " ");
			if (i == 0)
			if (SplittedString.size() != 3)
				return false;
			//string SplitBoundaryData[2];

			//SplitBoundaryData[0] = Data[i].substr(3,20);
			//SplitBoundaryData[1] = Data[i].substr(26,20);

			//for (unsigned int j = 0; j < 2; j++){
			for (unsigned int j = 0; j < this->SplittedString.size(); j++)
			{
				counter = counter + 1;
				//WTP number_cell = int(counter / number_values);
				number_cell = static_cast<int>(counter / number_values);
				number_value = counter - (number_cell * number_values);
				//search for the corresponding BC (the same boundary number and the index indentical to the counter
				for (unsigned int k = 0; k < this->BC.size(); k++)
				{
					if ((this->BC[k]->number == index_boundary) &&
						(this->BC[k]->index == number_cell))
					{
						this->BC[k]->value[number_value] = atof(
							this->SplittedString[j].data());
						//this->BC[k]->value[number_value] = atof(SplitBoundaryData[j].data());
						//WTP k = int(this->BC.size());
						k = static_cast<int>(this->BC.size());
					}
				}
			}
		}
	}

	return true;
}
/*-------------------------------------------------------------------------
   GeoSys - Function: ReadEclipseData
   Task: Read the data of the current time step of the Eclipse model from output file
   Return: nothing
   Programming: 09/2009 BG
   Modification: 04/2014 WTP minor changes, added support for 3 phases
   -------------------------------------------------------------------------*/
void CECLIPSEData::ReadEclipseData(std::string Filename)
{
	//WTP (void)Timestep; // unused
	vector<string> files = vector<string>();
	CReadTextfiles_ECL* TextFile;
	std::string tempstring;
	bool saturation_water, saturation_gas, saturation_oil;
	clock_t start, finish;
	double time;

	start = clock();

	saturation_water = saturation_gas = saturation_oil = false;

	if (verbosity > 2)
		std::cout << "        ReadEclipseData() ";
	//get memory for data structures
	if (Data == NULL) // allocate first time
	{
		Data = (double**)malloc((this->elements) * sizeof(double*));
		for (long i = 0; i < (this->elements); i++)
			Data[i] = (double*)malloc(this->numberOutputParameters * sizeof(double));
		//for(long i = 0; i < (this->elements); i++)
		//	for(long j = 0; j < this->times; j++)
		//		Data[i][j] = (double *) malloc(this->numberOutputParameters*sizeof(double));

		//Deallokieren
		/*Data2 = (double **) malloc(1000000 * sizeof(double*));
		   for(long i = 0; i < (1000000); i++)
		   Data2[i] = (double *) malloc(100 * sizeof(double));

		   for(long i=0;i<1000000;i++) free(Data2[i]);

		   free(Data2);*/
	}

	//initialise data structure
	for (long i = 0; i < (this->elements); i++)
	for (int k = 0; k < this->numberOutputParameters; k++)
		Data[i][k] = 0;

	//Reads the text file (usually 60 bytes in one row of the eclipse output
	bool Error;
	TextFile = new CReadTextfiles_ECL;
	Error = TextFile->Read_Text(Filename);
	if (Error == true)
	{
		std::cout << "\n";
		std::cout << " ERROR: The program is canceled in ReadEclipseData()" << "\n";
		//system("Pause");
		exit(0);
	}
	for (long j = 0; j < TextFile->NumberOfRows; j++)
	{
		//if (j==565 || j ==1384) {
		//	cout << TextFile->Daten[j] << "\n";
		//}

		double Multiplier = 1;
		if (TextFile->Data[j].length() > 0 && TextFile->Data[j].substr(1, 1) == "'")
		{
			std::string tempstring = TextFile->Data[j].substr(2, 8);
			tempstring = tempstring.substr(0, tempstring.find_first_of(" "));
			for (unsigned int k = 0; k < this->Variables.size(); k++)
			if (tempstring.compare(this->Variables[k]) == 0)
			{
				if (this->Variables[k].compare("ACAQNUM") == 0)
				{
					int index_boundary = atoi(
						TextFile->Data[j + 1].substr(0, 12).data());
					this->SplittedString.clear();
					this->SplitStrings(TextFile->Data[j + 2], " ");
					//WTP int temprows = int(atoi(this->SplittedString[2].data()) / 3) + 1;
					int temprows = static_cast<int>(atoi(this->SplittedString[2].data()) / 3) + 1;
					vector <string> temp_Daten;
					for (int l = 0; l < temprows; l++)
						temp_Daten.push_back(TextFile->Data[j + 3 +
						l]);
					if (this->ReadBoundaryData(index_boundary,
						temp_Daten) == false)
					{
						std::cout <<
							" ERROR: Problem reading boundary data." <<
							"\n";
						//system("Pause");
						exit(0);
					}
				}
				else
				{
					//convert pressure units from bar to Pa
					if ((this->Variables[k].compare("PRESSURE") ==
						0) ||
						(this->Variables[k].compare("PCOW") == 0) ||
						(this->Variables[k].compare("PCOG") == 0) ||
						(this->Variables[k].compare("PWAT") == 0) ||
						(this->Variables[k].compare("PGAS") == 0) ||
						(this->Variables[k].compare("POIL") == 0))
						Multiplier = 100000.0;

					//convert viscosity units from cP to Pa*s
					if ((this->Variables[k].compare("GAS_VISC") == 0) ||
						(this->Variables[k].compare("WAT_VISC") == 0) ||
						(this->Variables[k].compare("OIL_VISC") == 0) ||
						(this->Variables[k].compare("VWAT") == 0) ||
						(this->Variables[k].compare("VGAS") == 0) ||
						(this->Variables[k].compare("VOIL") == 0))
						Multiplier = 1.0E-3;

					// consider different orientation of the z-axis in Eclipse and GeoSys
					if (reservoir_conditions)
					{
						if ((this->Variables[k].compare("FLRWATK+") == 0) ||
							(this->Variables[k].compare("FLROILK+") == 0) ||
							(this->Variables[k].compare("FLRGASK+") == 0))
							Multiplier = -1.0;
					}
					else
					{
						if ((this->Variables[k].compare("FLOWATK+") == 0) ||
							(this->Variables[k].compare("FLOOILK+") == 0) ||
							(this->Variables[k].compare("FLOGASK+") == 0))
							Multiplier = -1.0;
					};
					
					if (this->Variables[k].compare("RS") == 0)
						Multiplier = 1.0;  // Unit: m^3 gas per m^3 oil

					if (this->Variables[k].compare("GAS_DEN") == 0)
						Multiplier = 1.0;  // Unit: kg/m^3

					//set a flag of which phase saturation was red
					if (this->Variables[k].compare("SWAT") == 0)
						saturation_water = true;
					if (this->Variables[k].compare("SGAS") == 0)
						saturation_gas = true;
					if (this->Variables[k].compare("SOIL") == 0)
						saturation_oil = true;

					// read number of datapoints
					double tempNumber =
						atoi(TextFile->Data[j].substr(12, 13).data());
					long temprows = 0;
					if (tempNumber == this->elements)
					{
						//WTP temprows = int(ceil(tempNumber / 4.0));
						temprows = static_cast<int>(ceil(tempNumber / 4.0));
						//Read data for identified variable
						long rowindex = 0;
						for (long l = j + 1; l < j + temprows + 1;
							l++)
						{
							if (l < (j + temprows))
							{
								std::string Dataline[4]; // Definition einen strings mit 4 Arrays, Arrays immer in eckigen Klammern definieren
								Dataline[0] =
									TextFile->Data[l].
									substr(2,
									15);
								Dataline[1] =
									TextFile->Data[l].
									substr(19,
									15);
								Dataline[2] =
									TextFile->Data[l].
									substr(36,
									15);
								Dataline[3] =
									TextFile->Data[l].
									substr(53,
									15);
								for (unsigned int m = 0;
									m < 4; m++)
								{
									rowindex =
										rowindex +
										1; //zaehlt die gesplitteten Zeilen in der Zeile, in Variante KB4 entspricht das 4
									this->Data[rowindex
										- 1][k]
										= atof(
										Dataline[m]
										.data())
										*
										Multiplier;
									//std::cout << "Element: " << rowindex << " Variable: " << this->Variables[k] 
									//<< " Value: " << this->Data[rowindex-1][k] << "\n";

								}
							}
							else
							{
								std::string delimiter = " ";
								SplittedString.clear();
								SplitStrings(
									TextFile->Data[l],
									delimiter);
								//std::cout << TextFile->Daten[l] << "\n";
								for (unsigned int m = 0;
									m <
									SplittedString.size();
								m++)
								{
									rowindex =
										rowindex +
										1;
									//std::cout << rowindex-1 << " " << k << " " << " " << this->eclgrid[rowindex-1]->x_barycentre <<  " " 
									//<< this->eclgrid[rowindex-1]->y_barycentre  << " " 
									//<< this->eclgrid[rowindex-1]->z_barycentre  << " " 
									//<< atof(SplittedString[m].data()) << "\n";
									this->Data[rowindex
										- 1][k]
										= atof(
										SplittedString
										[m].
										data()) *
										Multiplier;
									//std::cout << "Element: " << rowindex << " Variable: " << this->Variables[k] << " Value: " 
									//<< this->Data[rowindex-1][k] << "\n";

								}
							}
						}
						j = j + temprows;
						//a[k][1] = j + 1;
					}
					else
					{
						// there are inactive cells
						//std::cout << "Error while reading the data file. The number of data points doesn't fit to the grid." << "\n";
						//WTP temprows = int(ceil(tempNumber / 4.0));
						temprows = static_cast<int>(ceil(tempNumber / 4.0));
						//Read data for identified variable
						long rowindex = 0;
						for (long l = j + 1; l < j + temprows + 1;
							l++)
						{
							if (l < (j + temprows))
							{
								std::string Dataline[4]; // Definition einen strings mit 4 Arrays, Arrays immer in eckigen Klammern definieren
								Dataline[0] =
									TextFile->Data[l].
									substr(3,
									15);
								Dataline[1] =
									TextFile->Data[l].
									substr(20,
									15);
								Dataline[2] =
									TextFile->Data[l].
									substr(38,
									15);
								Dataline[3] =
									TextFile->Data[l].
									substr(55,
									15);
								//string delimiter=" ";
								//SplittedString.clear();
								//SplitStrings(TextFile->Daten[l], delimiter);
								//cout << TextFile->Daten[l] << "\n";
								for (unsigned int m = 0;
									m < 4; m++)
								{
									rowindex =
										rowindex +
										1;
									while (this->
										eclgrid[
											rowindex
												-
												1]->
												active ==
												0 &&
												rowindex <
												this->
												elements)
												rowindex =
												rowindex
												+ 1;
											if (rowindex <
												this->elements)
												this->Data[
													rowindex
														-
														1][
															k]
																=
																atof(
																Dataline
																[
																	m
																]
															.
																data())
																*
																Multiplier;
											else if (verbosity > 1)
												std::cout <<
												" WARNING: 1 data point couldn't be allocated to a grid point"
												<<
												"\n";
								}
							}
							else
							{
								std::string delimiter = " ";
								SplittedString.clear();
								SplitStrings(
									TextFile->Data[l],
									delimiter);
								//cout << TextFile->Daten[l] << "\n";
								for (unsigned int m = 0;
									m <
									SplittedString.size();
								m++)
								{
									rowindex =
										rowindex +
										1;
									//std::cout << rowindex-1 << " " << timestep << " " << k << " " << " " 
									//<< this->eclgrid[rowindex-1]->x_barycentre 
									//<<  " " << this->eclgrid[rowindex-1]->y_barycentre  << " " <
									//< this->eclgrid[rowindex-1]->z_barycentre  << " " 
									//<< atof(SplittedString[m].data()) << "\n";
									this->Data[rowindex
										- 1][k]
										= atof(
										SplittedString
										[m].
										data()) *
										Multiplier;
									//    std::cout << "Element: " << rowindex << " Variable: " << this->Variables[k] << " Value: " 
									//<< this->Data[rowindex-1][timestep][k] << "\n";

								}
							}
						}
						j = j + temprows;
					}
				}
			}
		}
	}

	// Release Textfile object - sb
	//delete(TextFile->Daten);
	delete (TextFile);

	if (this->E100 == true)
	{
		//--------------------------------------------------------------------------------------------
		// calculating phase pressure in the case of Eclipse E100 from pressure and capillary pressure
		//WTP if (int(this->Phases.size()) == 1)
		if (static_cast<int>(this->Phases.size()) == 1)
		{
			int index_pwat, index_pgas, index_poil, index_pressure;
			int index_aim = -1;
			index_pwat = this->GetVariableIndex("PWAT");
			index_pgas = this->GetVariableIndex("PGAS");
			index_poil = this->GetVariableIndex("POIL");
			index_pressure = this->GetVariableIndex("PRESSURE");

			if (index_pwat >= 0)
				index_aim = index_pwat;
			if (index_pgas >= 0)
				index_aim = index_pgas;
			if (index_poil >= 0)
				index_aim = index_poil;

			// Transfer the pressure
			for (int i = 0; i < this->elements; i++)
				this->Data[i][index_aim] = this->Data[i][index_pressure];
		}
		//assumption: if water and oil -> pressure = oil pressure
		//WTP if (int(this->Phases.size()) == 2)
		if (static_cast<int>(this->Phases.size()) == 2)
		{
			// water and oil
			if ((this->Phases[0] == "WATER") && (this->Phases[1] == "OIL"))
			{
				// Get the variable index
				int index_pwat, index_poil, index_pressure, index_pcap;
				index_pwat = this->GetVariableIndex("PWAT");
				index_poil = this->GetVariableIndex("POIL");
				index_pressure = this->GetVariableIndex("PRESSURE");
				index_pcap = this->GetVariableIndex("PCOW");

				// Transfer the pressure
				for (int i = 0; i < this->elements; i++)
				{
					this->Data[i][index_pwat] = this->Data[i][index_pressure];
					this->Data[i][index_poil] = this->Data[i][index_pwat] +
						this->Data[i][index_pcap];
				}
			}
			if ((this->Phases[0] == "WATER") && (this->Phases[1] == "GAS"))
			{
				std::cout << "\n";
				std::cout <<
					" ERROR: A GAS-WATER System can not be considered with E100 and OGS"
					<< "\n";
				//system("Pause");
				exit(0);
			}
			if ((this->Phases[0] == "OIL") && (this->Phases[1] == "GAS"))
			{
				// Get the variable index
				int index_pgas, index_poil, index_pressure, index_pcap;
				index_pgas = this->GetVariableIndex("PGAS");
				index_poil = this->GetVariableIndex("POIL");
				index_pressure = this->GetVariableIndex("PRESSURE");
				index_pcap = this->GetVariableIndex("PCOG");

				// Transfer the pressure
				for (int i = 0; i < this->elements; i++)
				{
					this->Data[i][index_poil] = this->Data[i][index_pressure];
					this->Data[i][index_pgas] = this->Data[i][index_poil] +
						this->Data[i][index_pcap];
				}
			}
		}
		if (static_cast<int>(this->Phases.size()) == 3)
		{
			//std::cout << "\n";
			//std::cout <<
			//	"Currently not more than 2 phases are considered for reading eclipse pressure"
			//	<< "\n";
			//std::cout <<
			//	"Three phases are only valid if a water saturation exists only at the boundaries of the model domain!!";
			//system("Pause");
			//exit(0);

			int index_pwat, index_poil, index_pgas, index_pressure,
				index_pcap_oil_water, index_pcap_oil_gas;
			index_pwat = this->GetVariableIndex("PWAT");
			index_poil = this->GetVariableIndex("POIL");
			index_pgas = this->GetVariableIndex("PGAS");
			index_pressure = this->GetVariableIndex("PRESSURE");
			index_pcap_oil_water = this->GetVariableIndex("PCOW");
			index_pcap_oil_gas = this->GetVariableIndex("PCOG");

			// Transfer the pressure
			for (int i = 0; i < this->elements; i++)
			{
				this->Data[i][index_pwat] = this->Data[i][index_pressure];
				this->Data[i][index_poil] = this->Data[i][index_pwat] +
					this->Data[i][index_pcap_oil_water];
				this->Data[i][index_pgas] = this->Data[i][index_poil] +
					this->Data[i][index_pcap_oil_gas];
			}
		}

		//--------------------------------------------------------------------------------------------
		// calculating phase saturation for all phases

		// 2 existing phases
		if (this->Phases.size() == 2)
		{
			int index_aim = -1, index_source = -1;
			int index_swat, index_sgas, index_soil;
			index_swat = this->GetVariableIndex("SWAT");
			index_sgas = this->GetVariableIndex("SGAS");
			index_soil = this->GetVariableIndex("SOIL");

			if ((index_swat >= 0) && (saturation_water == true))
			{
				index_source = index_swat;
				if (index_sgas >= 0)
					index_aim = index_sgas;
				if (index_soil >= 0)
					index_aim = index_soil;
			}
			if ((index_sgas >= 0) && (saturation_gas == true))
			{
				index_source = index_sgas;
				if (index_swat >= 0)
					index_aim = index_swat;
				if (index_soil >= 0)
					index_aim = index_soil;
			}
			if ((index_soil >= 0) && (saturation_oil == true))
			{
				index_source = index_soil;
				if (index_swat >= 0)
					index_aim = index_swat;
				if (index_sgas >= 0)
					index_aim = index_sgas;
			}

			// Calculating the phase saturation
			for (int i = 0; i < this->elements; i++)
				this->Data[i][index_aim] = 1 - this->Data[i][index_source];
		}
		// 3 existing phases
		if (this->Phases.size() == 3)
		{
			int index_aim = -1, index_source1 = -1, index_source2 = -1;
			int index_swat, index_sgas, index_soil;
			index_swat = this->GetVariableIndex("SWAT");
			index_sgas = this->GetVariableIndex("SGAS");
			index_soil = this->GetVariableIndex("SOIL");

			if (saturation_water == false)
			{
				index_aim = index_swat;
				index_source1 = index_sgas;
				index_source2 = index_soil;
			}
			if (saturation_gas == false)
			{
				index_aim = index_sgas;
				index_source1 = index_swat;
				index_source2 = index_soil;
			}
			if (saturation_oil == false)
			{
				index_aim = index_soil;
				index_source1 = index_sgas;
				index_source2 = index_swat;
			}
			// Calculating the phase saturation
			for (int i = 0; i < this->elements; i++)
				this->Data[i][index_aim] = 1 - this->Data[i][index_source1] -
				this->Data[i][index_source2];
			//std::cout << "Element: " << i << " RS: " << this->Data[i][this->GetVariableIndex("RS")] << " SWAT: " << this->Data[i][index_swat];
			//std::cout  << " SOIL: " << this->Data[i][index_soil] << " SGAS: " << this->Data[i][index_sgas] << "\n";
		}

		if (this->Phases.size() > 3)
		{
			std::cout <<
				" ERROR: More than 3 phases are not considered for reading eclipse data"
				<< "\n";
			//system("Pause");
			exit(0);
		}
	}
	//--------------------------------------------------------------------------------------------
	// output of vertical flow data
	//for (unsigned long i = 0; i < this->elements; i++) {
	//	cout << "Q(k) " << this->Data[i][0][2] << "\n";
	//}
	//std::cout << "Fertig";

	//Test output
	//vector <string> vec_string;
	//ostringstream temp;
	//tempstring = "Element; X; Y; Z";
	//for (unsigned int k = 0; k < this->Variables.size(); k++) {
	//	tempstring = tempstring + "; " + this->Variables[k];
	//}
	//vec_string.push_back(tempstring);
	// Loop over all elements
	//for (long i = 0; i < this->elements; i++){
	//	CECLIPSEBlock *elem = this->eclgrid[i];
	//	temp.str(""); temp.clear(); temp << i; tempstring = temp.str();
	//	temp.str(""); temp.clear(); temp << elem->x_barycentre; tempstring += "; " + temp.str();
	//	temp.str(""); temp.clear(); temp << elem->y_barycentre; tempstring += "; " + temp.str();
	//	temp.str(""); temp.clear(); temp << elem->z_barycentre; tempstring += "; " + temp.str();

	//	for (unsigned int k = 0; k < this->Variables.size(); k++) {
	//		temp.str(""); temp.clear(); temp.precision(12); temp << this->Data[i][k]; tempstring += "; " + temp.str();
	//	}
	//	//if (i == 470) {
	//	//	std::cout << i;
	//	//}
	//	vec_string.push_back(tempstring);
	//}  // end element loop

	//Test Output
	//int position = Filename.find_last_of("\\");
	//string path = Filename.substr(0,position);
	//position = path.find_last_of("\\");
	//path = path.substr(0,position);
	// //temp.str(""); temp.clear(); temp << timestep; tempstring = temp.str();
	//string aus_file = path + "\\CheckDataRedIn_0.csv";
	//ofstream aus;
	//aus.open(aus_file.data(),ios::out);
	//for (unsigned int i = 0; i < vec_string.size(); i++) {
	//	aus << vec_string[i] << "\n";
	//}
	//aus.close();

	finish = clock();
	time = (static_cast<double>(finish)-static_cast<double>(start)) / CLOCKS_PER_SEC;
	if (verbosity > 2)
		std::cout << "                    Time: " << time << " seconds." << "\n";

}

/*-------------------------------------------------------------------------
   GeoSys - Function: ConvertEclipseDataToUniformUnits
   Task: Converting the Eclipse units for the dissolved gas components to uniform units
   Return: nothing
   Programming: 01/2011 BG
   Modification: 02/2014 WTP: generalized for multi comps
   03/2014 WTP: renamed the function
   -------------------------------------------------------------------------*/
void CECLIPSEData::ConvertEclipseDataToUniformUnits(CRFProcess* m_pcs, long Timestep)
{
	clock_t start, finish;
	double time;
	//double comp_water;
	//double ConcWater;
	//int variable_index_ReservoirDensity_Water = -1;
	int variable_index_ReservoirDensity_Gas = -1;
	int variable_index_ReservoirDensity_Oil = -1;
	//int variable_index_FVF_Gas = -1;
	int variable_index_FVF_Water = -1;
	int variable_index_Saturation_Water = -1;
	int variable_index_Saturation_Gas = -1;
	int variable_index_Saturation_Oil = -1;
	int variable_index_reservoir_volume = -1;
	//const double epsilon = 1.0E-15;
//	CRFProcess* m_pcs; //KB0714

	start = clock();

	if (verbosity > 2)
		std::cout << "        ConvertingEclipseUnits()";

	// clear vectors
	vec_comps_total_indices.clear();
	vec_comps_gas_ECL_indicies.clear();
	vec_comps_oil_ECL_indicies.clear();

	// WTP: Maybe we should just do this once? 
	// Now get the indicies for the concentrations of all transported components
	if (this->E100 == true)
	{
		variable_index_FVF_Water = this->GetVariableIndex("BO");
		variable_index_ReservoirDensity_Gas = this->GetVariableIndex("GAS_DEN");

		for (unsigned int i = 0; i < this->vec_components_ECL_OGS.size(); i++)
		{
			if (vec_components_ECL_OGS[i][0] != -1)
			{
				comp_water_ECL_indicies = this->GetVariableIndex("RS");
				break;
			}
		}
	}
	else
	{  	// E300 case
		variable_index_ReservoirDensity_Gas = this->GetVariableIndex("DENG");
		variable_index_Saturation_Water = this->GetVariableIndex("SWAT");
		variable_index_Saturation_Gas = this->GetVariableIndex("SGAS");
		variable_index_Saturation_Oil = this->GetVariableIndex("SOIL");
		variable_index_ReservoirDensity_Oil = this->GetVariableIndex("DENO");
		variable_index_reservoir_volume = this->GetVariableIndex("RPORV");

		for (unsigned int i = 0; i < this->vec_components_ECL_OGS.size(); i++)
		{
			if (vec_components_ECL_OGS[i][0] != -1)
			{
				// get the correct component number to assemble the keyphrase
				std::stringstream ss;
				ss << vec_components_ECL_OGS[i][0];
				long dummy_index = this->GetVariableIndex("MLSC" + ss.str());
				vec_comps_total_indices.push_back(dummy_index);
				dummy_index = this->GetVariableIndex("YFW" + ss.str());
				vec_comps_gas_ECL_indicies.push_back(dummy_index);
				dummy_index = this->GetVariableIndex("XFW" + ss.str());
				vec_comps_oil_ECL_indicies.push_back(dummy_index);
				vec_ecl_comp_indicies.push_back(vec_components_ECL_OGS[i][0]);
			}
		}
	}

	// Get the indicies for the dissolved components in water
	if (this->vec_OGS_process_index_comps_water.size() == 0)
	{
		for (unsigned int l = 0; l < this->vec_components_ECL_OGS_pcs_names.size(); l++)
		{
			CRFProcess *n_pcs = NULL;            // if the component is transported get a process
			for (int k = 0; k < static_cast<int>(pcs_vector.size()); k++)  // get the pcs names
			{
				n_pcs = pcs_vector[k];
				if (n_pcs->nod_val_name_vector[0] == this->vec_components_ECL_OGS_pcs_names[l][1])
				{
					// if the process is found, store the index
					this->vec_OGS_process_index_comps_water.push_back(std::make_pair(vec_components_ECL_OGS_pcs_names[l][1], k));
					////wtp_debug_3003
					//std::cout << " OGS process index comps water: k: " << k << " pcs idx: " << vec_components_ECL_OGS_pcs_names[l][1] << "\n";
				}
			}
		}
	}
	// Get the indicies for the components in gas
	if (this->vec_OGS_process_index_comps_gas.size() == 0)
	{
		for (unsigned int l = 0; l < this->vec_components_ECL_OGS_pcs_names.size(); l++)
		{
			CRFProcess *n_pcs = NULL;            // if the component is transported get a process
			for (int k = 0; k < static_cast<int>(pcs_vector.size()); k++)  // get the pcs names
			{
				n_pcs = pcs_vector[k];
				if (n_pcs->nod_val_name_vector[0] == this->vec_components_ECL_OGS_pcs_names[l][2])
				{
					// if the process is found, store the index
					this->vec_OGS_process_index_comps_gas.push_back(std::make_pair(vec_components_ECL_OGS_pcs_names[l][2], k));
					////wtp_debug_3003
					//std::cout << " OGS process index comps gas: k: " << k << " pcs idx: " << vec_components_ECL_OGS_pcs_names[l][2] << "\n";
				}
			}
		}
	}
	// Get the indicies for the components in oil
	if (this->vec_OGS_process_index_comps_oil.size() == 0)
	{
		for (unsigned int l = 0; l < this->vec_components_ECL_OGS_pcs_names.size(); l++)
		{
			CRFProcess *n_pcs = NULL;            // if the component is transported get a process
			for (int k = 0; k < static_cast<int>(pcs_vector.size()); k++)  // get the pcs names
			{
				n_pcs = pcs_vector[k];
				if (n_pcs->nod_val_name_vector[0] == this->vec_components_ECL_OGS_pcs_names[l][3])
				{
					// if the process is found, store the index
					this->vec_OGS_process_index_comps_oil.push_back(std::make_pair(vec_components_ECL_OGS_pcs_names[l][3], k));
					////wtp_debug_3003
					//std::cout << " OGS process index comps oil: k: " << k << " pcs idx: " << vec_components_ECL_OGS_pcs_names[l][3] << "\n";
				}
			}
		}
	}

	// E100 case: Only one component can be in the water phase since the gasphase only consits of one single component. Furthermore, the oil phase 
	// must be used as water phase if dissolved components shall be incorporated
	if (this->E100 == true)
	{
		// first get the molar weight of the component
		double MolarWeight_Comp = 0.;

		for (unsigned int j = 0; j < cp_vec.size(); j++)
		{
			std::string dummy_str = cp_vec[j]->compname;    // get the component name
			if (dummy_str == vec_OGS_process_index_comps_water[0].first)
			{
				MolarWeight_Comp = cp_vec[j]->molar_weight;    // save the parameters and break the loop to save time
				////wtp_debug_3003
				//std::cout << " Molar weight of component: " << MolarWeight_Comp << "\n";
				break;
			}
		}
		////wtp_debug_3003
		//std::cout << " Element# \t FVF_Water \t comp_water \t ConcWater \t ConcGas " << "\n";
		for (long i = 0; i < this->elements; i++)
		{
			// dummy variable structures
			std::vector <double> dummy_vec_dbl;
			std::vector <double> dummy_vec_dbl_2;
			// get the formation volume factor of the liquid phase
			double FVF_Water = this->Data[i][variable_index_FVF_Water];

			////WTP DEBUG /////
			//std::cout << "hard coded FVF_water to 1 in ConvertEclipseToUniformUnits for E100!!!!" << "\n";
			//FVF_Water = 1.0;
			////WTP DEBUG /////

			// Get the value of dissolved component in water
			double comp_water = this->Data[i][comp_water_ECL_indicies];
			double ConcWater = comp_water * SurfaceDensity_Gas_E100 / (FVF_Water * MolarWeight_Comp);
			double ConcGas = this->Data[i][variable_index_ReservoirDensity_Gas] / MolarWeight_Comp; // WTP: Maybe also calculate the reservoir density from the surface density and the FVF?
			
			////wtp_debug_3003
			//std::cout << " " << i << " " << FVF_Water << " " << comp_water << " " << ConcWater << " " << ConcGas << "\n";

			if (Timestep > 1 || (m_pcs->iter_outer_cpl > 0))
			{
				vec_CompConc_Water_elements[i][0] = ConcWater;
				vec_CompConc_Gas_elements[i][0] = ConcGas;
			}
			else
			{
				dummy_vec_dbl.push_back(ConcWater);
				dummy_vec_dbl_2.push_back(ConcGas);
			}
			if (Timestep == 1 && m_pcs->iter_outer_cpl == 0)
			{
				vec_CompConc_Water_elements.push_back(dummy_vec_dbl);
				vec_CompConc_Gas_elements.push_back(dummy_vec_dbl_2);
			}
		} // end of elements loop
	} // end of E100 case

	// E300 case: Multiple Components in multiple phases...
	if (this->E100 != true)
	{
		// 1. since we definatly have dissolved components (-> check & adjust for only dissolved comps in oil & gas) 
		//    get the total molar density for each transported component
		// 2. Calculate the moles in the oil Phase for each component
		// 3. Calculate the moles in the gas Phase for each component
		// 4. Subtstract 3 and 4 from 1
		// Attention: This is ecl only and has nothing to do with OGS.
		// OGS only comes into play when we save the data

		for (long j = 0; j < this->elements; j++) // run over all elements and quantifiy the ammount of moles in each applicable phase
		{
			// dummy variable structures
			std::vector <double> dummy_vec_dbl;
			std::vector <double> dummy_vec_dbl_2;
			std::vector <double> dummy_vec_dbl_3;
			std::vector <double> dummy_vec_dbl_mol;
			int current_run = -1;
			for (unsigned int i = 0; i < this->vec_components_ECL_OGS.size(); i++)
			{
				if (vec_components_ECL_OGS[i][0] != -1) // if component is transported in any phase get the data
				{
					++current_run;
					// first get the molar weight of the current component
					double MolarWeight_Comp = this->vec_ecl_mw_comps[vec_ecl_comp_indicies[current_run] - 1];    // vec_ecl_mw_comps has the same length as the vec_components_ECL_OGS but not the same ORDER
					if (j == 0)
						this->Components.push_back(vec_cnames_ecl[vec_ecl_comp_indicies[current_run] - 1]);
					//double MolarWeight_Comp = this->vec_ecl_mw_comps[current_run];    // vec_ecl_mw_comps has the same length as the vec_components_ECL_OGS but not the same ORDER

					double CompConc_Gas = 0.;
					double CompConc_Oil = 0.;
					double CompConc_Water = 0.;
					double moles_in_gas = 0.;
					double moles_in_oil = 0.;
					double moles_in_water = 0.;
					double TotalMolarDensity = this->Data[j][vec_comps_total_indices[current_run]];
					double res_porevolume = this->Data[j][variable_index_reservoir_volume];
					double TotalMoles = TotalMolarDensity * res_porevolume * 1000;
					//if a gas phase exists, get the moles of the component in it
					if (this->Gas_phase_exists == true)
					{
						double weight_fraction_gas = this->Data[j][vec_comps_gas_ECL_indicies[current_run]];
						double saturation_gas = this->Data[j][variable_index_Saturation_Gas];
						double res_density_gas = this->Data[j][variable_index_ReservoirDensity_Gas];
						if (saturation_gas > 0.)
						{
							moles_in_gas = weight_fraction_gas * res_density_gas * saturation_gas * res_porevolume / MolarWeight_Comp;
							if (moles_in_gas > TotalMoles)
							{
								if ((moles_in_gas - TotalMoles) / TotalMoles > 0.0001 &&  verbosity > 1)
								{
									std::cout << "\n WARNING: Unrealistic data in element " << j << " for component: " << this->vec_cnames_ecl[vec_ecl_comp_indicies[current_run] - 1] << "\n";
									std::cout << " Moles in gas exceed total moles by " << moles_in_gas - TotalMoles << " (" << (moles_in_gas - TotalMoles) / TotalMoles << "%). \n";
									std::cout << " The value is corrected.";
								}
								moles_in_gas = TotalMoles;
							}
							CompConc_Gas = moles_in_gas / res_porevolume / saturation_gas;
						}
						else
						{
							moles_in_gas = 0.;
							CompConc_Gas = 0.;
						}

						// wtp debug
						//std::cout << moles_in_gas << "\n";

						if (Timestep > 1 || (m_pcs->iter_outer_cpl > 0))
						{
							vec_CompConc_Gas_elements[j][current_run] = CompConc_Gas;
							vec_mComp_gas_elements[j][current_run] = moles_in_gas;
						}
						else
						{
							dummy_vec_dbl_2.push_back(CompConc_Gas);
							dummy_vec_dbl_mol.push_back(moles_in_gas);
						}
					}
					//same for the oil phase
					if (this->Oil_phase_exists == true)
					{
						double weight_fraction_oil = this->Data[j][vec_comps_oil_ECL_indicies[current_run]];
						double saturation_oil = this->Data[j][variable_index_Saturation_Oil];
						double res_density_oil = this->Data[j][variable_index_ReservoirDensity_Oil];
						if (saturation_oil > 0.)
						{
							moles_in_oil = weight_fraction_oil * res_density_oil * saturation_oil * res_porevolume / MolarWeight_Comp;
							if (moles_in_oil > TotalMoles)
							{
								if ((moles_in_oil - TotalMoles) / TotalMoles > 0.01 && verbosity > 1)
								{
									std::cout << "\n WARNING: Unrealistic data in element " << j << " for component: " << this->vec_cnames_ecl[vec_ecl_comp_indicies[current_run] - 1] << "\n";
									std::cout << " Moles in oil exceed total moles by " << moles_in_oil - TotalMoles << " (" << (moles_in_oil - TotalMoles) / TotalMoles << "%). \n";
									std::cout << " The value is corrected.";
								}
								moles_in_oil = TotalMoles;
							}
							CompConc_Oil = moles_in_oil / res_porevolume / saturation_oil;
						}
						else
						{
							moles_in_oil = 0.;
							CompConc_Oil = 0.;
						}
						if (Timestep > 1 || (m_pcs->iter_outer_cpl > 0))
						{
							vec_CompConc_Oil_elements[j][current_run] = CompConc_Oil;
						}
						else
						{
							dummy_vec_dbl_3.push_back(CompConc_Oil);
						}
					}
					// the water phase consits of the "rest"
					if (Water_phase_exists == true)
					{
						double saturation_water = this->Data[j][variable_index_Saturation_Water];
						if (saturation_water > 0.)
						{
							moles_in_water = TotalMolarDensity * res_porevolume * 1000 - moles_in_gas - moles_in_oil;
							CompConc_Water = moles_in_water / res_porevolume / saturation_water;
						}
						else
						{
							moles_in_water = 0.;
							CompConc_Water = 0.;
						}
						if (Timestep > 1 || (m_pcs->iter_outer_cpl > 0))
						{
							vec_CompConc_Water_elements[j][current_run] = CompConc_Water;
						}
						else
						{
							dummy_vec_dbl.push_back(CompConc_Water);
						}
					}
					//current_run++;
				}
			}
			if (Timestep == 1 && m_pcs->iter_outer_cpl == 0)
			{
				if (Water_phase_exists == true)
					vec_CompConc_Water_elements.push_back(dummy_vec_dbl);
				if (Gas_phase_exists == true)
				{
					vec_CompConc_Gas_elements.push_back(dummy_vec_dbl_2);
					vec_mComp_gas_elements.push_back(dummy_vec_dbl_mol);
				}
				if (Oil_phase_exists == true)
					vec_CompConc_Oil_elements.push_back(dummy_vec_dbl_3);
			}
		}// end of elements loop
	}

	finish = clock();
	time = (static_cast<double>(finish)-static_cast<double>(start)) / CLOCKS_PER_SEC;
	if (verbosity > 2)
		std::cout << "              Time: " << time << " seconds." << "\n";
}

/*-------------------------------------------------------------------------
   GeoSys - Function: ReplaceASectionInFile
   Task: Determine corresponding nodes and elements between Geosys and Eclipse
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
bool CECLIPSEData::ReplaceASectionInFile(std::string Filename,
	string Keyword,
	vector <std::string> Data,
	bool CheckLengthOfSection)
{
	CReadTextfiles_ECL* TextFile;
	std::string tempstring;
	bool success = false;

	//Reads the text file (usually 60 bytes in one row of the eclipse output
	bool Error;
	TextFile = new CReadTextfiles_ECL;
	Error = TextFile->Read_Text(Filename);
	if (Error == true)
	{
		std::cout << " ERROR: Problem in ReplaceASectionInFile()." << "\n";
		//system("Pause");
		exit(0);
	}

	//scan file for the keyword and replace the following data
	for (long i = 0; i < TextFile->NumberOfRows; i++)
	if (TextFile->Data[i].substr(0, Keyword.length()) == Keyword)
	{
		//count data rows in the file and compare it to the rows in the new Data array, if equal than replace data
		success = true;
		long zeilen = 0;
		if (CheckLengthOfSection == true)
		{
			do {
				zeilen = zeilen + 1;
			} while (TextFile->Data[i + zeilen].length() > 1);
			zeilen = zeilen - 1;
		}
		else
			zeilen = long(Data.size());
		if (zeilen == long(Data.size()))
			//Replace Data in TextFile variable
		for (long j = 0; j < zeilen; j++)
			TextFile->Data[i + 1 + j] = Data[j];

		else
		{
			std::cout <<
				" ERROR: Replacing a section in the ECLIPSE input file is not possible because the section length doesn't fit!"
				<< "\n";
			success = false;
			return success;
		}
	}
	//Rewrite the file
	CWriteTextfiles_ECL* OutputFile;
	OutputFile = new CWriteTextfiles_ECL;
	OutputFile->Write_Text(Filename, TextFile->Data);

	//Release memory
	delete (TextFile);

	return success;
}
/*-------------------------------------------------------------------------
   GeoSys - Function: ReplaceASectionInData
   Task: Replaces a given section in the ECL *.Fxx file data array
   Return: nothing
   Programming: 02/2014 WTP
   Modification:
   -------------------------------------------------------------------------*/
bool CECLIPSEData::ReplaceASectionInData(CReadTextfiles_ECL* eclFFile,
	string Keyword,
	vector <std::string> Data,
	bool CheckLengthOfSection)
{
	std::string tempstring;
	bool success = false;
	//WTP std::cout << pathECLFFile << "\n"; //KB 

	//scan data structure for the keyword and replace the following data
	//for (long i = 0; i < eclFFile->NumberOfRows; i++)
	for (size_t i = 0; i < eclFFile->Data.size(); i++)
		if (eclFFile->Data[i].substr(0, Keyword.length()) == Keyword)
		{
			//count data rows in the file and compare it to the rows in the new Data array, if equal than replace data
			success = true;
			long zeilen = 0;
			if (CheckLengthOfSection == true)
			{
				do {
					zeilen = zeilen + 1;
				}
				//while (TextFile->Data[i + zeilen] != "");
				while (eclFFile->Data[i + zeilen].length() > 1);
				zeilen = zeilen - 1;
			}
			else
				zeilen = long(Data.size());
			if (zeilen == long(Data.size()))
				//Replace Data in TextFile variable
			for (long j = 0; j < zeilen; j++)
				eclFFile->Data[i + 1 + j] = Data[j];
			else
			{
				std::cout <<
					" ERROR: Replacing a section in the ECLIPSE input file is not possible because the section length doesn't fit!"
					<< "\n";
				success = false;
				return success;
			}
		}
	return success;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: WriteIncludeFile
   Task: Generates or adds a section in an include file
   Return: nothing
   Programming: 05/2011 CB
   Modification: re-added 03/2013 wtp
   -------------------------------------------------------------------------*/
bool CECLIPSEData::WriteIncludeFile(std::string Filename, std::string Keyword, vector <std::string> Data, bool append) {
	bool success = true;
	ofstream aus;

	if (append)
		aus.open(Filename.c_str(), ios::app);
	else
		aus.open(Filename.c_str());

	aus << Keyword << "\n";
	//	for (long i = 0; i < Data.size(); i++) {
	for (std::size_t i = 0; i < Data.size(); i++) { // WTP
		aus << Data[i] << "\n";
	}
	aus << "/" << "\n";
	aus.close();
	return success;
};
/*-------------------------------------------------------------------------
   GeoSys - Function: ReplaceWellRate
   Task: Replace Injection Rate with Data from .well-File
   Return: bool flag (not used properly yet)
   Programming: 02/2011 KB
   Modification: 11/2013 WTP: added support for multiple comps
   -------------------------------------------------------------------------*/
bool CECLIPSEData::ReplaceWellRate(CReadTextfiles_ECL* eclDataFile)
{
	// General function description:
	// checks each timecurve of each well and adjuststs the injection or production rate/target according to the data found in the *.well file

	clock_t start, finish;
	double time;

	bool success = false;
	bool rewrite_file = false;
	bool injection = false;
	bool production = false;

	// new timestep
	std::vector <std::string> rate_new;
	std::vector <std::string> limit_new;
	std::vector <std::string> phase_new;
	std::vector <std::string> control_mode_new;

	// vector strucutre for component mix in E300
	std::vector <std::vector <std::string> > comp_mix_new;

	// vector of strings for output
	std::vector <std::string> new_schedule_section;

	string::size_type   str_start = 0;             // substring coordinates
	string::size_type   str_end;
	double epsilon = 0.0000001;

	// resize all vectors to the number of wells
	rate_new.resize(ecl_well.size());
	limit_new.resize(ecl_well.size());
	phase_new.resize(ecl_well.size());
	control_mode_new.resize(ecl_well.size());

	// for convinience get the component mix at the same 
	comp_mix_new.resize(ecl_well.size());

	if (verbosity > 2)
		std::cout << "        Updating well rates";
	start = clock();

	// First run over all available wells and check if the current time corresponds to any change in injection/production pattern
	for (unsigned int i = 0; i < this->ecl_well.size(); i++)
	{
		//now run over the time loop
		for (unsigned int j = 0; j < this->ecl_well[i]->time.size(); j++)
		{
			// check the current timestep

			//if(this->ecl_well[i]->time[j] < this->actual_time)
			if (this->ecl_well[i]->time[j] <= this->actual_time)
			{
				// if the timestep corresponds to an entry in the *.well file access the data & set a bool flag 
				rewrite_file = true;

				rate_new[i] = ecl_well[i]->rate[j];
				limit_new[i] = this->ecl_well[i]->limit[j];
				phase_new[i] = this->ecl_well[i]->phase[j];
				control_mode_new[i] = this->ecl_well[i]->control_mode[j];

				comp_mix_new[i].clear();
				for (unsigned int k = 0; k < this->ecl_well[i]->comp_mix[j].size(); k++)
					comp_mix_new[i].push_back(this->ecl_well[i]->comp_mix[j][k]);

				// test for injection or production and set the corresponding flags
				double rate_value_new = atof(rate_new[i].c_str());
				if (fabs(rate_value_new) < epsilon)
					rate_value_new = 0.0;
				if (rate_value_new >= 0.0)
				{
					injection = true;
					production = false;
				}
				if (rate_value_new < 0.0)
				{
					production = true;
					injection = false;
				}
				//break;
			}
			else
				break;// end if actual time < threshold
		} // end time loop
	} // end well loop

	// now assemble the new schedule section from scratch
	if (rewrite_file == true)
	{
		// first get the headline right
		new_schedule_section.push_back("SCHEDULE");
		new_schedule_section.push_back("");
		// now write the MESSAGES and the RPTRST  keyword
		if (this->ECL_keyword_messages.size() > 0)
		{
			new_schedule_section.push_back("MESSAGES");
			for (unsigned int i = 0; i<this->ECL_keyword_messages.size(); i++)
				new_schedule_section.push_back(ECL_keyword_messages[i]);
			new_schedule_section.push_back("");
		}

		if (this->ECL_keyword_rptrst.size() > 0)
		{
			new_schedule_section.push_back("RPTRST");
			for (unsigned int i = 0; i<this->ECL_keyword_rptrst.size(); i++)
				new_schedule_section.push_back(ECL_keyword_rptrst[i]);
			new_schedule_section.push_back("");
		}
		// write tuning or zippy2
		if (this->ECL_keyword_zippy2.size() > 0)
		{
			new_schedule_section.push_back("ZIPPY2");
			for (unsigned int i = 0; i<this->ECL_keyword_zippy2.size(); i++)
				new_schedule_section.push_back(ECL_keyword_zippy2[i]);
			new_schedule_section.push_back("");
		}
		if (this->ECL_keyword_tuning.size() > 0)
		{
			new_schedule_section.push_back("TUNING");
			for (unsigned int i = 0; i<this->ECL_keyword_tuning.size(); i++)
				new_schedule_section.push_back(ECL_keyword_tuning[i]);
			new_schedule_section.push_back("");
		}
		// now create the welpsecs keyword. For that each well has to be checked for any changes concerning the preffered phase
		if (this->ECL_keyword_welspecs.size() > 0)
		{
			new_schedule_section.push_back("WELSPECS");
			for (unsigned int i = 0; i < this->ecl_well.size(); i++)
			{
				// read string
				std::string well_data_old = this->ECL_keyword_welspecs[i];
				std::string substring_1 = "";
				std::string substring_2 = "";
				std::string token = "";
				str_start = 0;

				// set the correct token for the current phase
				if (well_data_old.find("GAS") != string::npos)
					token = "GAS";
				if (well_data_old.find("WATER") != string::npos)
					token = "WATER";
				if (well_data_old.find("OIL") != string::npos)
					token = "OIL";
				if (well_data_old.find("LIQ") != string::npos)
					token = "LIQ";

				// now assemble the correct line
				str_end = well_data_old.find(token, str_start);
				substring_1 = well_data_old.substr(str_start, str_end - str_start);
				str_start = str_end + token.length();
				substring_2 = well_data_old.substr(str_start, well_data_old.length() - str_start);
				// final record of keyword assemlbed    
				new_schedule_section.push_back(substring_1 + phase_new[i] + substring_2);
			}
			new_schedule_section.push_back("/");
			new_schedule_section.push_back("");
		}
		// now add the compdat keyword
		if (this->ECL_keyword_compdat.size() > 0)
		{
			new_schedule_section.push_back("COMPDAT");
			for (unsigned int i = 0; i < this->ECL_keyword_compdat.size(); i++)
				new_schedule_section.push_back(ECL_keyword_compdat[i]);
			new_schedule_section.push_back("/");
			new_schedule_section.push_back("");
		}

		// now assemble the injection and production data
		// first the injection stuff
		// change the component mixture and select the correnct comp for each well if needed
		if (injection == true)
		{
			// first check if it is an E100 case or not since E300 demands more keywords
			if (this->E100 != true)
			{
				std::vector <string> local_winjgas_key;        // local variables for temp storage
				std::string local_data_dummy;

				// first write the WELLSTRE keyword
				new_schedule_section.push_back("WELLSTRE");
				// the keyword expects one line of record for each injection well -> depends on the injection rate
				for (unsigned int i = 0; i < this->ecl_well.size(); i++)
				{
					// check injection rate
					//double rate_value_new = atof(rate_new[i].c_str());
					/*    if(rate_value_new > epsilon)
						{*/
					std::string dummy_name = "";
					std::string current_mix = "";
					// if the well name is given with ' ' delete them for the mixture name
					if (this->ecl_well[i]->name.substr(0, 1) == "'")
					{
						dummy_name = ecl_well[i]->name;
						dummy_name = dummy_name.substr(1, dummy_name.length() - 2);
					}
					else
						dummy_name = this->ecl_well[i]->name;
					// add up the string consisting of the phase mixture
					for (unsigned int j = 0; j < comp_mix_new[i].size(); j++)
						current_mix += comp_mix_new[i][j] + " ";
					// assemble the whole record
					new_schedule_section.push_back(dummy_name + "_mix" + " " + current_mix + "/");
					// save the data for writing the winjgas keyword afterwards
					local_winjgas_key.push_back(local_data_dummy = ecl_well[i]->name + " STREAM " + dummy_name + "_mix/");
					//}
				}
				new_schedule_section.push_back("/");
				new_schedule_section.push_back("");

				// quickly add the winjgas keyword based on the data saved during the wellstre keyword
				new_schedule_section.push_back("WINJGAS");
				for (unsigned int i = 0; i < local_winjgas_key.size(); i++)
					new_schedule_section.push_back(local_winjgas_key[i]);
				new_schedule_section.push_back("/");
				new_schedule_section.push_back("");
			}

			// now change the actual injection rate, keyword does not depend on E100 / E300 selection
			new_schedule_section.push_back("WCONINJE");
			for (unsigned int i = 0; i < this->ecl_well.size(); i++)
			{
				// check injection rate
				double rate_value_new = atof(rate_new[i].c_str());
				if (fabs(rate_value_new) < epsilon)
					rate_value_new = 0.;
				if (rate_value_new >= 0.)
				{
					// depending on the selected control mode the data has to be put together differently
					// please refer to the ECLIPSE manual for information
					if (control_mode_new[i] == "RATE")
						new_schedule_section.push_back(ecl_well[i]->name + " " + phase_new[i] + " OPEN " + control_mode_new[i]
						+ " " + rate_new[i] + " " + "1*" + " " + limit_new[i] + "/");
					if (control_mode_new[i] == "RESV")
						new_schedule_section.push_back(ecl_well[i]->name + " " + phase_new[i] + " OPEN " + control_mode_new[i]
						+ " " + "1*" + " " + rate_new[i] + " " + limit_new[i] + "/");
					if (control_mode_new[i] == "BHP")
						new_schedule_section.push_back(ecl_well[i]->name + " " + phase_new[i] + " OPEN " + control_mode_new[i]
						+ " " + "1*" + " " + "1*" + " " + limit_new[i] + "/");
				}
				//// if the injetion rate is 0 close the well
				//if(rate_value_new == 0.)
				//    new_schedule_section.push_back(ecl_well[i]->name + " " + phase_new[i] + " SHUT /");

			}
			new_schedule_section.push_back("/");
			new_schedule_section.push_back("");
		}

		// now do the same for the production wells
		if (production == true)
		{
			new_schedule_section.push_back("WCONPROD");
			for (unsigned int i = 0; i < this->ecl_well.size(); i++)
			{
				// check injection rate
				double rate_value_new = atof(rate_new[i].c_str());
				if (rate_value_new < -epsilon)
				{
					// adjust for the - sign
					std::string adj_rate = rate_new[i].substr(1, rate_new[i].length() - 1);

					// again the keyword has to be assembled differently depending on the selected control mode
					if (control_mode_new[i] == "ORAT")
						new_schedule_section.push_back(ecl_well[i]->name + " OPEN " + control_mode_new[i]
						+ " " + adj_rate + " " + "1*" + " " + "1*" + " " + "1*" + " " + "1*" + " " + limit_new[i] + "/");
					if (control_mode_new[i] == "GRAT")
						new_schedule_section.push_back(ecl_well[i]->name + " OPEN " + control_mode_new[i]
						+ " " + "1*" + " " + "1*" + " " + adj_rate + " " + "1*" + " " + "1*" + " " + limit_new[i] + "/");
					if (control_mode_new[i] == "LRAT")
						new_schedule_section.push_back(ecl_well[i]->name + " OPEN " + control_mode_new[i]
						+ " " + "1*" + " " + "1*" + " " + adj_rate + " " + "1*" + " " + "1*" + " " + limit_new[i] + "/");
					if (control_mode_new[i] == "RESV")
						new_schedule_section.push_back(ecl_well[i]->name + " OPEN " + control_mode_new[i]
						+ " " + "1*" + " " + "1*" + " " + "1*" + " " + adj_rate + " " + "1*" + " " + limit_new[i] + "/");
					if (control_mode_new[i] == "BHP")
						new_schedule_section.push_back(ecl_well[i]->name + " OPEN " + control_mode_new[i]
						+ " " + "1*" + " " + "1*" + " " + "1*" + " " + "1*" + " " + "1*" + " " + limit_new[i] + "/");
				}
			}
			new_schedule_section.push_back("/");
			new_schedule_section.push_back("");
		}

		// write TSTEP
		new_schedule_section.push_back("TSTEP");
		new_schedule_section.push_back("1*1 /");
		// write END
		new_schedule_section.push_back("");
		new_schedule_section.push_back("");
		new_schedule_section.push_back("END");


		// Now assemble the new *.DATA file (only virtual, will be written to the file later)
		std::string start_key = "SCHEDULE";
		//scan file for the keyword and cap the data structure at it
		for (long i = 0; i < eclDataFile->NumberOfRows; i++)
		if (eclDataFile->Data[i].substr(0, start_key.length()) == start_key)
		{
			eclDataFile->Data.resize(i);
			break;
		}
		// now add the newly made schedule section to it
		for (unsigned int i = 0; i < new_schedule_section.size(); i++)
			eclDataFile->Data.push_back(new_schedule_section[i]);
	}

	finish = clock();
	time = (static_cast<double>(finish)-static_cast<double>(start)) / CLOCKS_PER_SEC;
	if (verbosity > 2)
		std::cout << "           Time: " << time << " seconds." << "\n";

	// not used yet
	success = true;
	return success;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: CorrespondingElements
   Task: Determine corresponding nodes and elements between Geosys and Eclipse
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
bool CECLIPSEData::CorrespondingElements()
{
	MeshLib::CElem* m_ele = NULL;
	CFEMesh* m_msh = fem_msh_vector[0];

	//check if number of elements is equal
	if (long(m_msh->ele_vector.size()) != this->elements)
	{
		std::cout << " Error: The number of elements is not equal between Geosys and Eclipse." <<
			"\n";
		return false;
	}

	//set element vectors to the start value -1
	CorrespondingEclipseElement.resize(m_msh->ele_vector.size());
	for (unsigned long i = 0; i < m_msh->ele_vector.size(); i++)
		CorrespondingEclipseElement[i] = -1;
	CorrespondingGeosysElement.resize(m_msh->ele_vector.size());
	for (long j = 0; j < this->elements; j++)
		CorrespondingGeosysElement[j] = -1;

	for (unsigned long i = 0; i < m_msh->ele_vector.size(); i++)
	for (long j = 0; j < this->elements; j++)
	{
		m_ele = m_msh->ele_vector[i];
		double const* grav_c(m_ele->GetGravityCenter());
		//check if coordinates of the gravity centre are equal
		if ((grav_c[0] == this->eclgrid[j]->x_barycentre) &&
			(grav_c[1] == this->eclgrid[j]->y_barycentre) &&
			(grav_c[2] == -this->eclgrid[j]->z_barycentre))
		{
			CorrespondingEclipseElement[i] = j;
			CorrespondingGeosysElement[j] = i;
		}
	}

	//check if all values in the correspondingElementVector are larger than -1
	for (unsigned long i = 0; i < m_msh->ele_vector.size(); i++)
	{
/*		std::cout << "   CorrespondingEclipseElement[" << i << "] " <<
			CorrespondingEclipseElement[i] << "\n"*/;
		if (CorrespondingEclipseElement[i] < 0)
		{
			std::cout << " ERROR: No ECLIPSE element linked to the geosys element!" << "\n";
			return 0;
		}
	}
	for (long j = 0; j < this->elements; j++)
	if (CorrespondingGeosysElement[j] < 0)
	{
		std::cout << " ERROR: No geosys element linked to the ECLIPSE element!" << "\n";
		return 0;
	}
	return 1;
}

/*-------------------------------------------------------------------------
GeoSys - Function: CompareElementsGeosysEclipse
Task: Re-worked function of CompareElementsGeosysEclipse by SB
Return: bool value if there occured an error
Programming: 09/2012 SB
Modification: re-added 03/2013 wtp
compare for hex,pris,tet 06/2014 WB
-------------------------------------------------------------------------*/
bool CECLIPSEData::CompareElementsGeosysEclipse()
{
	MeshLib::CElem* m_ele = NULL;
	CFEMesh* m_msh = fem_msh_vector[0];
	//Math_Group::vec<MeshLib::CNode*> ele_nodes(8);
	//std::vector <MeshLib::CNode*> ele_nodes;
	MeshLib::CNode * m_cnode;
	clock_t start, finish;
	double time;
	/*double epsilon = 1e-7;*/
	double epsilon = 1e-2; // WB: differ starts at 3 digits.
	const double * pnt;
	//    const double * gc;
	start = clock();

	if (verbosity > 2)
		std::cout << "        CompareElements() ";

    if (CompareElementsGeosysEclipse_called == true)
    {
		if (verbosity > 2)
			std::cout << " - Function already called once, skip this time" << "\n";
    }
    else
    {
        //check if number of elements is equal
        //WB: number of OGS mesh == total number of corresponding OGS element regarding Eclipse Cell 
        for (size_t i = 0; i < long(this->eclgrid.size()); i++)
        {
            this->eclgridelenum += this->eclgrid[i]->correspondingelenum;
        }
        if (long(m_msh->ele_vector.size()) != this->eclgridelenum)
        {
            std::cout << " ERROR: The number of elements between Geosys and Eclipse is not equal!" << "\n";
            return false;
        }

        //Compare element/cell-wise data:node index,node coordinates
		for (size_t i = 0; long(i < this->eclgrid.size()); i++)
        {
            bool ElementIsEqual = true;
            // loop over all corresponding ele
            for (int j = 0; j < this->eclgrid[i]->correspondingelenum; j++)
            {
                m_ele = m_msh->ele_vector[this->eclgrid[i]->correspondingeleindex[j]];
                //loop over all nodes in one ele
				for (size_t k = 0; k < m_ele->GetNodesNumber(0); k++)
                {
                    m_cnode = m_ele->GetNode(k);
                    int equalnodeindex = -1;
                    //find the same index node
					for (size_t h = 0; h < this->eclgrid[i]->Nodeindex.size(); h++)
                    {
                        if (this->eclgrid[i]->Nodeindex[h] == static_cast<long>(m_cnode->GetIndex()))
                        {
                            equalnodeindex = h;
                            break;
                        }
                    }
                    //if no same index node or the same index nodes have different coordinates, return false
                    if (equalnodeindex == -1)
                        ElementIsEqual = false;
                    else
                    {
                        pnt = m_cnode->getData();
                        if (fabs(pnt[0] - this->eclgrid[i]->x_coordinates[equalnodeindex]) > epsilon ||
                            fabs(pnt[1] - this->eclgrid[i]->y_coordinates[equalnodeindex]) > epsilon ||
                            fabs(pnt[2] - -this->eclgrid[i]->z_coordinates[equalnodeindex]) > epsilon)
                            ElementIsEqual = false;
                    }
                }

                if (ElementIsEqual == false)
                    return 0;
                else
                {
                    //Get the gravity centre of the element from geosys
                    double const* gc(m_ele->GetGravityCenter());
                    //gc = m_ele->GetGravityCenter();
                    ////BW:Gravity Center is averaged by all corresponding elements gravity center
                    this->eclgrid[i]->x_barycentre += gc[0] / this->eclgrid[i]->correspondingelenum;
                    this->eclgrid[i]->y_barycentre += gc[1] / this->eclgrid[i]->correspondingelenum;
                    this->eclgrid[i]->z_barycentre += gc[2] / this->eclgrid[i]->correspondingelenum;
                    //BW:ECLGRID Block volume is added up by the volumet of all corresponding elements
                    this->eclgrid[i]->volume += m_ele->GetVolume();
                }
            }
        }

        this->CompareElementsGeosysEclipse_called = true;
        finish = clock();
        time = (static_cast<double>(finish)-static_cast<double>(start)) / CLOCKS_PER_SEC;
		if (verbosity > 2)
			std::cout << "                    Time: " << time << " seconds." << "\n";
        return 1;
    }

	return true;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: CreateFaces
   Task: Re-worked function of CreateFaces by SB
   Return: bool value if there occured an error
   Programming: 09/2012 SB
   Modification: re-added 03/2013 wtp
   add collapse grid 04/2014 BW
   -------------------------------------------------------------------------*/
bool CECLIPSEData::CreateFaces(string Projectname)

{
	CFaces* m_face = NULL;
	//MeshLib::CElem* m_element = NULL;
	vector <long> element_indices;
	CFEMesh* m_msh = fem_msh_vector[0];
	//Math_Group::vec <MeshLib::CNode*> element_nodes(8);
	//MeshLib::CNode* Node0 = NULL;
	//MeshLib::CNode* Node1 = NULL;
	//MeshLib::CNode* Node2 = NULL;
	//MeshLib::CNode* Node3 = NULL;
	//MeshLib::CNode* Node4 = NULL;
	//MeshLib::CNode* Node5 = NULL;
	//MeshLib::CNode* Node6 = NULL;
	//MeshLib::CNode* Node7 = NULL;
	vector <MeshLib::CNode*> vec_face_nodes;
	bool threenodesflag = false;
	long faceindex;
	long Blockindexnumber;
	clock_t start, finish;
	double time;
	bool bfaces_exist = true;
	bool generate_bfaces = false;

	start = clock();//TICK
	if (verbosity > 2)
		std::cout << "        CreateFaces() ";

    if (this->CreateFaces_called == true)
    {
		if (verbosity > 2)
			std::cout << " - Function already called once, skip this time" << "\n";
    }
    else
    {
        string Filename;
        string tempstring;
        CReadTextfiles_ECL* TextFile;
        std::stringstream in;

        bool error = false;
        //check if file *.bface exists, which stores the neighbours of each element
        Filename = Projectname + ".bfaces";
        TextFile = new CReadTextfiles_ECL;
		bfaces_exist = CheckIfFileExists(Filename);
		if (bfaces_exist == true)
		{
			error = TextFile->Read_Text(Filename);

			long tempindentification;

			//If it is a radial Model
			if (error == true)
				generate_bfaces = true;
			else
			{
				// Check Indentification Number, if false, give warning information
				in.str((string)TextFile->Data[1]);
				in >> tempindentification;
				in.clear();
				//if (this->IdentificationNumber != tempindentification)
				//{
				//	std::cout << "Warning: '*.bfaces' file is not generated at the same time as the '*.list' and '*.neighbours' file!" << "\n";
				//	std::cout << "         It might cause wrong information on the geometry of ECLIPSE Grid!" << "\n";
				//	//system("Pause");
				//	//exit(0);
				//}

				//Total number of faces attached to ECLIPSE BLOCKS
				in.str((string)TextFile->Data[2]);
				in >> this->FlowFaceNum;
				in.clear();
				//Initialize faces vector
				for (long tempi = 0; tempi < this->FlowFaceNum; tempi++)
				{
					m_face = new CFaces(static_cast<int>(this->Phases.size()));
					faces.push_back(m_face);
				}

				//Read the neighbours of each element
				for (long i = 3; i < TextFile->NumberOfRows; i++)
				{

					in.str((string)TextFile->Data[i]);
					in >> faceindex;
					faces[faceindex]->index = faceindex;
					in >> faces[faceindex]->nfnodes >> faces[faceindex]->category;

					//Nodal Index
					faces[faceindex]->NodeIndex.clear();
					faces[faceindex]->NodeIndex.resize(faces[faceindex]->nfnodes);
					for (int j = 0; j < faces[faceindex]->nfnodes; j++)
					{
						in >> faces[faceindex]->NodeIndex[j];
					}

                //Connected Block Index
                in >> Blockindexnumber;
                faces[faceindex]->connected_blocks.resize(Blockindexnumber);
                for (size_t j = 0; j < faces[faceindex]->connected_blocks.size(); j++)
                {
                    in >> faces[faceindex]->connected_blocks[j];
                    this->eclgrid[faces[faceindex]->connected_blocks[j]]->FaceIndex[faces[faceindex]->category - 1 + j];
                }

					in.clear();

				}
				delete (TextFile);
			}
		}
		else//Generate Face info based on geometry
			generate_bfaces = true;

		if (generate_bfaces == true)
		{
			//std::cout << "        -> Creating new *.bfaces file" << "\n";
            long int firsteclgrid = 0;
            //Find the first active cell index
            for (int m = 0; m < long(this->eclgrid.size()); m++)
            {
                if (this->eclgrid[m]->active != 0){
                    firsteclgrid = m;
                    break;
                }
            }

            // Attach first 6 Faces for the first ECIPSE active Cell
            // 'I+'
            int FaceIplus[] = { 1, 5, 3, 7 };
            this->FaceGenerator(firsteclgrid, FaceIplus, 1);
            // 'I-'
            int FaceIminus[] = { 0, 4, 2, 6 };
            this->FaceGenerator(firsteclgrid, FaceIminus, 2);
            // 'J+'
            int FaceJplus[] = { 2, 3, 6, 7 };
            this->FaceGenerator(firsteclgrid, FaceJplus, 3);
            // 'J-'
            int FaceJminus[] = { 0, 1, 4, 5 };
            this->FaceGenerator(firsteclgrid, FaceJminus, 4);
            // 'K+'
            int FaceKplus[] = { 4, 5, 6, 7 };
            this->FaceGenerator(firsteclgrid, FaceKplus, 5);
            // 'K-'
            int FaceKminus[] = { 0, 1, 2, 3 };
            this->FaceGenerator(firsteclgrid, FaceKminus, 6);

            // Attach all faces in the other ECLIPSE Cell
            for (long i = firsteclgrid + 1; i < static_cast<long>(this->eclgrid.size()); i++)
            {
                //if (i == 15879)
                //cout << i << '\n';
                // 'I+'
                int FaceIplus[] = { 1, 5, 3, 7 };
                this->FaceGenerator(i, FaceIplus, 1);
                // 'I-'
                int FaceIminus[] = { 0, 4, 2, 6 };
                this->FaceGenerator(i, FaceIminus, 2);
                // 'J+'
                int FaceJplus[] = { 2, 3, 6, 7 };
                this->FaceGenerator(i, FaceJplus, 3);
                // 'J-'
                int FaceJminus[] = { 0, 1, 4, 5 };
                this->FaceGenerator(i, FaceJminus, 4);
                // 'K+'
                int FaceKplus[] = { 4, 5, 6, 7 };
                this->FaceGenerator(i, FaceKplus, 5);
                // 'K-'
                int FaceKminus[] = { 0, 1, 2, 3 };
                this->FaceGenerator(i, FaceKminus, 6);
            }
            this->FlowFaceNum = long(this->faces.size());
        }

        //Generate More information about faces used for interface
        for (long i = 0; i < this->FlowFaceNum; i++)
        {
            // std::cout << " i : " << i << "\n";
            // Loop over 'I-'/'J-'/'K-' face which are not on the boundary
            //if (facecategory == 2 || facecategory == 4 || facecategory == 6)
            //if (this->eclgrid[Blockindex]->NeighbourElement[facecategory - 2] > -1)
            //	continue;
            //if (faceindex==130)
            //cout << faceindex << '\n';
            if (faces[i]->index != i)
            {
                std::cout << " ERROR: Index of faces does not match!!!" << '\n';
                exit(0);
            }

            //facedirection assign 

            switch (faces[i]->category)
            {
            case 1: faces[i]->model_axis = "I+";
                break;
            case 2: faces[i]->model_axis = "I-";
                break;
            case 3: faces[i]->model_axis = "J+";
                break;
            case 4: faces[i]->model_axis = "J-";
                break;
            case 5: faces[i]->model_axis = "K+";
                break;
            case 6: faces[i]->model_axis = "K-";
                break;
			default:  if (verbosity > 1) std::cout << " WARNING: Unknown face appeared for block " << faces[faceindex]->connected_blocks[0] << '\n';
                break;
            }
            vec_face_nodes.clear();
            for (size_t j = 0; j < faces[i]->NodeIndex.size(); j++)
                vec_face_nodes.push_back(m_msh->nod_vector[faces[i]->NodeIndex[j]]);

            if (faces[i]->nfnodes == 4)
            {
                if (faces[i]->CheckIfPointFormPlane(vec_face_nodes[0], vec_face_nodes[1], vec_face_nodes[2],
                    vec_face_nodes[3]) == true)
                { //Four nodes form a plane
                    faces[i]->CreateFace(vec_face_nodes[0], vec_face_nodes[1], vec_face_nodes[2], vec_face_nodes[3]);
                    for (unsigned long j = 0; j < vec_face_nodes.size(); j++)
                    {
                        //Connect nodes with faces
                        m_msh->nod_vector[vec_face_nodes[j]->GetIndex()]->connected_faces.push_back(faces[i]->index);
                        //Calculate distance between node and gravity centre of the face and store it in a vector
                        //m_msh->nod_vector[vec_face_nodes[j]->GetIndex()]->distance_to_connected_faces.push_back(
                        //this->CalculateDistanceBetween2Points(vec_face_nodes[j]->getData(), 
                        //m_face->GetFaceGravityCentre()));
                        m_msh->nod_vector[vec_face_nodes[j]->GetIndex()]->
                            distance_to_connected_faces.push_back(sqrt(MathLib::sqrDist(vec_face_nodes[j]->getData(),
                            faces[i]->GetFaceGravityCentre())));
                    }
                }
                else
                { //Four nodes cannot form a plane needs to splitted
                    threenodesflag = true;
                    //Face 1:node0,node1,node2
                    m_face = new CFaces(static_cast<int>(this->Phases.size()));

                    m_face->index = faces[i]->index;
                    faces[i]->splittedfaceindex.push_back(faces.size());
                    m_face->model_axis = faces[i]->model_axis;

                    vec_face_nodes.clear();
                    vec_face_nodes.push_back(m_msh->nod_vector[faces[i]->NodeIndex[0]]);
                    vec_face_nodes.push_back(m_msh->nod_vector[faces[i]->NodeIndex[1]]);
                    vec_face_nodes.push_back(m_msh->nod_vector[faces[i]->NodeIndex[2]]);

                    m_face->CreateFace(vec_face_nodes[0], vec_face_nodes[1], vec_face_nodes[2],
                        vec_face_nodes[2], threenodesflag);
                    for (unsigned long j = 0; j < vec_face_nodes.size(); j++)
                    {
                        m_msh->nod_vector[vec_face_nodes[j]->GetIndex()]->connected_faces.push_back(this->faces.size());
                        m_msh->nod_vector[vec_face_nodes[j]->GetIndex()]->
                            distance_to_connected_faces.push_back(sqrt(MathLib::sqrDist(vec_face_nodes[j]->getData(),
                            m_face->GetFaceGravityCentre())));
                    }
                    faces.push_back(m_face);

                    //Face 2:node3,node2,node1
                    m_face = new CFaces(static_cast<int>(this->Phases.size()));

                    m_face->index = faces[i]->index;
                    faces[i]->splittedfaceindex.push_back(faces.size());
                    m_face->model_axis = faces[i]->model_axis;

                    vec_face_nodes.clear();
                    vec_face_nodes.push_back(m_msh->nod_vector[faces[i]->NodeIndex[3]]);
                    vec_face_nodes.push_back(m_msh->nod_vector[faces[i]->NodeIndex[2]]);
                    vec_face_nodes.push_back(m_msh->nod_vector[faces[i]->NodeIndex[1]]);

                    m_face->CreateFace(vec_face_nodes[0], vec_face_nodes[1], vec_face_nodes[2],
                        vec_face_nodes[2], threenodesflag);
                    for (unsigned long j = 0; j < vec_face_nodes.size(); j++)
                    {
                        m_msh->nod_vector[vec_face_nodes[j]->GetIndex()]->connected_faces.push_back(this->faces.size());
                        m_msh->nod_vector[vec_face_nodes[j]->GetIndex()]->
                            distance_to_connected_faces.push_back(sqrt(MathLib::sqrDist(vec_face_nodes[j]->getData(),
                            m_face->GetFaceGravityCentre())));
                    }
                    faces.push_back(m_face);
                    // WTP debug
                    //double test = faces[i]->splittedfaceindex[0];
                    //double test2 = faces[i]->splittedfaceindex[1];
                    // Add up the area for non-splitted face
                    faces[i]->face_area = faces[faces[i]->splittedfaceindex[0]]->face_area +
                        faces[faces[i]->splittedfaceindex[1]]->face_area;
                }
            }
            else if (vec_face_nodes.size() == 3)
            {
                threenodesflag = true;
                faces[i]->CreateFace(vec_face_nodes[0], vec_face_nodes[1], vec_face_nodes[2],
                    vec_face_nodes[2], threenodesflag);
                for (unsigned long j = 0; j < vec_face_nodes.size(); j++)
                {
                    m_msh->nod_vector[vec_face_nodes[j]->GetIndex()]->connected_faces.push_back(faces[i]->index);
                    m_msh->nod_vector[vec_face_nodes[j]->GetIndex()]->
                        distance_to_connected_faces.push_back(sqrt(MathLib::sqrDist(vec_face_nodes[j]->getData(),
                        faces[i]->GetFaceGravityCentre())));
                }
            }
        }
        ////data output
        //vector <string> vec_string;
        //ostringstream temp;
        //tempstring = " $FACES";
        //vec_string.push_back(tempstring);
        //// Loop over all elements
        //for (long i = 0; i < this->FlowFaceNum; i++)
        //{
        //	temp.precision(12);
        //	temp.str("");
        //	temp.clear();
        //	temp << i;
        //	tempstring = temp.str();
        //	temp.str("");
        //	temp.clear();
        //	temp << this->faces[i]->nfnodes;
        //	tempstring += "  " + temp.str();
        //	temp.str("");
        //	temp.clear();
        //	temp << this->faces[i]->category;
        //	tempstring += "  " + temp.str();
        //	//FaceNodeIndex
        //	for (size_t j = 0; j < this->faces[i]->nfnodes; j++){
        //		temp.str(""); temp.clear(); temp << this->faces[i]->NodeIndex[j]; tempstring += " " + temp.str();
        //	}
        //	//FaceElementIndex
        //	temp.str(""); temp.clear(); temp << this->faces[i]->connected_blocks.size(); tempstring += " " + temp.str();
        //	for (size_t j = 0; j < this->faces[i]->connected_blocks.size(); j++){
        //		temp.str(""); temp.clear(); temp << this->faces[i]->connected_blocks[j]; tempstring += " " + temp.str();
        //	}
        //	vec_string.push_back(tempstring);
        //} // 
        //// "#STOP": OGS Keyword of ending for a section
        //tempstring = "#STOP";
        //vec_string.push_back(tempstring);

        //// Test Output
        //Filename = Projectname + ".bfaces";
        //ofstream aus;
        //aus.open(Filename.data(), ios::out);
        //for (unsigned int i = 0; i < vec_string.size(); i++)
        //	aus << vec_string[i] << "\n";
        //aus.close();

        CreateFaces_called = true;
        finish = clock();
        time = (static_cast<double>(finish)-static_cast<double>(start)) / CLOCKS_PER_SEC;
		if (verbosity > 2)
			std::cout << "                        Time: " << time << " seconds." << "\n";

    }
	return true;

}
/*-------------------------------------------------------------------------
GeoSys - Function: FaceGenerator
Task: For each ECLIPSE Block, using NodeIndex construct faces
Return: Nothing
Programming: 09/2014 BW
Modification:
-------------------------------------------------------------------------*/
void CECLIPSEData::FaceGenerator(long eclcellindex, int facenodeindex[], int category)
{
	CFaces *m_Face = NULL;
	vector <MeshLib::CNode*> vec_face_nodes;
	//CFEMesh* m_msh = fem_msh_vector[0];
	int samenodecounts = 0;
	long tempblockindex = -9999;

	//temporarily saving the nodes for constructing the face
	int facerecord[] = { 0, 0, 0, 0 };
	for (int i = 0; i < 4; i++)
	{
		bool samenodeflag = false;
		for (unsigned int j = i + 1; j < 4; j++)
		{
			if (this->eclgrid[eclcellindex]->Nodeindex[facenodeindex[j]] ==
				this->eclgrid[eclcellindex]->Nodeindex[facenodeindex[i]])
			{
				samenodeflag = true;
				samenodecounts += 1;
				break;
			}
		}
		if (samenodeflag == false)
			facerecord[i - samenodecounts] = facenodeindex[i];
	}

	// category group
	int categorymatch = 0;
	if (category == 1)	categorymatch = 2;
	if (category == 2)	categorymatch = 1;
	if (category == 3)	categorymatch = 4;
	if (category == 4)	categorymatch = 3;
	if (category == 5)	categorymatch = 6;
	if (category == 6)	categorymatch = 5;

	if (samenodecounts < 2)
	{
		bool existedface = false;

		//I+ Face definitely needs to build
		if (category == 1)
		{
			//if (this->eclgrid[eclcellindex]->column == this->columns){
			m_Face = new CFaces(static_cast<int>(this->Phases.size()));
			m_Face->nfnodes = 4 - samenodecounts;  // Set Facenodes Number: 4 or 3

			m_Face->NodeIndex.resize(m_Face->nfnodes);	  // Set Facenodes Number: 4 or 3
			m_Face->category = category;  // Set Nodes Category

			for (int k = 0; k < m_Face->nfnodes; k++)
			{	 // Record Nodes Number
				m_Face->NodeIndex[k] = this->eclgrid[eclcellindex]->Nodeindex[facerecord[k]];
			}

			m_Face->connected_blocks.push_back(eclcellindex); //Record Block Index
			m_Face->index = long(this->faces.size());
			faces.push_back(m_Face);
			this->eclgrid[eclcellindex]->FaceIndex[category - 1] = long(this->faces.size()) - 1;	 //Record Face index
			return;
			//}
		}

		//I- Face
		if (category == 2)
		{
			//At first column, definitely needs to build
			if (this->eclgrid[eclcellindex]->column == 1)
			{
				m_Face = new CFaces(static_cast<int>(this->Phases.size()));
				m_Face->nfnodes = 4 - samenodecounts;  // Set Facenodes Number: 4 or 3

				m_Face->NodeIndex.resize(m_Face->nfnodes);	  // Set Facenodes Number: 4 or 3
				m_Face->category = category;  // Set Nodes Category

				for (int k = 0; k < m_Face->nfnodes; k++)
				{	 // Record Nodes Number
					m_Face->NodeIndex[k] = this->eclgrid[eclcellindex]->Nodeindex[facerecord[k]];
				}
				m_Face->connected_blocks.push_back(eclcellindex); //Record Block Index
				m_Face->index = long(this->faces.size());
				faces.push_back(m_Face);
				this->eclgrid[eclcellindex]->FaceIndex[category - 1] = long(faces.size()) - 1;	 //Record Face index
				return;
			}

			//Not first column, first find the active block beside it, figure the face is identical or not

			tempblockindex = this->IndexArray[this->eclgrid[eclcellindex]->column - 2][this->eclgrid[eclcellindex]->row - 1][this->eclgrid[eclcellindex]->layer - 1];
			int ii = 3;
			while (tempblockindex == -1 && (this->eclgrid[eclcellindex]->column - ii) >= 0)
			{
				tempblockindex = this->IndexArray[this->eclgrid[eclcellindex]->column - ii][this->eclgrid[eclcellindex]->row - 1][this->eclgrid[eclcellindex]->layer - 1];
				ii++;
			}
		}

		//J+ Face definitely needs to build
		if (category == 3)
		{
			m_Face = new CFaces(static_cast<int>(this->Phases.size()));
			m_Face->nfnodes = 4 - samenodecounts;  // Set Facenodes Number: 4 or 3

			m_Face->NodeIndex.resize(m_Face->nfnodes);	  // Set Facenodes Number: 4 or 3
			m_Face->category = category;  // Set Nodes Category

			for (int k = 0; k < m_Face->nfnodes; k++)
			{	 // Record Nodes Number
				m_Face->NodeIndex[k] = this->eclgrid[eclcellindex]->Nodeindex[facerecord[k]];
			}
			m_Face->connected_blocks.push_back(eclcellindex); //Record Block Index
			m_Face->index = long(this->faces.size());
			faces.push_back(m_Face);
			this->eclgrid[eclcellindex]->FaceIndex[category - 1] = long(faces.size()) - 1;	 //Record Face index
			return;
			/*}*/
		}

		//j- Face
		if (category == 4)
		{
			//At first row, definitely needs to build
			if (this->eclgrid[eclcellindex]->row == 1)
			{
				m_Face = new CFaces(static_cast<int>(this->Phases.size()));
				m_Face->nfnodes = 4 - samenodecounts;  // Set Facenodes Number: 4 or 3

				m_Face->NodeIndex.resize(m_Face->nfnodes);	  // Set Facenodes Number: 4 or 3
				m_Face->category = category;  // Set Nodes Category

				for (int k = 0; k < m_Face->nfnodes; k++)
				{	 // Record Nodes Number
					m_Face->NodeIndex[k] = this->eclgrid[eclcellindex]->Nodeindex[facerecord[k]];
				}
				m_Face->connected_blocks.push_back(eclcellindex); //Record Block Index
				m_Face->index = long(this->faces.size());
				faces.push_back(m_Face);
				this->eclgrid[eclcellindex]->FaceIndex[category - 1] = long(faces.size()) - 1;	 //Record Face index
				return;
			}

			//Not first row, first find the active block along the same column, figure the face is identical or not
			tempblockindex = this->IndexArray[this->eclgrid[eclcellindex]->column - 1][this->eclgrid[eclcellindex]->row - 2][this->eclgrid[eclcellindex]->layer - 1];
			int ii = 3;
			//cout << this->eclgrid[eclcellindex]->row - ii << '\n';
			while (tempblockindex == -1 && (this->eclgrid[eclcellindex]->row - ii) >= 0)
			{
				tempblockindex = this->IndexArray[this->eclgrid[eclcellindex]->column - 1][this->eclgrid[eclcellindex]->row - ii][this->eclgrid[eclcellindex]->layer - 1];
				ii++;
			}

		}

		//K+ Face at first layer, definitely needs to build
		if (category == 5)
		{
			m_Face = new CFaces(static_cast<int>(this->Phases.size()));
			m_Face->nfnodes = 4 - samenodecounts;  // Set Facenodes Number: 4 or 3

			m_Face->NodeIndex.resize(m_Face->nfnodes);	  // Set Facenodes Number: 4 or 3
			m_Face->category = category;  // Set Nodes Category

			for (int k = 0; k < m_Face->nfnodes; k++)
			{	 // Record Nodes Number
				m_Face->NodeIndex[k] = this->eclgrid[eclcellindex]->Nodeindex[facerecord[k]];
			}
			m_Face->connected_blocks.push_back(eclcellindex); //Record Block Index
			m_Face->index = long(this->faces.size());
			faces.push_back(m_Face);
			this->eclgrid[eclcellindex]->FaceIndex[category - 1] = long(faces.size()) - 1;	 //Record Face index
			return;
		}

		//K- Face
		if (category == 6)
		{
			//At first layer, definitely needs to build
			if (this->eclgrid[eclcellindex]->layer == 1)
			{
				m_Face = new CFaces(static_cast<int>(this->Phases.size()));
				m_Face->nfnodes = 4 - samenodecounts;  // Set Facenodes Number: 4 or 3

				m_Face->NodeIndex.resize(m_Face->nfnodes);	  // Set Facenodes Number: 4 or 3
				m_Face->category = category;  // Set Nodes Category

				for (int k = 0; k < m_Face->nfnodes; k++)
				{	 // Record Nodes Number
					m_Face->NodeIndex[k] = this->eclgrid[eclcellindex]->Nodeindex[facerecord[k]];
				}
				m_Face->connected_blocks.push_back(eclcellindex); //Record Block Index
				m_Face->index = long(this->faces.size());
				faces.push_back(m_Face);
				this->eclgrid[eclcellindex]->FaceIndex[category - 1] = long(faces.size()) - 1;	 //Record Face index
				return;
			}

			//Not first layer, first find the active block along the same column, figure the face is identical or not
			tempblockindex = this->IndexArray[this->eclgrid[eclcellindex]->column - 1][this->eclgrid[eclcellindex]->row - 1][this->eclgrid[eclcellindex]->layer - 2];
			int ii = 3;
			while (tempblockindex == -1 && (this->eclgrid[eclcellindex]->layer - ii) >= 0)
			{
				tempblockindex = this->IndexArray[this->eclgrid[eclcellindex]->column - 1][this->eclgrid[eclcellindex]->row - 1][this->eclgrid[eclcellindex]->layer - ii];
				ii++;
			}
		}

		//double Check!!
		if (tempblockindex != -1)
		{
			long tempindex = this->eclgrid[tempblockindex]->FaceIndex[categorymatch - 1];
			if (tempindex != -1)
			{//face existed
				size_t nodesize = faces[tempindex]->NodeIndex.size();
				if (faces[tempindex]->NodeIndex[0] == this->eclgrid[eclcellindex]->Nodeindex[facerecord[0]])
				{
					if (faces[tempindex]->NodeIndex[1] == this->eclgrid[eclcellindex]->Nodeindex[facerecord[1]])
					{
						if (faces[tempindex]->NodeIndex[2] == this->eclgrid[eclcellindex]->Nodeindex[facerecord[2]])
						{
							if (nodesize == 4)
							{
								if (faces[tempindex]->NodeIndex[3] == this->eclgrid[eclcellindex]->Nodeindex[facerecord[3]])
								{
									this->eclgrid[eclcellindex]->FaceIndex[category - 1] = tempindex;
									faces[tempindex]->connected_blocks.push_back(eclcellindex);
									existedface = true;
									return;
								}
							}
							else if (nodesize == 3)
							{
								this->eclgrid[eclcellindex]->FaceIndex[category - 1] = tempindex;
								faces[tempindex]->connected_blocks.push_back(eclcellindex);
								existedface = true;
								return;
							}
						}
					}
				}
			}
		}

		//std::cout << "Warning: Block " << eclcellindex << " Face " << category << " can not create following I,J,K, it needs iteration to create." << '\n';
		{
			unsigned long m;
			for (m = 0; m < faces.size() && existedface == false; m++)
			{
				size_t nodesize = faces[m]->NodeIndex.size();
				int category_local = faces[m]->category;

				if (static_cast<int>(nodesize) == 4 - samenodecounts)
				{
					if (categorymatch == category_local)
					{
						if (faces[m]->NodeIndex[0] == this->eclgrid[eclcellindex]->Nodeindex[facerecord[0]])
						{
							if (faces[m]->NodeIndex[1] == this->eclgrid[eclcellindex]->Nodeindex[facerecord[1]])
							{
								if (faces[m]->NodeIndex[2] == this->eclgrid[eclcellindex]->Nodeindex[facerecord[2]])
								{
									if (nodesize == 4)
									{
										if (faces[m]->NodeIndex[3] == this->eclgrid[eclcellindex]->Nodeindex[facerecord[3]])
										{
											this->eclgrid[eclcellindex]->FaceIndex[category - 1] = m;
											faces[m]->connected_blocks.push_back(eclcellindex);
											existedface = true;
											return;
										}
									}
									else if (nodesize == 3)
									{
										this->eclgrid[eclcellindex]->FaceIndex[category - 1] = m;
										faces[m]->connected_blocks.push_back(eclcellindex);
										existedface = true;
										return;
									}
								}
							}
						}
					}
				}
				if (existedface == true)
					break;
			}
		}// end omp

		if (existedface == false)
		{
			m_Face = new CFaces(static_cast<int>(this->Phases.size()));
			m_Face->nfnodes = 4 - samenodecounts;  // Set Facenodes Number: 4 or 3

			m_Face->NodeIndex.resize(m_Face->nfnodes);	  // Set Facenodes Number: 4 or 3
			m_Face->category = category;  // Set Nodes Category

			for (int k = 0; k < m_Face->nfnodes; k++)
			{	 // Record Nodes Number
				m_Face->NodeIndex[k] = this->eclgrid[eclcellindex]->Nodeindex[facerecord[k]];
			}
			m_Face->connected_blocks.push_back(eclcellindex); //Record Block Index
			m_Face->index = long(this->faces.size());
			faces.push_back(m_Face);
			this->eclgrid[eclcellindex]->FaceIndex[category - 1] = long(faces.size()) - 1;	 //Record Face index
		}
	}
	else
	{
		this->eclgrid[eclcellindex]->FaceIndex[category - 1] = -1;	   //two nodes will not have a face
	}
};

/*-------------------------------------------------------------------------
   GeoSys - Function: CreateFacesAtElements
   Task: For each element, find all faces and store them in connected_faces
   Return: bool value if there occured an error
   Programming: 09/2009 SB
   Modification: 04/2014 WB: for prism cell collapsed faces
   -------------------------------------------------------------------------*/
bool CECLIPSEData::ConnectFacesToElements(void)
{
	//CECLIPSEBlock*m_block=NULL;
	CFaces* m_face = NULL;
	long ind_block, ind_face;
	vector <long> faces_at_block;
	clock_t start, finish;
	double time;

	start = clock();

	if (verbosity > 2)
		std::cout << "        ConnectFacesToElements() ";
    if (this->ConnectFacesToElements_called == true)
    {
		if (verbosity > 2)
			std::cout << " - Function already called once, skip this time" << "\n";
    }
    else
    {
        // go through all faces and add them to the eclipse blocks
        for (unsigned long i = 0; i < this->faces.size(); i++)
        {
            m_face = this->faces[i];
            ind_face = m_face->index;
            for (unsigned long j = 0; j < m_face->connected_blocks.size(); j++)
            {
                ind_block = m_face->connected_blocks[j];
                this->eclgrid[ind_block]->connected_faces.push_back(ind_face);
            }
        }
        // Test output
        //for(unsigned long i=0;i<eclgrid.size();i++){
        //    std::cout << "    Faces at block " << eclgrid[i]->index << " :  ";
        //    for(int j=0;j<eclgrid[i]->connected_faces.size();j++) std::cout << eclgrid[i]->connected_faces[j] << ", ";
        //    std::cout << "\n";
        //}

        for (long i = 0; i < long(eclgrid.size()); i++)
        {
            // Test output: all faces at one block
			if (eclgrid[i]->connected_faces.size() != 6  && verbosity > 1)
                std::cout << "  WARNING: Block " << i << " has " << eclgrid[i]->connected_faces.size() << " faces." << "\n";

            // initialize vector faces_at_block
            faces_at_block.clear();
            for (unsigned long j = 0; j < 6; j++)
                faces_at_block.push_back(-1);

            for (unsigned long j = 0; j < eclgrid[i]->connected_faces.size(); j++)
            {
                //std::cout << eclgrid[i]->connected_faces[j] << ", ";
                // Order faces +i, -i, +j. -j, +k , -k
                m_face = faces[eclgrid[i]->connected_faces[j]]; // get face
                if (m_face->model_axis == "I+")
                {
                    if (m_face->connected_blocks[0] == i)
                        faces_at_block[0] = m_face->index;  // this is +i for this element
                    else
                        faces_at_block[1] = m_face->index;  // this is -i for this element
                }
                if (m_face->model_axis == "I-")
                {
                    if (m_face->connected_blocks[0] == i)
                        faces_at_block[1] = m_face->index;  // for this element it is i-
                    else
                        faces_at_block[0] = m_face->index;  // for this element it is i+
                }
                if (m_face->model_axis == "J+")
                {
                    if (m_face->connected_blocks[0] == i)
                        faces_at_block[2] = m_face->index;
                    else
                        faces_at_block[3] = m_face->index;
                }
                if (m_face->model_axis == "J-")
                {
                    if (m_face->connected_blocks[0] == i)
                        faces_at_block[3] = m_face->index;
                    else
                        faces_at_block[2] = m_face->index;
                }
                if (m_face->model_axis == "K+")
                {
                    if (m_face->connected_blocks[0] == i)
                        faces_at_block[4] = m_face->index;
                    else
                        faces_at_block[5] = m_face->index;
                }
                if (m_face->model_axis == "K-")
                {
                    if (m_face->connected_blocks[0] == i)
                        faces_at_block[5] = m_face->index;
                    else
                        faces_at_block[4] = m_face->index;
                }
            }
            // copy sorted faces to ecl block
            long face_index_adjust = 0;
            for (unsigned long j = 0; j < eclgrid[i]->connected_faces.size(); j++)
            {
                // loop over index of collapsed faces:WB
                if (faces_at_block[j] == -1)
                    face_index_adjust += 1;
                eclgrid[i]->connected_faces[j] = faces_at_block[j + face_index_adjust];

                //// Test output
                //std::cout << "    Faces at block " << eclgrid[i]->index << " :  ";
                //for(j=0;j<eclgrid[i]->connected_faces.size();j++) std::cout << eclgrid[i]->connected_faces[j] << ", ";
                //std::cout << "\n";
            }
        }

        ConnectFacesToElements_called = true;
        finish = clock();
        time = (static_cast<double>(finish)-static_cast<double>(start)) / CLOCKS_PER_SEC;
		if (verbosity > 2)
			std::cout << "             Time: " << time << " seconds." << "\n";

    }
	return true;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: GetFlowForFaces
   Task: Get's the flow from the corresponding elements
   Return: true or false
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
bool CECLIPSEData::GetFlowForFaces(int phase_index)
{
	//CFEMesh* m_msh = fem_msh_vector[0];
	CFaces* m_face = NULL;
	long number_elements;
	long element_index;
	long variable_index;
	//WW long time_index;
	double vz = 1.0; // Sign
	clock_t start, finish;
	double time;

	start = clock();
	if (verbosity > 2)
		std::cout << "        GetFlowForFaces()";

	// WTP Determine current phase
	int phase_flag = 0;
	if (this->Phases[phase_index] == "WATER")
		phase_flag = 1;
	else if (Phases[phase_index] == "GAS")
		phase_flag = 2;
	else if (Phases[phase_index] == "OIL")
		phase_flag = 3;
	else
	{
		std::cout << " ERROR: " << Phases[phase_index] << " is not yet considered in GetFlowForFaces" << "\n";
		std::cout << flush;
		//system("Pause");
		exit(0);
	}


	//for (long i = 0; i < long(this->faces.size()); i++)
	for (long i = 0; i < this->FlowFaceNum; i++)
	{
		m_face = this->faces[i];
		//set initial values
		m_face->phases[phase_index]->q_norm = 0.0;

		//number_elements = long(m_face->connected_elements.size());
		number_elements = long(m_face->connected_blocks.size());

		if ((number_elements < 1) || (number_elements > 2))//Faces could only connected to 2 elements
		{
			std::cout << " ERROR: Problem connecting faces and elements" << "\n";
            std::cout << "        number of connected elements to face " << i << " is < 1 or > 2" << "\n";
            std::cout << flush;
			//system("Pause");
			exit(0);
		}

		// right face (x)
		if (m_face->model_axis == "I+")
		{
			//element_index = m_face->connected_elements[0];
			element_index = m_face->connected_blocks[0];
			//WW		time_index = 0;
			if (this->RadialModellJpos == true)
				variable_index = GetVariableIndex(GetECLFlowVariableName(phase_flag, 2));
			else
				variable_index = GetVariableIndex(GetECLFlowVariableName(phase_flag, 1));

			if (variable_index < 0)
			{
				std::cout << " ERROR: There are no variables!" << "\n";
                std::cout << flush;
				//system("Pause");
				exit(0);
			}
			//calculate q [m3/m2.s] for Geosys from Q [m3/d] from Eclipse
			m_face->phases[phase_index]->q_norm =
				(this->Data[element_index][variable_index] /
				m_face->face_area) / 86400.0;
			//BW:if face is splitted into two face, due to distortion
			if (m_face->splittedfaceindex.size() > 0)
				for (size_t tempi = 0; tempi < m_face->splittedfaceindex.size(); tempi++)
				{
					this->faces[m_face->splittedfaceindex[tempi]]->phases[phase_index]->q_norm = m_face->phases[phase_index]->q_norm;
					this->faces[m_face->splittedfaceindex[tempi]]->Calculate_components_of_a_vector(0, phase_index, Radialmodell);
				}				
			else
				m_face->Calculate_components_of_a_vector(0, phase_index, Radialmodell);
			//std::cout << " Flow through face " << m_face->index << " element " << element_index << " :  " 
			//<< this->Data[element_index][time_index][variable_index] << " ,  " << m_face->q_norm << "\n";
		}
		// upper face (y)
		if (m_face->model_axis == "J+")
		{
			element_index = m_face->connected_blocks[0];
			if (this->RadialModellJpos == true)
				variable_index = GetVariableIndex(GetECLFlowVariableName(phase_flag, 1));
			else
				variable_index = GetVariableIndex(GetECLFlowVariableName(phase_flag, 2));

			if (variable_index < 0)
			{
				std::cout << " ERROR: There are no variables!" << "\n";
                std::cout << flush;
				//system("Pause");
				exit(0);
			}
			//calculate q [m3/m2.s] for Geosys from Q [m3/d] from Eclipse
			m_face->phases[phase_index]->q_norm =
				(this->Data[element_index][variable_index] /
				m_face->face_area) / 86400.0;
			//BW:if face is splitted into two face, due to distortion
			if (m_face->splittedfaceindex.size() > 0)
				for (size_t tempi = 0; tempi < m_face->splittedfaceindex.size(); tempi++)
				{
					this->faces[m_face->splittedfaceindex[tempi]]->phases[phase_index]->q_norm = m_face->phases[phase_index]->q_norm;
					this->faces[m_face->splittedfaceindex[tempi]]->Calculate_components_of_a_vector(0, phase_index, Radialmodell);
				}
			else
				m_face->Calculate_components_of_a_vector(0, phase_index, Radialmodell);
			//std::cout << " Flow through face " << m_face->index <<  " element " << element_index << " :  " << 
			//this->Data[element_index][time_index][variable_index] << " ,  " 
			//<< m_face->q_norm << "\n";
		}
		// bottom face (z) (Eclipse has the opposite z direction -> k = +
		if (m_face->model_axis == "K+")
		{
			//element_index = m_face->connected_elements[0];
			element_index = m_face->connected_blocks[0];
			variable_index = GetVariableIndex(GetECLFlowVariableName(phase_flag, 3)); // K+ direction

			//WW time_index = 0;
			if (variable_index < 0)
			{
				std::cout << " ERROR: There are no variables!" << "\n";
                std::cout << flush;
				//system("Pause");
				exit(0);
			}
			//calculate q [m3/m2.s] for Geosys from Q [m3/d] from Eclipse
			m_face->phases[phase_index]->q_norm =
				(this->Data[element_index][variable_index] /
				m_face->face_area) / 86400.0;
			//BW:if face is splitted into two face, due to distortion
			if (m_face->splittedfaceindex.size() > 0)
				for (size_t tempi = 0; tempi < m_face->splittedfaceindex.size(); tempi++)
				{
					this->faces[m_face->splittedfaceindex[tempi]]->phases[phase_index]->q_norm = m_face->phases[phase_index]->q_norm;
					this->faces[m_face->splittedfaceindex[tempi]]->Calculate_components_of_a_vector(0, phase_index, Radialmodell);
				}
			else
				m_face->Calculate_components_of_a_vector(0, phase_index, Radialmodell);
			//std::cout << " Flow through face " << m_face->index <<  " element " << element_index << " :  " 
			//<< this->Data[element_index][time_index][variable_index] << " ,  " << m_face->q_norm << "\n";
		}
		// left, lower and top faces -> now flow over these faces
		if ((m_face->model_axis == "I-") || (m_face->model_axis == "J-") ||
			(m_face->model_axis == "K-"))
		{
			m_face->phases[phase_index]->q_norm = 0;
			//BW:if face is splitted into two face, due to distortion
			if (m_face->splittedfaceindex.size() > 0)
				for (size_t tempi = 0; tempi < m_face->splittedfaceindex.size(); tempi++)
				{
					this->faces[m_face->splittedfaceindex[tempi]]->phases[phase_index]->q_norm = m_face->phases[phase_index]->q_norm;
					this->faces[m_face->splittedfaceindex[tempi]]->Calculate_components_of_a_vector(0, phase_index, Radialmodell);
				}
				else
					m_face->Calculate_components_of_a_vector(0, phase_index, Radialmodell);
		}

		//Get additional flow for border cells that contain boundary conditions
		//if (m_face->connected_elements.size() == 1)
		if (m_face->connected_blocks.size() == 1)
		{
			//long ele_index = m_face->connected_elements[0];
			long ele_index = m_face->connected_blocks[0];
			//check if there is one BC in the element
			if (this->eclgrid[ele_index]->ConnectedBoundaryCondition.size() > 1)
			{
				std::cout <<
					" ERROR: There is more than 1 boundary condition assigned to the cell " <<
					ele_index << "\n";
                std::cout << flush;
				//system("Pause");
				exit(0);
			}
			if (this->eclgrid[ele_index]->ConnectedBoundaryCondition.size() == 1)
			{
				long bc_index =
					this->eclgrid[ele_index]->ConnectedBoundaryCondition[0];
				//check if the boundary inflow fits to the current face
				if (m_face->model_axis == this->BC[bc_index]->boundary_position)
				{
					//check if q is 0
					if (m_face->phases[phase_index]->q_norm == 0)
					{
						//calculate q [m3/m2.s] for Geosys from Q [m3/d] from Eclipse
						//cout << " BC Flow through face " << m_face->index << " :  " << this->BC[bc_index]->value[0] << "\n";
						vz = 1.0;
						//if(this->BC[bc_index]->value[0] < 0){ // outflow out of model area
						//The values in I+, J+ and K+ direction have to be counted with the opposit direction
						// because the flow from BC is given positive if fluid flows from the BC into the model
						if (m_face->model_axis == "I+")
							vz = -1.0;
						if (m_face->model_axis == "J+")
							vz = -1.0;
						if (m_face->model_axis == "K+")
							vz = -1.0;
						//}
						m_face->phases[phase_index]->q_norm =
							(this->BC[bc_index]->value[0] * vz /
							m_face->face_area) / 86400.0;
						//BW:if face is splitted into two face, due to distortion
						if (m_face->splittedfaceindex.size() > 0)
							for (size_t tempi = 0; tempi < m_face->splittedfaceindex.size(); tempi++)
							{
								this->faces[m_face->splittedfaceindex[tempi]]->phases[phase_index]->q_norm = m_face->phases[phase_index]->q_norm;
								this->faces[m_face->splittedfaceindex[tempi]]->Calculate_components_of_a_vector(0, phase_index, Radialmodell);
							}
						else
							m_face->Calculate_components_of_a_vector(0, phase_index, Radialmodell);
						//std::cout << " Flow through face " << m_face->index <<  " element " << element_index << " :  " 
						//<< this->BC[bc_index]->value[0] << " ,  " << m_face->q_norm << "\n";
					}
					else
					{
						std::cout <<
							" ERROR: There is already a flow assigned to the boundary face, which shouldn't happen "
							<< m_face->index << "\n";
                        std::cout << flush;
						//system("Pause");
						exit(0);
					}
				}
			}
		}
		//std::cout << " Face " << m_face->index << ", " << m_face->connected_blocks[0] ;
		//if(m_face->connected_blocks.size() > 1) std::cout << ", " << m_face->connected_blocks[1];
		//std::cout << ":           " << m_face->q_norm << ":   " << m_face->q[0] << ", " << m_face->q[1] << ", " << m_face->q[2] << "\n";
	}

	finish = clock();
	time = (static_cast<double>(finish)-static_cast<double>(start)) / CLOCKS_PER_SEC;
	if (verbosity > 2)
	{
		std::cout << "                     Time: " << time << " seconds." << "\n";
		std::cout << flush;
	}
	return true;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: GetVelForFaces
   Task: Get's the velocities from the corresponding elements
   derived from GetFlowForFaces()
   Return: true or false
   Programming: 09/2009 SB
   Modification:
   -------------------------------------------------------------------------*/
bool CECLIPSEData::GetVelForFaces(void)
{
	//CFEMesh* m_msh = fem_msh_vector[0];
	CFaces* m_face = NULL;
	long number_elements;
	long element_index;
	long variable_index;
	//WW long time_index;

	for (unsigned long i = 0; i < this->faces.size(); i++)
	{
		m_face = this->faces[i];
		number_elements = long(m_face->connected_blocks.size());

		if ((number_elements < 1) || (number_elements > 2))
			return false;

		// right face (x)
		if (m_face->model_axis == "I+")
		{
			element_index = m_face->connected_blocks[0];
			//WW		time_index = 0;
			variable_index = this->GetVariableIndex("VELWATI+");
			if (variable_index < 0)
				return false;
			m_face->v_norm = this->Data[element_index][variable_index];
			//m_face->Calculate_components_of_a_vector(1);
		}
		// upper face (y)
		if (m_face->model_axis == "J+")
		{
			element_index = m_face->connected_blocks[0];
			variable_index = this->GetVariableIndex("VELWATJ+");
			//WW		time_index = 0;
			if (variable_index < 0)
				return false;
			m_face->v_norm = this->Data[element_index][variable_index];
			//m_face->Calculate_components_of_a_vector(1);
		}
		// bottom face (z) (Eclipse has the opposite z direction -> k = +
		if (m_face->model_axis == "K+")
		{
			element_index = m_face->connected_blocks[0];
			variable_index = this->GetVariableIndex("VELWATK+");
			//WW		time_index = 0;
			if (variable_index < 0)
				return false;
			m_face->v_norm = this->Data[element_index][variable_index];
			//m_face->Calculate_components_of_a_vector(1);
		}
		// left, lower and top faces -> now flow over these faces
		if ((m_face->model_axis == "I-") || (m_face->model_axis == "J-") ||
			(m_face->model_axis == "K-"))
			m_face->v_norm = 0;
		//m_face->Calculate_components_of_a_vector(1);
		if (verbosity > 1)
		{
			std::cout << " Face " << m_face->index << ", " << m_face->connected_blocks[0];
			if (m_face->connected_blocks.size() > 1)
				std::cout << ", " << m_face->connected_blocks[1];
		}

		//if(fabs(m_face->v_norm) > MKleinsteZahl) std::cout << " q/v: " << m_face->q_norm/m_face->face_area/0.5 /m_face->v_norm ;
		//std::cout << "\n";
	}
	return true;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: CalcBlockBudegt
   Task: Summ all budget terms from the faces for each block - corresponds to a mass balance in steady state
   This function is for checking purposes only
   Return: true or false
   Programming: 09/2009 SB
   Modification:
   -------------------------------------------------------------------------*/
bool CECLIPSEData::CalcBlockBudget(int phase_index)
{
	CECLIPSEBlock* m_block = NULL;
	CFaces* m_face = NULL;
	double flow, flow_face, flow_max;
	double max_error = 0;

	// summ budget flow terms at feces for all blocks
	for (unsigned i = 0; i < this->eclgrid.size(); i++)
	{
		m_block = this->eclgrid[i];
		flow = flow_max = 0.0;
		//std::cout << "  Budget for block " << i << ": " << "\n" ;
		for (unsigned int j = 0; j < m_block->connected_faces.size(); j++)
		{
			m_face = faces[m_block->connected_faces[j]];
			flow_face = m_face->phases[phase_index]->q_norm * m_face->face_area;
			//std::cout << "flow_face " << flow_face << "\n";
			if (fabs(flow_face) > flow_max)
				flow_max = fabs(flow_face);
			if ((j == 1) || (j == 3) || (j == 5)) // minus faces: -i, -j or -k
				flow_face *= -1.0;
			//std::cout << " face " << m_face->index << ": " << flow_face << " ;   " << "\n";
			flow += flow_face;
		}
		//std::cout << " total flow: " << flow << "\n";
		if ((fabs(flow / flow_max) > 1.0e-3) && (fabs(flow) > 1.0e-10))
		{
			if (max_error < abs(flow / flow_max))
				max_error = abs(flow / flow_max);
			std::cout << " ERROR in budget for block " << i << " :  Sum_flow: " << flow <<
				", max_flow : " << flow_max << ", rel_error: " << flow / flow_max <<
				", max_error: " << max_error << "\n";
		}
		//std::cout << " max_error: " << max_error << "\n";
		//if (max_error > 0.0043)
		//    std::cout << "\n";
	}
	return true;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: MakeNodeVector
   Task: XXX
   Return: true or false
   Programming: XXX
   Modification:
   -------------------------------------------------------------------------*/
//bool CECLIPSEData::MakeNodeVector(CRFProcess *m_pcs, std::string path, int timestep, int phase_index)
bool CECLIPSEData::MakeNodeVector(void)
{
	CFEMesh* m_msh = fem_msh_vector[0]; //SB: ToDo hart gesetzt
	//CFaces *m_face=NULL;
	//WW double weights_xyz[3];
	CPointData_ECL* m_NodeData = NULL;
	m_NodeData = new CPointData_ECL;
	vector <double> temp_q;

	//double coord_v_x[3] = {1.0,0.0,0.0};
	//double coord_v_y[3] = {0.0,1.0,0.0};
	//double coord_v_z[3] = {0.0,0.0,1.0};
	//double * normal_vec_face;
	//double distance;
	//double weight;
	//double sum_weights;
	//double val;

	for (int i = 0; i < 3; i++)
		temp_q.push_back(-1.0E+99);

	/* Go through GeoSys mesh, and creates CPointData - then Nodes in
	   CECLIPSEData::NodeData have the same order as in the Geosys mesh,
	   and can be used directly for interpolation to gauss points*/
	//this->NodeData.clear();
	if (this->NodeData.size() < 1)
	{
		for (unsigned long i = 0; i < m_msh->nod_vector.size(); i++)
		{
			// create new instance of CPointData
			// set m_NodeData (x y z) as the same geosys node order's (x y z)
			m_NodeData = new CPointData_ECL(m_msh->nod_vector[i]->getData());
			//Get the node
			// TF	m_node = m_msh->nod_vector[i];
			// TF	m_NodeData->x = m_node->X();
			// TF	m_NodeData->y = m_node->Y();
			// TF	m_NodeData->z = m_node->Z();
			m_NodeData->phase_pressure.resize(3);
			m_NodeData->phase_saturation.resize(3);
			//WTP m_NodeData->phase_density.resize(3);
			m_NodeData->phase_reservoir_density.resize(3);	// WTP 04/2014
			m_NodeData->phase_surface_density.resize(3);
			m_NodeData->phase_viscosity.resize(3);
			m_NodeData->phase_fvf.resize(3);

			//for (int j = 0; j < int(this->Phases.size()); j++)
			//	m_NodeData->q.push_back(temp_q);

			//Set variable to zero
			m_NodeData->pressure = -1.0E+99;
			//m_NodeData->CO2inLiquid = -1.0E+99;
			//m_NodeData->deltaDIC = -1.0E+99;
			for (long k = 0; k < 3; k++)
			{
				//m_NodeData->Flow[k] = 0.0;
				m_NodeData->phase_pressure[k] = -1.0E+99;
				m_NodeData->phase_saturation[k] = -1.0E+99;
				//WTP m_NodeData->phase_density[k] = -1.0E+99;
				m_NodeData->phase_reservoir_density[k] = -1.0E+99;	//WTP 04/2014
				m_NodeData->phase_surface_density[k] = -1.0E+99;
				m_NodeData->phase_viscosity[k] = -1.0E+99;
				m_NodeData->phase_fvf[k] = -1.0E+99;
				//WW		weights_xyz[k] = 0.0;
			}

			// transfer Data to node
			this->NodeData.push_back(m_NodeData);

			// Test output
			// std::cout << " Node " << i << " (X,Y,Z): (" << m_NodeData->x << ", " << m_NodeData->y << ", "<< m_NodeData->z << ") ";
			// for (long k = 0; k < 3; k++){
			//     std::cout << m_NodeData->phase_density[k] << " " << m_NodeData->phase_pressure[k] << " " << m_NodeData->phase_saturation[k] << " ";
			// }
			// std::cout << "\n";
			// // Test output
			// std::cout << " Node " << i << " (X,Y,Z): (" << this->NodeData[this->NodeData.size()-1]->x << ", " << 
			//this->NodeData[this->NodeData.size()-1]->y << ", "
			//<< this->NodeData[this->NodeData.size()-1]->z << ") ";
			// for (long k = 0; k < 3; k++){
			//     std::cout << this->NodeData[this->NodeData.size()-1]->phase_density[k] << " " << 
			//this->NodeData[this->NodeData.size()-1]->phase_pressure[k] << " " 
			//<< this->NodeData[this->NodeData.size()-1]->phase_saturation[k] << " ";
			// }
			// std::cout << "\n";
		}
		//BW: Assign Connected Blocks to Nodevector
		for (size_t i = 0; i < this->eclgrid.size(); i++)
		{
			for (size_t j = 0; j < this->eclgrid[i]->Nodeindex.size(); j++)
			{
				bool findblock = false;
				if (this->NodeData[this->eclgrid[i]->Nodeindex[j]]->connect_block.size() == 0)
					this->NodeData[this->eclgrid[i]->Nodeindex[j]]->connect_block.push_back(this->eclgrid[i]->index);
				else
				{
					for (size_t k = 0; k < this->NodeData[this->eclgrid[i]->Nodeindex[j]]->connect_block.size(); k++)
						if (static_cast<long>(this->NodeData[this->eclgrid[i]->Nodeindex[j]]->connect_block[k]) == this->eclgrid[i]->index)
							findblock = true;
					if (findblock == false)
						this->NodeData[this->eclgrid[i]->Nodeindex[j]]->connect_block.push_back(this->eclgrid[i]->index);
				}
			}
		}
	}
	else
	{
		for (unsigned long i = 0; i < m_msh->nod_vector.size(); i++)
		{
			// Set variable to zero
			this->NodeData[i]->pressure = -1.0E+99;
			//WTP this->NodeData[i]->CO2inLiquid = -1.0E+99;
			for (long k = 0; k < 3; k++)
			{
				//this->NodeData[i]->Flow[k] = 0.0;
				this->NodeData[i]->phase_pressure[k] = -1.0E+99;
				this->NodeData[i]->phase_saturation[k] = -1.0E+99;
				//WTP this->NodeData[i]->phase_density[k] = -1.0E+99;
				this->NodeData[i]->phase_reservoir_density[k] = -1.0E+99;	// WTP 04/2014
				this->NodeData[i]->phase_surface_density[k] = -1.0E+99;
				this->NodeData[i]->phase_viscosity[k] = -1.0E+99;
				this->NodeData[i]->phase_fvf[k] = -1.0E+99;
			}
		}
	}
	// for (unsigned long i = 0; i < m_msh->nod_vector.size(); i++) {
	//                             //    //Get the node
	//     m_NodeData = this->NodeData[i];
	//     m_node = m_msh->nod_vector[i];
	//     m_NodeData->x = m_node->X();
	//     m_NodeData->y = m_node->Y();
	//     m_NodeData->z = m_node->Z();
	//     //sum_weights = 0;
	//     // sum the distance weighted data from each connected face
	//     for (unsigned int j = 0; j < m_node->connected_faces.size(); j++) {
	//       distance = weight = 0;
	//       m_face = this->faces[m_node->connected_faces[j]];
	//       distance = m_node->distance_to_connected_faces[j];
	//       //Weight of each face depending on distance
	//       weight = (1.0 / distance);
	//       // Sum of weights
	//       //sum_weights += weight;
	//       normal_vec_face = m_face->PlaneEquation->GetNormalVector();
	//       // Go through all three coordinates x, y, z and check, if face is perpendicular to axis
	//       for (int k = 0; k < 3; k++){
	//           if(k == 0) val = fabs(PointProduction(normal_vec_face,coord_v_x));
	//           if(k == 1) val = fabs(PointProduction(normal_vec_face,coord_v_y));
	//           if(k == 2) val = fabs(PointProduction(normal_vec_face,coord_v_z));
	//           if(val > MKleinsteZahl){ // face not perpendicular to ccordinate axis k
	//             //m_NodeData->Flow[k] += m_face->phases[phase_index]->q[k] * weight;
	//             weights_xyz[k] += weight;
	//             //std::cout << " Node " << i << " contributed by face " << m_face->index << ", " << m_face->model_axis << "  , " 
	//             //<< m_face->q_norm << ": ";
	//             //for(int mm=0;mm<3;mm++) std::cout << m_face->q[mm] << ", ";
	//             //std::cout << "\n";
	//           }
	//       }
	//     }
	// 
	//     // normalize weighted sum by sum_of_weights sum_w
	//     for (int k = 0; k < 3; k++) {
	//      // if(weights_xyz[k] > 0.0)
	//         ////m_NodeData->Flow[k] = m_NodeData->Flow[k] / weights_xyz[k];
	//      // else{
	//         ////m_NodeData->Flow[k] = m_NodeData->Flow[k];
	//         //std::cout << " Warning - no faces for direction (I=0,J=1,K=2): " << k << " at node " << i << 
	//         //" with non-zero contributions! " << "\n";
	//      // }
	//     }
	// }
	// // TEst Output
	// std::string tempstring;
	// ostringstream temp;
	// vector <string> vec_string;
	// double v_geosys[3];
	// vec_string.push_back("Node; X; Y; Z; v_x_ECL; v_y_ECL; v_z_ECL; v_ECL; v_x_Geos; v_y_Geos; v_z_Geos; v_Geos");
	// 
	// for(unsigned long i=0; i< this->NodeData.size();i++){
	//     m_NodeData = this->NodeData[i];
	//     //std::cout << " NodeData["<< i << "]:  (X,Y,Z): (" << m_NodeData->x << ", " << m_NodeData->y << ", " << 
	//     //m_NodeData->z << "):        " ;
	//     tempstring="";
	//     temp.str(""); temp.clear(); temp << i; tempstring = temp.str();
	//     temp.str(""); temp.clear(); temp << m_NodeData->x; tempstring += "; " +temp.str();
	//     temp.str(""); temp.clear(); temp << m_NodeData->y; tempstring += "; " + temp.str();
	//     temp.str(""); temp.clear(); temp << m_NodeData->z; tempstring += "; " + temp.str();
	//     for(int k=0;k< 3;k++)  {
	//         //temp.str(""); temp.clear(); temp << m_NodeData->Flow[k]; tempstring += "; " + temp.str();
	//         //std::cout << m_NodeData->Data_separated[k] << " ";
	//     }
	//     //std::cout << "\n";
	//     //calculate v
	///      temp.str(""); temp.clear(); temp << sqrt(pow(m_NodeData->Flow[0],2) + pow(m_NodeData->Flow[1],2) 
	//         //+ pow(m_NodeData->Flow[2], 2)); tempstring += "; " + temp.str();
	// 
	//     //get node velocity of geosys
	//     InterpolateGeosysVelocitiesToNodes(m_pcs, v_geosys, i);
	//     for(int k=0;k< 3;k++)  {
	//         temp.str(""); temp.clear(); temp << v_geosys[k]; tempstring += "; " + temp.str();
	//     }
	//     //calculate v
	//     temp.str(""); temp.clear(); temp << sqrt(pow(v_geosys[0], 2) + pow(v_geosys[1], 2) + pow(v_geosys[2], 2)); 
	//     //tempstring += "; " + temp.str();
	//     //write string
	//     vec_string.push_back(tempstring);
	// }
	// int position = int(path.find_last_of("\\"));
	// path = path.substr(0,position);
	// position = int(path.find_last_of("\\"));
	// path = path.substr(0,position);
	// temp.str(""); temp.clear(); temp << timestep; tempstring = temp.str();
	// std::string aus_file = path + "\\CheckVelocityAtNodes_" + tempstring + ".csv";
	// ofstream aus;
	// aus.open(aus_file.data(),ios::out);
	// for (unsigned long i = 0; i < vec_string.size(); i++) {
	//     aus << vec_string[i] << "\n";
	// }
	// aus.close();
	return true;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: InterpolateDataFromFacesToNodes
   Task: Interpolates phase velocity from faces (Eclipse) to nodes
   Return: true or false
   Programming: 10/2009 BG / SB
   Modification: 04/2014 WTP: fixed interpolation error due to additional nodes in deformation process
   -------------------------------------------------------------------------*/
void CECLIPSEData::InterpolateDataFromFacesToNodes(long ele_nr,
	double* n_vel_x,
	double* n_vel_y,
	double* n_vel_z,
	int phase_index)
{
	CFEMesh* m_msh = fem_msh_vector[0]; //SB: ToDo hart gesetzt
	MeshLib::CElem* m_ele = NULL;
	CFaces* m_face = NULL;
	MeshLib::CNode* m_node = NULL;
	double distance;
	double weight, weights_xyz[3];
	//double sum_weights;
	double val = 0;
	double coord_v_x[3] = { 1.0, 0.0, 0.0 };
	double coord_v_y[3] = { 0.0, 1.0, 0.0 };
	double coord_v_z[3] = { 0.0, 0.0, 1.0 };
	double data_separated[3];
	double* normal_vec_face;
	//WTP Math_Group::vec <long>nod_index(8);
	//WTP long nod_ind;
	bool choose;
	int face_flag = 1; // option: 0...take all faces at one node, 1...take 2 faces: only the faces which belong to the element
	//todo eingabeoption

	//cout << " InterpolateDataFromFacesToNodes " << "\n" << "\n";

	/* Go through GeoSys mesh, and get each element (=Eclipse block). For each block, look at all corner points
	   (=Nodes) of the element, as phase velocities are needed there. For each node, look at all faces connected to the node.
	   From each face, get the flow accross the face, and make a weighted sum in the nodes.
	   Faces are not considered, if the are perpendicular to the flow direction, as they then would always contribute zero
	   Faces can be choosen with two options:
	   - All faces at one node, then also information of phase velocities from outside the element considered is
	   taken into account.
	   - Or only faces which are directly connected to the element considered.
	   The resulting phase velocities are stored component wise in the vectors n_vel_xyz and passed back.
	   */

	// Get element for which to calculate velocities
	m_ele = m_msh->ele_vector[ele_nr]; // get element
	//WTP m_ele->GetNodeIndeces(nod_index); // get node indexes of nodes on element

	//if ((ele_nr == 646) || (ele_nr == 647) || (ele_nr == 648) || (ele_nr == 683) || (ele_nr == 684) 
	//|| (ele_nr == 685)  || (ele_nr == 720)  ||  (ele_nr == 721)  ||  (ele_nr == 722))
	//	cout << "Zentrum" << "\n";

	//WTP for (long i = 0; i < long(nod_index.Size()); i++) // go through list of connected nodes of the element
	// maybe long is the problem. Why not double? WTP: Why should it be a prob?
	for (long i = 0; i < static_cast<long>(m_ele->GetNodesNumber(false)); i++) // go through list of connected nodes of the element
	{ //Get the connected node
		//WTP nod_ind = nod_index[i];
		long nod_ind = m_ele->GetNodeIndex(i); // WTP: get the node index in one move
		m_node = m_msh->nod_vector[nod_ind]; // m_node = Knoten, der gerade behandelt wird--> siehe Koordinaten
		//cout << " Get Node " << nod_ind << " on element " << ele_nr << "\n";
		// Set variable to zero, initialization
		for (long k = 0; k < 3; k++)
		{
			data_separated[k] = 0.0;
			weights_xyz[k] = 0.0;
		}

		// sum the distance weighted data from each connected face
		//WTP for (unsigned int j = 0; j < m_node->connected_faces.size(); j++)
		for (unsigned int j = 0; j < static_cast<unsigned int>(m_node->connected_faces.size()); j++)
		{
			distance = weight = 0.0;
			m_face = this->faces[m_node->connected_faces[j]]; // Get face on node
			// Consider only if face is on element
			choose = false;

			/*if (m_face->connected_elements[0] == ele_nr)*/
			for (size_t i = 0; i < this->eclgrid[this->faces[m_face->index]->connected_blocks[0]]->correspondingelenum; i++)
			if (this->eclgrid[this->faces[m_face->index]->connected_blocks[0]]->correspondingeleindex[i] == ele_nr)
				choose = true;  //check if face is part of the element


			//if (m_face->connected_elements.size() > 1)
			if (this->faces[m_face->index]->connected_blocks.size() > 1)
			for (size_t i = 0; i < this->eclgrid[this->faces[m_face->index]->connected_blocks[1]]->correspondingelenum; i++)
			if (this->eclgrid[this->faces[m_face->index]->connected_blocks[1]]->correspondingeleindex[i] == ele_nr)
				choose = true;  //check if face is part of the element

			if (face_flag == 0)
				choose = true;

			// for radial flow model
			if (this->Radial_I == true)
				// skip all J faces
				//if (m_face->model_axis.find("J") == 0)
			if (this->faces[m_face->index]->model_axis.find("J") == 0)
				choose = false;
			if (this->Radial_J == true)
				// skip all I faces
				//if (m_face->model_axis.find("I") == 0)
			if (this->faces[m_face->index]->model_axis.find("I") == 0)
				choose = false;

			if (!choose)
				continue;
			else
			{
				//Testoutput
				//std::cout << "\n";
				//std::cout << " Node " << nod_ind << " contributed by face " << m_face->index << ", " 
				//<< m_face->model_axis << "  , " << m_face->phases[phase_index]->q_norm << ": ";
				//for(int mm=0;mm<3;mm++)std::cout << m_face->phases[phase_index]->q[mm] << ", ";
				//std::cout << "\n";

				distance = m_node->distance_to_connected_faces[j];
				//Weight of each face depending on distance
				weight = (1.0 / distance);
				// Sum of weights
				//sum_weights += weight;
				normal_vec_face = m_face->PlaneEquation->GetNormalVector();
				// Go through all three coordinates x, y, z and check, if face is perpendicular to axis
				for (int k = 0; k < 3; k++)
				{
					if (k == 0)
						val =
						fabs(PointProduction(normal_vec_face,
						coord_v_x));
					if (k == 1)
						val =
						fabs(PointProduction(normal_vec_face,
						coord_v_y));
					if (k == 2)
						val =
						fabs(PointProduction(normal_vec_face,
						coord_v_z));
					if (val > MKleinsteZahl) // face not perpendicular to ccordinate axis k
					{
						data_separated[k] +=
							m_face->phases[phase_index]->q[k] * weight;
						weights_xyz[k] += weight;
					}
				}
			} // end if face is on element
		} // end loop faces at one node of element considered

		// normalize weighted sum by sum_of_weights sum_w
		for (int k = 0; k < 3; k++)
		{
			if (weights_xyz[k] > 0.0)
				data_separated[k] = data_separated[k] / weights_xyz[k];
			else
			{
				data_separated[k] = data_separated[k];
				if ((k == 1) && (this->Radial_I == true))
				{
				}          // do nothing, Radial model perpendicular x axis
				else if ((k == 0) && (this->Radial_J == true))
				{
				}                // do nothing, Radial model perpendicular y axis
				else if (verbosity > 1)
				{
					std::cout << " WARNING: no faces for direction (I=0,J=1,K=2): " << k << "at node " << nod_ind << "\n";
					std::cout << " from element " << ele_nr << " with non-zero contributions! " << "\n";
				}
			}
		}
		// transfer Data to node
		n_vel_x[i] = data_separated[0];
		n_vel_y[i] = data_separated[1];
		n_vel_z[i] = data_separated[2];

		//Test output
		//std::cout << " Node " << nod_ind << " (X,Y,Z): (" << m_node->X() << ", " << m_node->Y() << ", "<< m_node->Z() << ") " ;
		//cout << " Node " << nod_ind << " ";
		//std::cout << data_separated[0] << " " << data_separated[1] << " " << data_separated[2] << " " << "\n";
	} // end loop connected nodes

	// Test Output
	//for(long i=0; i< nod_index.Size();i++){
	//	std::cout << " Node "<< nod_index[i] << ":  (X,Y,Z): (" << m_msh->nod_vector[nod_index[i]]->X() << ", " 
	//  << m_msh->nod_vector[nod_index[i]]->Y() << ", "<< m_msh->nod_vector[nod_index[i]]->Z() << "):        " ;
	//	std::cout << n_vel_x[i] << " " << n_vel_y[i] << " " << n_vel_z[i] << " " << "\n";
	//}/* */
}

/*-------------------------------------------------------------------------
   GeoSys - Function: InterpolateDataFromBlocksToNodes
   Task: Interpolates data like phase pressure or saturation from blocks (Eclipse) to nodes (GeoSys)
   Return: true or false
   Programming: 10/2009 SB
   Modification: 04/2014 WTP: added support for multi comp
   -------------------------------------------------------------------------*/
void CECLIPSEData::InterpolateDataFromBlocksToNodes(CRFProcess* m_pcs,
	std::string path,
	int phase_index,
	long Timestep)
{
	(void)path; // unused
	clock_t start, finish;
	double time;
	CFEMesh* m_msh = fem_msh_vector[0]; //SB: ToDo hart gesetzt
	CECLIPSEBlock* m_block = NULL;
	MeshLib::CNode* m_node = NULL;
	//CFaces *m_face=NULL;
	//WW double distance;
	double volume;
	double weight;
	double sum_weights, sum_weights_density, sum_weights_viscosity; // , sum_weights_fvf;
	double saturation, phase_pressure, phase_press;//WTP gas_dissolved, gas_dis;
	// double vapor_component_mass_fraction; // WTP
	double res_den, reservoir_density;
	double surface_density;
	double phase_visc = 0.0, phase_viscosity;

	//CPointData_ECL* m_NodeData = NULL;
	long variable_index_phase_pressure = -1, variable_index_saturation = -1;
	//long variable_index_Gas_dissolved = -1;
	//long variable_index_Water_FVF = -1, variable_index_Gas_FVF = -1, variable_index_Oil_FVF = -1;
	//WTP long variable_index_vapor_mass_fraction = -1;    // SB redo wtp
	long variable_index_Water_density = -1, variable_index_Gas_density = -1, variable_index_Oil_density = -1;
	long variable_index_Water_viscosity = -1, variable_index_Gas_viscosity = -1, variable_index_Oil_viscosity = -1;

	int phase_flag = 0;

	std::vector <double> water_conc_frac;
	std::vector <double> gas_conc_frac;
	std::vector <double> gas_mol_frac;
	std::vector <double> oil_conc_frac;

	start = clock();

	if (verbosity > 2)
		std::cout << "        InterpolateDataFromBlocksToNodes()";

	//get saturation index if there are more than 1 phases
	// the indicies are stored within a vector since the number of variables varies depending on the components in the system
	if (this->Phases[phase_index] == "WATER")
	{
		phase_flag = 1;
		if (static_cast<int>(this->Phases.size()) > 1)
			variable_index_saturation = this->GetVariableIndex("SWAT");
		variable_index_phase_pressure = this->GetVariableIndex("PWAT");
		if (this->E100 == true)
		{
			variable_index_Water_density = this->GetVariableIndex("WAT_DEN");
			variable_index_Water_viscosity = this->GetVariableIndex("WAT_VISC");
		}
		else
		{
			variable_index_Water_density = this->GetVariableIndex("DENW");
			variable_index_Water_viscosity = this->GetVariableIndex("VWAT");
		}
	}
	else
	{
		if (this->Phases[phase_index] == "GAS")
		{
			phase_flag = 2;
			if (static_cast<int>(this->Phases.size()) > 1)
				variable_index_saturation = this->GetVariableIndex("SGAS");
			variable_index_phase_pressure = this->GetVariableIndex("PGAS");
			if (this->E100 == true)
			{
				variable_index_Gas_density = this->GetVariableIndex("GAS_DEN");
				variable_index_Gas_viscosity = this->GetVariableIndex("GAS_VISC");
			}
			else
			{
				variable_index_Gas_density = this->GetVariableIndex("DENG");
				variable_index_Gas_viscosity = this->GetVariableIndex("VGAS");
			}
		}
		else
		{
			if (this->Phases[phase_index] == "OIL")
			{
				phase_flag = 3;
				if (static_cast<int>(this->Phases.size()) > 1)
					variable_index_saturation = this->GetVariableIndex("SOIL");
				variable_index_phase_pressure = this->GetVariableIndex("POIL");
				if (this->E100 == true)
				{
					variable_index_Oil_density = this->GetVariableIndex("OIL_DEN");
					variable_index_Oil_viscosity = this->GetVariableIndex("OIL_VISC");
				}
				else
				{
					variable_index_Oil_density = this->GetVariableIndex("DENO");
					variable_index_Oil_viscosity = this->GetVariableIndex("VOIL");
				}
				//variable_index_Oil_FVF = this->GetVariableIndex("BO");
			}
			else
			{
				std::cout << " ERROR: " << Phases[phase_index] << " is not considered yet!" << "\n";
				//system("Pause");
				exit(0);
			}
		}
	}

	// Run over all nodes, collect the data and interpolate it to element values
	for (unsigned long i = 0; i < m_msh->nod_vector.size(); i++)
	{
		//Get the node
		m_node = m_msh->nod_vector[i];
		// set all local variables to defaults
		phase_pressure = 0.0;
		phase_visc = 0.0;
		phase_viscosity = 0.0;
		saturation = 0.0;
		sum_weights = 0.0;
		//WTP gas_dissolved = 0.0;
		reservoir_density = 0.0;
		surface_density = 0.0;
		sum_weights_density = 0.0;
		sum_weights_viscosity = 0.0;
		//sum_weights_fvf = 0.0;
		//vapor_component_mass_fraction = 0.0; // SB redo wtp

		// wtp clear the vectors
		water_conc_frac.clear();
		gas_conc_frac.clear();
		oil_conc_frac.clear();
		gas_conc_frac.clear();

		//WTP 04/2014 initialize the vector needed for the gas composition but only if the values are given
		if (vec_CompConc_Water_elements.size() > 0)
		for (unsigned int k = 0; k<vec_CompConc_Water_elements[0].size(); k++)
			water_conc_frac.push_back(0.0);
		if (vec_CompConc_Gas_elements.size() > 0)
		for (unsigned int k = 0; k<vec_CompConc_Gas_elements[0].size(); k++)
			gas_conc_frac.push_back(0.0);
		if (vec_CompConc_Oil_elements.size() > 0)
		for (unsigned int k = 0; k<vec_CompConc_Oil_elements[0].size(); k++)
			oil_conc_frac.push_back(0.0);
		if (vec_mComp_gas_elements.size() > 0)
		for (unsigned int k = 0; k < vec_mComp_gas_elements[0].size(); k++)
			gas_mol_frac.push_back(0.0);

		// sum the distance weighted data from each connected block
		/*for (unsigned int j = 0; j < m_node->getConnectedElementIDs().size(); j++)*/
		for (size_t j = 0; j < this->NodeData[i]->connect_block.size(); j++)  //BW: go from ECLIPSE BLOCK side
		{
			weight = 0.0;
			/*m_block = this->eclgrid[m_node->getConnectedElementIDs()[j]];*/
			m_block = this->eclgrid[this->NodeData[i]->connect_block[j]];
			//calculate representive volume of the considered node in each connected element for weighting
			//volume = m_block->volume / 8.0; // ToDo volume weighting doesn't work if element is not simple!!!!!!!!! BG
			volume = m_block->volume / this->NodeData[i]->connect_block.size(); //BW: to count as how many blocks are connnected
			weight = (1.0 / volume);

			// Sum of weights
			sum_weights += weight;

			//Phase pressure
			phase_press = this->Data[m_block->index][variable_index_phase_pressure];
			phase_pressure += phase_press * weight;

			// Phase saturation
			if (static_cast<int>(this->Phases.size()) > 1)
			{
				const double sat = this->Data[m_block->index][variable_index_saturation];
				saturation += sat * weight;
			}

			/* WTP 04/2014
			//dissolved gas
			if (variable_index_Gas_dissolved > -1)
			{
			gas_dis = this->Data[m_block->index][variable_index_Gas_dissolved];
			gas_dissolved += gas_dis * weight;
			den = this->Data[m_block->index][variable_index_Gas_density];
			if (den > 0)
			{
			density += den * weight;
			sum_weights_density += weight;
			}
			}
			// KB total molare density , redo wtp
			if (variable_index_vapor_mass_fraction > -1) {
			vap_comp_mass_frac = this->Data[m_block->index][variable_index_vapor_mass_fraction];
			vapor_component_mass_fraction += vap_comp_mass_frac *weight;
			}
			*/
			// first do the stuff for the water phase
			if (phase_flag == 1)
			{
				// first get the dissolved gas comps in water
				if (vec_CompConc_Water_elements.size() > 0)
				{
					for (unsigned int k = 0; k < vec_CompConc_Water_elements[m_block->index].size(); k++)
					{
						double dummy_dbl = water_conc_frac[k];
						double dummy_dbl_add = vec_CompConc_Water_elements[m_block->index][k];
						// weigh it and add it up
						dummy_dbl += dummy_dbl_add * weight;
						// save the total element value
						water_conc_frac[k] = dummy_dbl;
					}
				}
				// Now get the phase densities --> WTP: Are these used yet? 
				if (variable_index_Water_density > -1)
				{
					// get the reservoir density of the phase...
					res_den = this->Data[m_block->index][variable_index_Water_density];
					if (res_den > 0.0)
					{
						reservoir_density += res_den * weight;
						sum_weights_density += weight;
					}
				}
				// Phase viscosities
				if (variable_index_Water_viscosity > -1)
				{
					phase_visc = this->Data[m_block->index][variable_index_Water_viscosity];
					if (phase_visc > 0.0)
					{
						phase_viscosity += phase_visc * weight;
						sum_weights_viscosity += weight;
					}
				}

			}
			else if (phase_flag == 2) // Now get it for the gas phase
			{
				if (vec_CompConc_Gas_elements.size() > 0)
				{
					for (unsigned int k = 0; k < vec_CompConc_Gas_elements[m_block->index].size(); k++)
					{
						double dummy_dbl = gas_conc_frac[k];
						//double dummy_dbl_add = vec_RS_like_values_elements[m_block->index][k];
						double dummy_dbl_add = vec_CompConc_Gas_elements[m_block->index][k];
						// weigh it and add it up
						dummy_dbl += dummy_dbl_add * weight;
						// save the total element value
						gas_conc_frac[k] = dummy_dbl;
					}
				}
				if (vec_mComp_gas_elements.size() > 0)
				{
					for (unsigned int k = 0; k < vec_mComp_gas_elements[m_block->index].size(); k++)
					{
						double dummy_dbl = gas_mol_frac[k];
						//double dummy_dbl_add = vec_RS_like_values_elements[m_block->index][k];
						double dummy_dbl_add = vec_mComp_gas_elements[m_block->index][k];
						// weigh it and add it up
						dummy_dbl += dummy_dbl_add * weight;
						// save the total element value
						gas_mol_frac[k] = dummy_dbl;
					}
				}
				if (variable_index_Gas_density > -1)
				{
					// get the reservoir density
					res_den = this->Data[m_block->index][variable_index_Gas_density];
					if (res_den > 0.0)
					{
						reservoir_density += res_den * weight;
						sum_weights_density += weight;
					}
				}
				if (variable_index_Gas_viscosity > -1)
				{
					phase_visc = this->Data[m_block->index][variable_index_Gas_viscosity];
					if (phase_visc > 0.0)
					{
						phase_viscosity += phase_visc * weight;
						sum_weights_viscosity += weight;
					}
				}
			}
			// Now the oil phase
			else if (phase_flag == 3)
			{ 
				if(vec_CompConc_Oil_elements.size() > 0)
				{
					for (unsigned int k = 0; k < vec_CompConc_Oil_elements[m_block->index].size(); k++)
					{
						double dummy_dbl = oil_conc_frac[k];
						//double dummy_dbl_add = vec_RS_like_values_elements[m_block->index][k];
						double dummy_dbl_add = vec_CompConc_Oil_elements[m_block->index][k];
						// weigh it and add it up
						dummy_dbl += dummy_dbl_add * weight;
						// save the total element value
						oil_conc_frac[k] = dummy_dbl;
					}
				}
				if (variable_index_Oil_density > -1)
				{
					// first the reservoir density
					res_den = this->Data[m_block->index][variable_index_Oil_density];
					if (res_den > 0.0)
					{
						reservoir_density += res_den * weight;
						sum_weights_density += weight;
					}
				}
				if (variable_index_Oil_viscosity > -1)
				{
					phase_visc = this->Data[m_block->index][variable_index_Oil_viscosity];
					if (phase_visc > 0.0)
					{
						phase_viscosity += phase_visc * weight;
						sum_weights_viscosity += weight;
					}
				}
			}

			// Test output
			//if ((m_block->z_barycentre == -2001.8)) {
			//    //std::cout << " Node " << i << "\n";
			//    //std::cout << " :  connected block: " << m_node->connected_elements[j] << " (X,Y,Z): " << m_block->x_barycentre << ", " 
			//<< m_block->y_barycentre << ", " << m_block->z_barycentre << " with saturation " << sat;
			//    //std::cout << "  distance: " << distance <<  " weight: " << weight << " sum of weights: " << sum_weights << "\n";
			//    std::cout <<  "      Node " << i << " Element: " << m_node->connected_elements[j] << " SGAS: " << sat << "\n";
			//}
		} // for j=connected_elements

		// Now weigh the values according to their distribution
		// phase pressure
		phase_pressure = phase_pressure / sum_weights;
		this->NodeData[i]->phase_pressure[phase_index] = phase_pressure;
		this->NodeData[i]->pressure = phase_pressure; // the last phase pressure is the highest pressure

		// in a multiphase case: phase saturation
		if (static_cast<int>(this->Phases.size()) > 1)
		{
			saturation = saturation / sum_weights;
			if ((saturation >= 0.0) && (saturation <= 1.0))
				this->NodeData[i]->phase_saturation[phase_index] = saturation;
			else
			{
				if (verbosity > 1)
					std::cout << " WARNING: The saturation is not between 0 and 1! Values are corrected!" << "\n";
				if (saturation < 0.0)
					this->NodeData[i]->phase_saturation[phase_index] = 0.0;
				if (saturation > 1.0)
					this->NodeData[i]->phase_saturation[phase_index] = 1.0;
			}
			//// Testoutput BW
			//if (std::fabs(saturation - 0.0) != 0 && this->Phases[phase_index] == "GAS")
			//if (i==22958)
			//{
			//	cout << "Node " << i << " Connectedblocknumber " << this->NodeData[i]->connect_block.size() << " Sum_Weights: " << sum_weights << " saturaion " << saturation << '\n';
			//	for (int jjj = 0; jjj < this->NodeData[i]->connect_block.size(); jjj++)
			//	{
			//		cout << "Connected Blocks: " << this->NodeData[i]->connect_block[jjj]<<" ";
			//		cout << "Saturation: " << this->Data[this->eclgrid[this->NodeData[i]->connect_block[jjj]]->index][variable_index_saturation] << " ";
			//		cout << "Volume: " << this->eclgrid[this->NodeData[i]->connect_block[jjj]]->volume << '\n';
			//	}
			//}

		}

		// phase density
		if (sum_weights_density > 0)
		{
			reservoir_density = reservoir_density / sum_weights_density;
			surface_density = surface_density / sum_weights_density;
		}
		this->NodeData[i]->phase_reservoir_density[phase_index] = reservoir_density;
		this->NodeData[i]->phase_surface_density[phase_index] = surface_density;

		// phase viscosity
		if (sum_weights_viscosity > 0)
		{
			phase_viscosity = phase_viscosity / sum_weights_viscosity;
		}
		this->NodeData[i]->phase_viscosity[phase_index] = phase_viscosity;


		/* WTP 04/2014
		//dissolved gas in oil (black oil mode)
		if (variable_index_Gas_dissolved > -1)
		{
		gas_dissolved = gas_dissolved / sum_weights;
		this->NodeData[i]->CO2inLiquid = gas_dissolved;
		}

		// KB total molare density, redo wtp
		if (variable_index_vapor_mass_fraction > -1){
		vapor_component_mass_fraction = vapor_component_mass_fraction / sum_weights;
		this->NodeData[i]->VaporComponentMassFraction = vapor_component_mass_fraction;
		}
		*/

		// WTP Store the data (write it back to OGS)
		if (phase_flag == 1)
		{
			if (vec_CompConc_Water_elements.size() > 0)
			{
				//for(size_t k=0; k<vec_gas_comps_inLiquid_indicies.size(); k++) 
				for (size_t k = 0; k < vec_CompConc_Water_elements[m_block->index].size(); k++)
				{
					double dummy_dbl = water_conc_frac[k] / sum_weights;

					// store the data
					if (Timestep == 1 && m_pcs->iter_outer_cpl == 0)
					{
						this->NodeData[i]->CompConcInWater_Nodes.push_back(dummy_dbl);
						this->NodeData[i]->delta_CompConcInWater_Nodes.push_back(0.0);
					}
					else
					{
						this->NodeData[i]->CompConcInWater_Nodes[k] = dummy_dbl;
					}

				}
			}
		}
		else if (phase_flag == 2)
		{
			if (this->vec_CompConc_Gas_elements.size() > 0)
			{
				for (size_t k = 0; k < vec_CompConc_Gas_elements[m_block->index].size(); k++)
				{
					double dummy_dbl = gas_conc_frac[k] / sum_weights;
					// store the data
					if (Timestep == 1 && m_pcs->iter_outer_cpl == 0)
					{
						this->NodeData[i]->CompConcInGas_Nodes.push_back(dummy_dbl);
						this->NodeData[i]->delta_CompConcInGas_Nodes.push_back(0.0);
					}
					else
					{
						this->NodeData[i]->CompConcInGas_Nodes[k] = dummy_dbl;
					}
				}
			}
			if (this->vec_mComp_gas_elements.size() > 0)
			{
				for (size_t k = 0; k < vec_mComp_gas_elements[m_block->index].size(); k++)
				{
					double dummy_dbl = gas_mol_frac[k] / sum_weights;
					if (Timestep == 1 && m_pcs->iter_outer_cpl == 0)
						this->NodeData[i]->mCompInGas_Nodes.push_back(dummy_dbl);
					else
						this->NodeData[i]->mCompInGas_Nodes[k] = dummy_dbl;
				}
			}
		}
		else if (phase_flag == 3)
		{
			if (this->vec_CompConc_Oil_elements.size() > 0)
			{
				for (size_t k = 0; k < vec_CompConc_Oil_elements[m_block->index].size(); k++)
				{
					double dummy_dbl = oil_conc_frac[k] / sum_weights;
					if (Timestep == 1 && m_pcs->iter_outer_cpl == 0)
					{
						this->NodeData[i]->CompConcInOil_Nodes.push_back(dummy_dbl);
						this->NodeData[i]->delta_CompConcInOil_Nodes.push_back(0.0);
					}
					else
					{
						this->NodeData[i]->CompConcInOil_Nodes[k] = dummy_dbl;
					}
				}
			}
		}
	}

	//Test output
	//if ((m_block->z_barycentre == -2001.8)) {
	//	//cout << " Node " << i << " (X,Y,Z): (" << m_node->X() << ", " << m_node->Y() << ", "<< m_node->X() << ") ";
	//	//cout << " Pressure at Node " << this->NodeData[i]->pressure << "\n";
	//	//cout << " Node " << i << " (X,Y,Z): (" << m_node->X() << ", " << m_node->Y() << ", "<< m_node->X() << ") ";
	//	//cout << " Saturation at Node " << this->NodeData[i]->phase_saturation[phase_index] << "\n";
	//	cout << " Node " << i << " SGAS: " <<  this->NodeData[i]->phase_saturation[phase_index] << "\n";
	//}

	//// Test Output Nodes
	////if ((m_pcs->Tim->step_current==0) || (m_pcs->Tim->step_current==1) || (m_pcs->Tim->step_current==2) || (m_pcs->Tim->step_current==5) ||
	//(m_pcs->Tim->step_current==10) || (m_pcs->Tim->step_current==20) || (m_pcs->Tim->step_current==30) || (m_pcs->Tim->step_current==40) || 
	//(m_pcs->Tim->step_current==50) || (m_pcs->Tim->step_current==100) || (m_pcs->Tim->step_current==999) || (m_pcs->Tim->step_current==1000)) 
	//{
	//vector <string> vec_string;
	//std::string tempstring;
	//ostringstream temp;
	//double Val;
	//if (this->Phases.size() == 1)
	//{
	//    int nidx1 = m_pcs->GetNodeValueIndex("PRESSURE1") + 1; //+1... new time level
	//    vec_string.push_back("Node; X; Y; Z; P-Geosys; P-Eclipse");
	//    for(unsigned long i = 0; i < this->NodeData.size(); i++)
	//    {
	//        m_NodeData = this->NodeData[i];
	//        //Get the Pressure from the Geosys run
	//        Val = m_pcs->GetNodeValue(i, nidx1);
	//        //std::cout << " NodeData[" << i << "]:  (X,Y,Z): (" << m_NodeData->x << ", " << m_NodeData->y << ", " << m_NodeData->z << "):        " ;
	//        //std::cout << ", Geosys-P: " << Val << ", Eclipse: " << this->NodeData[i]->pressure << "\n";;
	//        temp.str("");
	//        temp.clear();
	//        temp << i;
	//        tempstring = temp.str();
	//        temp.str("");
	//        temp.clear();
	//        temp << m_NodeData->x;
	//        tempstring += "; " + temp.str();
	//        temp.str("");
	//        temp.clear();
	//        temp << m_NodeData->y;
	//        tempstring += "; " + temp.str();
	//        temp.str("");
	//        temp.clear();
	//        temp << m_NodeData->z;
	//        tempstring += "; " + temp.str();
	//        temp.str("");
	//        temp.clear();
	//        temp.precision(12);
	//        temp << Val;
	//        tempstring += "; " + temp.str();
	//        temp.str("");
	//        temp.clear();
	//        temp.precision(12);
	//        temp << this->NodeData[i]->pressure;
	//        tempstring += "; " + temp.str();
	//        vec_string.push_back(tempstring);
	//    }
	//}
	//if (this->Phases.size() > 1)
	//    if (phase_index > 0)
	//    {
	//        int nidx1 = m_pcs->GetNodeValueIndex("PRESSURE1") + 1; //+1... new time level
	//        int nidx2 = m_pcs->GetNodeValueIndex("PRESSURE2") + 1;
	//        int nidx3 = m_pcs->GetNodeValueIndex("SATURATION1") + 1;
	//        //vec_string.push_back("Node; X; Y; Z; P-Geosys_phase1; P-Geosys_phase2; Saturation_Geosys_phase1; P-Eclipse_phase1; 
	//        //P-Eclipse_phase2; Saturation_Eclipse_phase1");
	//        vec_string.push_back(
	//                "Node; X; Y; Z; P-Geosys_phase1; P-Geosys_phase2; Saturation_Geosys_phase1; P-Eclipse_phase1; P-Eclipse_phase2; 
	//                //Saturation_Eclipse_phase1; 
	//q_x1; q_y1; q_z1; q_x2; q_y2; q_z2");
	//        for(unsigned long i = 0; i < this->NodeData.size(); i++)
	//        {
	//            m_NodeData = this->NodeData[i];
	//            //std::cout << " NodeData[" << i << "]:  (X,Y,Z): (" << m_NodeData->x << ", " << m_NodeData->y << ", "
	//            // << m_NodeData->z << "):        " ;
	//            //std::cout << ", Geosys-P: " << Val << ", Eclipse: " << this->NodeData[i]->pressure << "\n";;
	//            temp.str("");
	//            temp.clear();
	//            temp << i;
	//            tempstring = temp.str();
	//            temp.str("");
	//            temp.clear();
	//            temp << m_NodeData->x;
	//            tempstring += "; " + temp.str();
	//            temp.str("");
	//            temp.clear();
	//            temp << m_NodeData->y;
	//            tempstring += "; " + temp.str();
	//            temp.str("");
	//            temp.clear();
	//            temp << m_NodeData->z;
	//            tempstring += "; " + temp.str();

	//            //Get the Pressure1 from the Geosys run
	//            Val = m_pcs->GetNodeValue(i, nidx1);
	//            temp.str("");
	//            temp.clear();
	//            temp.precision(12);
	//            temp << Val;
	//            tempstring += "; " + temp.str();
	//            //Get the Pressure2 from the Geosys run
	//            Val = m_pcs->GetNodeValue(i, nidx2);
	//            temp.str("");
	//            temp.clear();
	//            temp.precision(12);
	//            temp << Val;
	//            tempstring += "; " + temp.str();
	//            //Get the Saturation1 from the Geosys run
	//            Val = m_pcs->GetNodeValue(i, nidx3);
	//            temp.str("");
	//            temp.clear();
	//            temp.precision(12);
	//            temp << Val;
	//            tempstring += "; " + temp.str();

	//            //Eclipse data
	//            temp.str("");
	//            temp.clear();
	//            temp.precision(12);
	//            temp << this->NodeData[i]->phase_pressure[0];
	//            tempstring += "; " + temp.str();
	//            temp.str("");
	//            temp.clear();
	//            temp.precision(12);
	//            temp << this->NodeData[i]->phase_pressure[1];
	//            tempstring += "; " + temp.str();
	//            temp.str("");
	//            temp.clear();
	//            temp.precision(12);
	//            temp << this->NodeData[i]->phase_saturation[0];
	//            tempstring += "; " + temp.str();

	//            ////Velocities
	//            //temp.str(""); temp.clear(); temp.precision(12); temp << this->NodeData[i]->q[0][0]; tempstring += "; " + temp.str();
	//            //temp.str(""); temp.clear(); temp.precision(12); temp << this->NodeData[i]->q[0][1]; tempstring += "; " + temp.str();
	//            //temp.str(""); temp.clear(); temp.precision(12); temp << this->NodeData[i]->q[0][2]; tempstring += "; " + temp.str();
	//            //temp.str(""); temp.clear(); temp.precision(12); temp << this->NodeData[i]->q[1][0]; tempstring += "; " + temp.str();
	//            //temp.str(""); temp.clear(); temp.precision(12); temp << this->NodeData[i]->q[1][1]; tempstring += "; " + temp.str();
	//            //temp.str(""); temp.clear(); temp.precision(12); temp << this->NodeData[i]->q[1][2]; tempstring += "; " + temp.str();

	//            vec_string.push_back(tempstring);
	//        }
	//    }
	//    
	//CWriteTextfiles *TextFile;
	//TextFile = new CWriteTextfiles;
	//int position = int(path.find_last_of("\\"));
	//std::string path_new;
	//path_new = path.substr(0,position);
	//position = int(path_new.find_last_of("\\"));
	//path_new = path_new.substr(0,position);
	//temp.str(""); temp.clear(); temp << m_pcs->Tim->step_current; tempstring = temp.str();
	//TextFile->Write_Text(path_new + "\\CheckNodeValues_" + tempstring + ".csv", vec_string);
	//} // End Test Output


	finish = clock();
	time = (static_cast<double>(finish)-static_cast<double>(start)) / CLOCKS_PER_SEC;
	if (verbosity > 2)
	{
		std::cout << "    Time: " << time << " seconds." << "\n";
		std::cout << flush;
	}
}

/*-------------------------------------------------------------------------
   GeoSys - Function: InterpolateGeosysVelocitiesToNodes
   Task: Returns the Nodevelocity of Geosys; v at the node as a inverse distance
   weighted mean of the connecting elements velocities
   Return: Velocity
   Programming: 11/2009 BG based on CB
   Modification:
   09/2011	TF changed access to coordinates of mesh node,
   - substituted access to mesh_element from pointer to direct access into the vector
   - made the mesh node a const pointer
   - made the pointer to the mesh const, made the mesh itself const
   - substituted pow(x,2) by x*x
   -------------------------------------------------------------------------*/
void CECLIPSEData::InterpolateGeosysVelocitiesToNodes(CRFProcess* m_pcs, double* vel_nod, long node_number)
{
	CFEMesh const* const mesh = fem_msh_vector[0]; //SB: ToDo hart gesetzt
	size_t elem;
	double vel_ele[3];
	double distance, weight, sum_w(0.0);
	double PoreVel(0.0); //WW , theta;

	//double *vel_nod;
	//WW theta = m_pcs->m_num->ls_theta;

	// initialize data structures
	for (size_t i = 0; i < 3; i++)
		vel_nod[i] = vel_ele[i] = 0;

	// Get node coordinates
	MeshLib::CNode const* node(mesh->nod_vector[node_number]);
	double const* const coord(node->getData()); //Coordinates(coord);
	// get the indices of velocity of flow process
	long idxVx = m_pcs->GetElementValueIndex("VELOCITY1_X") + 1;
	long idxVy = m_pcs->GetElementValueIndex("VELOCITY1_Y") + 1;
	long idxVz = m_pcs->GetElementValueIndex("VELOCITY1_Z") + 1;

	for (size_t el = 0; el < node->getConnectedElementIDs().size(); el++)
	{
		distance = weight = 0; // initialize for each connected element
		elem = node->getConnectedElementIDs()[el];

		// get the velocity components of element elem
		// if idxVx = 0 it implies that the index for this parameter doesn't exist
		if (idxVx != 0)
			vel_ele[0] = m_pcs->GetElementValue(elem, idxVx);
		if (idxVy != 0)
			vel_ele[1] = m_pcs->GetElementValue(elem, idxVy);
		if (idxVz != 0)
			vel_ele[2] = m_pcs->GetElementValue(elem, idxVz);
		// calculate distance node <-> element center of gravity
		double const* grav_c(mesh->ele_vector[elem]->GetGravityCenter());
		for (size_t i = 0; i < 3; i++)
			distance += (coord[i] - grav_c[i]) * (coord[i] - grav_c[i]);  // TF pow((coord[i]-grav_c[i]),2);
		// linear inverse distance weight = 1/(distance)
		distance = sqrt(distance); // for quadratic interpolation uncomment this line
		weight = (1 / distance);
		sum_w += weight;
		for (size_t i = 0; i < 3; i++)
			// the 3 velocity components
			vel_nod[i] += vel_ele[i] * weight;
	}
	// normalize weighted sum by sum_of_weights sum_w
	const double sum_w_inv(1.0 / sum_w);
	for (size_t i = 0; i < 3; i++)
		vel_nod[i] *= sum_w_inv;
	// absolute value of velocity vector
	for (size_t i = 0; i < 3; i++)
		PoreVel += vel_nod[i] * vel_nod[i];  // TF pow(vel_nod[i],2);
	PoreVel = sqrt(PoreVel);
}

/*-------------------------------------------------------------------------
   GeoSys - Function: WriteDataToGeoSys
   Task: Uses the node velocities to calculate the velocity at gauss points
   Return: true or false
   Programming: 09/2009 SB
   Modification: 04/2014 WTP: added support for multi comps & 3 phases
   -------------------------------------------------------------------------*/
void CECLIPSEData::WriteDataToGeoSys(CRFProcess* m_pcs, const std::string path)
{
	//Datenuebergabe: pressure, saturation (ev. kf, n)
	long index_pressure1, index_pressure2, index_pressure3;
	long index_saturation1, index_saturation2, index_saturation3;
	long index_water_density, index_gas_density, index_oil_density;
	long index_water_viscosity, index_gas_viscosity, index_oil_viscosity;
	//long index_water_fvf, index_gas_fvf;
	int phase1, phase2, phase3 = 0; // , transportphase = 0; // WTP
	double value = 0.0; // , value_gas = 0.0; // WTP
	//double saturation_w, 
	//double porosity = 0.1; // SB redo wtp
	//double Molweight_CO2;
	//double epsilon = 1E-15;

	// Sb redo wtp
	// CB_merge_0513
	REACTINT *m_rei = NULL;
	if (REACTINT_vec.size() > 0)
		m_rei = REACTINT_vec[0];

	const clock_t start = clock();
	if (verbosity > 2)
		std::cout << "        WriteDataToGeoSys()" << "\n";

	//WTP_CB: get no of nodes in mesh
	long nnodes = fem_msh_vector[0]->GetNodesNumber(false);

	// Pressure 1 is the wetting phase, Pressure 2 the non wetting phase
	//Pressure
	//WTP if (int(this->Phases.size()) == 1)
	if (static_cast<int>(this->Phases.size()) == 1)
	{
		index_pressure1 = m_pcs->GetNodeValueIndex("PRESSURE1") + 1; //+1... new time level
		index_water_density = m_pcs->GetNodeValueIndex("DENSITY1");  //DENSITY1 is not a variable for some flow process, e.g. LIQUID_FLOW
		//WTP_CB: new way to get the total number of nodes in the mesh since nod_val_vector is not suitable anymore
		//for(unsigned long i = 0; i < m_pcs->nod_val_vector.size(); i++)
		for (long i = 0; i < nnodes; i++)
		{
			value = this->NodeData[i]->pressure;
			// WTP: this should be changed to a more versatile approach
			//fstream datei("K:\\Geomechanics\\Benchmark_Problem1\\OGS_Benchmarks\\Therzagi\\OGS-ECL\\klein\\two-way\\pressure_ECLIPSE.txt", ios::out | ios::app);
			//datei << value << endl;
			m_pcs->SetNodeValue(i, index_pressure1, value);
			//WTP value = this->NodeData[i]->phase_density[0];
            value = this->NodeData[i]->phase_reservoir_density[0];
			m_pcs->SetNodeValue(i, index_water_density, value);
			//datei.close();
		}
	}
	//WTP else if (int(this->Phases.size()) > 1)
	else if (static_cast<int>(this->Phases.size()) > 1)
	{
		if ((this->Phases[0] == "WATER") && (this->Phases[1] == "GAS") &&
			(this->E100 == true))
		{
			std::cout <<
				" ERROR: GAS-WATER systems can not be considered with E100 and GeoSys" <<
				"\n";
            std::cout << flush;
			//system("Pause");
			exit(0);
		}
		//WTP switch(int(this->Phases.size()))
		switch (static_cast<int>(this->Phases.size()))
		{
		case 2:
			phase1 = 0;
			phase2 = 1;
			break;
		case 3:
			// Assumption that there are 3 phases but water is only used for the boundaries -> the oil and gas phase are the relevant 
			//one for the exchange with OGS
			if (this->E100 == true)
			{
				phase1 = 1;
				phase2 = 2;
			}
			else
			{
				phase1 = 0;
				phase2 = 1;
				phase3 = 2;
			}
			break;
		default:
			std::cout << " ERROR: There are not more than 3 phases possible!" << "\n";
            std::cout << flush;
			//system("Pause");
			exit(0);
			break;
		}
		index_pressure1 = m_pcs->GetNodeValueIndex("PRESSURE1") + 1; //+1... new time level
		index_pressure2 = m_pcs->GetNodeValueIndex("PRESSURE2") + 1; //+1... new time level
		index_pressure3 = m_pcs->GetNodeValueIndex("PRESSURE3") + 1; //+1... new time level
		index_saturation1 = m_pcs->GetNodeValueIndex("SATURATION1") + 1; //+1... new time level
		index_saturation2 = m_pcs->GetNodeValueIndex("SATURATION2") + 1; //+1... new time level
		index_saturation3 = m_pcs->GetNodeValueIndex("SATURATION3") + 1; //+1... new time level
		index_water_density = m_pcs->GetNodeValueIndex("DENSITY1");// + 1; // WTP: Also new time level??
		index_gas_density = m_pcs->GetNodeValueIndex("DENSITY2");// + 1;
		index_oil_density = m_pcs->GetNodeValueIndex("DENSITY3");// + 1;
		index_water_viscosity = m_pcs->GetNodeValueIndex("VISCOSITY1");
		index_gas_viscosity = m_pcs->GetNodeValueIndex("VISCOSITY2");
		index_oil_viscosity = m_pcs->GetNodeValueIndex("VISCOSITY3");

		// WTP: getting the index for the temperature coupling (not yet completed),
		//index_temperature1 = m_pcs->GetNodeValueIndex("TEMPERATURE1") + 1;

		// If 3 Phases are applicable more variables have to be set. 
		if (static_cast<int>(Phases.size()) < 3 || this->E100 == true)
		{
			//WTP_CB: new way to get the total number of nodes in the mesh since nod_val_vector is not suitable anymore
			//for(unsigned long i = 0; i < m_pcs->nod_val_vector.size(); i++)
			for (long i = 0; i < nnodes; i++)
			{
				//// WTP: Saving the temperature data
				//value = this->NodeData[i]->temperature;
				//m_pcs->SetNodeValue(i, index_temperature1, value);

				//calculating capillary pressure (difference of the two phase pressures)
				//WTP: why take the data from pcow and/or pcog ??
				value = this->NodeData[i]->phase_pressure[phase2] -
					this->NodeData[i]->phase_pressure[phase1];
				m_pcs->SetNodeValue(i, index_pressure1, value);

				value = this->NodeData[i]->phase_pressure[phase2];
				m_pcs->SetNodeValue(i, index_pressure2, value);

				value = this->NodeData[i]->phase_reservoir_density[phase1];
				if (index_water_density > 0)
					m_pcs->SetNodeValue(i, index_water_density, value);

				value = this->NodeData[i]->phase_reservoir_density[phase2];
				if (index_gas_density > 0)
					m_pcs->SetNodeValue(i, index_gas_density, value);

				value = this->NodeData[i]->phase_viscosity[phase1];
				if (index_water_viscosity > 0)
					m_pcs->SetNodeValue(i, index_water_viscosity, value);

				value = this->NodeData[i]->phase_viscosity[phase2];
				if (index_gas_viscosity > 0)
					m_pcs->SetNodeValue(i, index_gas_viscosity, value);

				if (this->Phases.size() < 3)
				{
					//Saturation 1 is the wetting phase
					value = this->NodeData[i]->phase_saturation[phase1];
					//cout << " Node: " << i << " saturation1: " << value << "\n";
					m_pcs->SetNodeValue(i, index_saturation1, value);
					value = this->NodeData[i]->phase_saturation[phase2];
					//if (i == 22958)
					//	cout << " Node: " << i << " saturation2: " << value << "\n";
					m_pcs->SetNodeValue(i, index_saturation2, value);
				}
				else
				{
					// the water phase is only used for boundary conditions, the wetting phase is phase 2 and the non-wetting phase is phase 3
					value = 1 - this->NodeData[i]->phase_saturation[phase2];
					//cout << " Node: " << i << " saturation1: " << value << "\n";
					m_pcs->SetNodeValue(i, index_saturation1, value);
					m_pcs->SetNodeValue(i, index_saturation2, this->NodeData[i]->phase_saturation[phase2]);
				}
				////Saturation 1 is the wetting phase
				//value = this->NodeData[i]->phase_saturation[phase1];
				////std::cout << " Node: " << i << " saturation1: " << value << "\n";
				//m_pcs->SetNodeValue(i, index_saturation1, value);
			}
		}
		else if (this->E100 != true) // WTP: ESome special case for 3Phases with E100 was already implemented...
		{
			for (long i = 0; i < nnodes; i++)
			{
				//calculating capillary pressure (difference of the two phase pressures)
				//WTP: why to take the data from pcow and/or pcog ??
				value = this->NodeData[i]->phase_pressure[phase2] -
					this->NodeData[i]->phase_pressure[phase1];
				m_pcs->SetNodeValue(i, index_pressure1, value);

				value = this->NodeData[i]->phase_pressure[phase2];
				m_pcs->SetNodeValue(i, index_pressure2, value);

				// we can only do this if we have 3 phases
				value = this->NodeData[i]->phase_pressure[phase3];
				m_pcs->SetNodeValue(i, index_pressure3, value);

				value = this->NodeData[i]->phase_reservoir_density[phase1];
				if (index_water_density > 0)
					m_pcs->SetNodeValue(i, index_water_density, value);

				value = this->NodeData[i]->phase_reservoir_density[phase2];
				if (index_gas_density > 0)
					m_pcs->SetNodeValue(i, index_gas_density, value);

				value = this->NodeData[i]->phase_reservoir_density[phase3];
				if (index_oil_density > 0)
					m_pcs->SetNodeValue(i, index_oil_density, value);

				value = this->NodeData[i]->phase_viscosity[phase1];
				if (index_water_viscosity > 0)
					m_pcs->SetNodeValue(i, index_water_viscosity, value);

				value = this->NodeData[i]->phase_viscosity[phase2];
				if (index_gas_viscosity > 0)
					m_pcs->SetNodeValue(i, index_gas_viscosity, value);

				value = this->NodeData[i]->phase_viscosity[phase3];
				if (index_oil_viscosity > 0)
					m_pcs->SetNodeValue(i, index_oil_viscosity, value);

				// the water phase is only used for boundary conditions, the wetting phase is phase 2 and the non-wetting phase is phase 3
				value = 1 - this->NodeData[i]->phase_saturation[phase2] - this->NodeData[i]->phase_saturation[phase3];
				//std::cout << " Node: " << i << " saturation1: " << value << "\n";
				m_pcs->SetNodeValue(i, index_saturation1, value);
				//std::cout << " Node: " << i << " saturation1: " << value << "\n";
				value = this->NodeData[i]->phase_saturation[phase2];
				m_pcs->SetNodeValue(i, index_saturation2, value);

				value = this->NodeData[i]->phase_saturation[phase3];
				m_pcs->SetNodeValue(i, index_saturation3, value);
			}
		}
		//WTP 04/2014: Multi comps
		// Alternate way: First get the indicies no matter what version (E100 or E300) is running.
		// Then assign the correct values to the nodes. This way only one single vector structure is needed.

		// DISSOLVED COMPONENTS IN WATER
		// Get the indicies for the dissolved Gas in Water (or rather its components in E300 speak)
		if (this->vec_OGS_process_index_comps_water.size() == 0)
		{
			for (unsigned int l = 0; l < this->vec_components_ECL_OGS_pcs_names.size(); l++)
			{
				CRFProcess *n_pcs = NULL;            // if the component is transported get a process
				for (int k = 0; k < static_cast<int>(pcs_vector.size()); k++)  // get the pcs names
				{
					n_pcs = pcs_vector[k];
					if (n_pcs->nod_val_name_vector[0] == this->vec_components_ECL_OGS_pcs_names[l][1])
						// if the process is found, store the index
						this->vec_OGS_process_index_comps_water.push_back(std::make_pair(vec_components_ECL_OGS_pcs_names[l][1], k));
				}
			}
		};
		// Now run over all the nodes and collect the concentrations
		// Try to use iterator over components instead array + index.
		for (unsigned int j = 0; j < vec_OGS_process_index_comps_water.size(); j++)
		{
			// store the current process index (could be changed)
			const int ProcessIndexCompWater = vec_OGS_process_index_comps_water[j].second;

			// get index of species concentration in nodevaluevector of this process
			const int indexConcentration_water =
				pcs_vector[ProcessIndexCompWater]->GetNodeValueIndex(
				pcs_vector[ProcessIndexCompWater]->pcs_primary_function_name[0]) + 1;    // +1: new timelevel

			// Now run over all nodes and collect the values        
			for (long i = 0; i < nnodes; i++)
			{
				// Store the calculated value for the dissolved gas concentration
				pcs_vector[ProcessIndexCompWater]->SetNodeValue(i, indexConcentration_water, this->NodeData[i]->CompConcInWater_Nodes[j]);

				// check for unphysical value
				if (this->NodeData[i]->CompConcInWater_Nodes[j] < 0. && verbosity > 1)
				{
					std::cout << " WARNING at Node: " << i << " Pressure: " << m_pcs->GetNodeValue(i, index_pressure2) << " Gas in water: " <<
						this->NodeData[i]->CompConcInWater_Nodes[j] << "\n";
					std::cout <<  "        Error in calculation for dissolved gas: " << this->NodeData[i]->CompConcInWater_Nodes[j] << "\n";
				}
			}
		} // end vec_OGS_process_index_comps_water

		// COMPONENTS IN GAS
		// Get the indicies for the components in Gas
		if (this->vec_OGS_process_index_comps_gas.size() == 0)
		{
			for (unsigned int l = 0; l < this->vec_components_ECL_OGS_pcs_names.size(); l++)
			{
				CRFProcess *n_pcs = NULL;            // if the component is transported get a process
				for (int k = 0; k < static_cast<int>(pcs_vector.size()); k++)  // get the pcs names
				{
					n_pcs = pcs_vector[k];
					if (n_pcs->nod_val_name_vector[0] == this->vec_components_ECL_OGS_pcs_names[l][2])
						// if the process is found, store the index
						this->vec_OGS_process_index_comps_gas.push_back(std::make_pair(vec_components_ECL_OGS_pcs_names[l][2], k));
				}
			}
		};

		// Sort the mole fraction data for heat updates

		// Now run over all the nodes and collect the concentrations
		// Try to use iterator over components instead array + index.
		for (unsigned int j = 0; j < vec_OGS_process_index_comps_gas.size(); j++)
		{
			// store the current process index (could be changed)
			const int ProcessIndexCompGas = vec_OGS_process_index_comps_gas[j].second;

			// get index of species concentration in nodevaluevector of this process
			const int indexConcentration_gas =
				pcs_vector[ProcessIndexCompGas]->GetNodeValueIndex(
				pcs_vector[ProcessIndexCompGas]->pcs_primary_function_name[0]) + 1;    // +1: new timelevel

			// Now run over all nodes and collect the values        
			for (long i = 0; i < nnodes; i++)
			{
				// Store the calculated value for the dissolved gas concentration
				pcs_vector[ProcessIndexCompGas]->SetNodeValue(i, indexConcentration_gas, this->NodeData[i]->CompConcInGas_Nodes[j]);

				// check for unphysical value
				if (this->NodeData[i]->CompConcInGas_Nodes[j] < 0. && verbosity > 1)
				{
					std::cout << " WARNING at Node: " << i << " Pressure: " << m_pcs->GetNodeValue(i, index_pressure2) << " Gas in gas: " <<
						this->NodeData[i]->CompConcInGas_Nodes[j] << "\n";
				   std::cout << "          Error in calculation for gas comps : " << this->NodeData[i]->CompConcInGas_Nodes[j] << "\n";
				}
			}
		} // end vec_OGS_process_index_comps_gas

		// COMPONENTS IN OIL
		// Get the indicies for the components in oil
		if (this->vec_OGS_process_index_comps_oil.size() == 0)
		{
			for (unsigned int l = 0; l < this->vec_components_ECL_OGS_pcs_names.size(); l++)
			{
				CRFProcess *n_pcs = NULL;            // if the component is transported get a process
				for (int k = 0; k < static_cast<int>(pcs_vector.size()); k++)  // get the pcs names
				{
					n_pcs = pcs_vector[k];
					if (n_pcs->nod_val_name_vector[0] == this->vec_components_ECL_OGS_pcs_names[l][3])
						// if the process is found, store the index
						this->vec_OGS_process_index_comps_oil.push_back(std::make_pair(vec_components_ECL_OGS_pcs_names[l][3], k));
				}
			}
		};
		// Now run over all the nodes and collect the concentrations
		// Try to use iterator over components instead array + index.
		for (unsigned int j = 0; j < vec_OGS_process_index_comps_oil.size(); j++)
		{
			// store the current process index (could be changed)
			const int ProcessIndexCompOil = vec_OGS_process_index_comps_oil[j].second;

			// get index of species concentration in nodevaluevector of this process
			const int indexConcentration_oil =
				pcs_vector[ProcessIndexCompOil]->GetNodeValueIndex(
				pcs_vector[ProcessIndexCompOil]->pcs_primary_function_name[0]) + 1;    // +1: new timelevel

			// Now run over all nodes and collect the values        
			for (long i = 0; i < nnodes; i++)
			{
				// Store the calculated value for the dissolved gas concentration
				pcs_vector[ProcessIndexCompOil]->SetNodeValue(i, indexConcentration_oil, this->NodeData[i]->CompConcInOil_Nodes[j]);

				// check for unphysical value
				if (this->NodeData[i]->CompConcInOil_Nodes[j] < 0. && verbosity > 1)
				{
					std::cout << " WARNING at Node: " << i << " Pressure: " << m_pcs->GetNodeValue(i, index_pressure2) << " Gas in oil: " <<
						this->NodeData[i]->CompConcInOil_Nodes[j] << "\n";
					std::cout << "         Error in calculation for comps: " << this->NodeData[i]->CompConcInOil_Nodes[j] << "\n";
				}
			}
		} // end vec_OGS_process_index_comps_gas
	} // END of else clause for more than one phase

	// Interpolate the velocity at the Gauss points
	for (std::size_t i = 0; i < this->Phases.size(); i++)
		m_pcs->CalGPVelocitiesfromECLIPSE(path, m_pcs->Tim->step_current,
		i, this->Phases[i]);

	//#ifdef DEBUG_ECLIPSE_OGS_CONVERTER
	// Add cmake var DEBUG_ECLIPSE_OGS_CONVERTER:BOOL = OFF
	//debug(this);
	//#endif    // DEBUG_ECLIPSE_OGS_CONVERTER

	//// Test Output Elements
	//CFEMesh* m_msh = fem_msh_vector[0]; //SB: ToDo hart gesetzt
	//CECLIPSEBlock* m_block = NULL;
	//if(1 < 0) // SB // WTP: hu?
	//if ((m_pcs->Tim->step_current == 0) || (m_pcs->Tim->step_current == 1) ||
	//    (m_pcs->Tim->step_current == 2) || (m_pcs->Tim->step_current == 5) ||
	//    (m_pcs->Tim->step_current == 10) || (m_pcs->Tim->step_current == 20) ||
	//    (m_pcs->Tim->step_current == 30) || (m_pcs->Tim->step_current == 40) ||
	//    (m_pcs->Tim->step_current == 50) || (m_pcs->Tim->step_current == 100) ||
	//    (m_pcs->Tim->step_current == 999) || (m_pcs->Tim->step_current == 1000))
	//{
	//    vector <string> vec_string;
	//    std::string tempstring;
	//    ostringstream temp;
	//    MeshLib::CElem* elem = NULL;
	//    Math_Group::vec<MeshLib::CNode*> Nodes(8);
	//    double Val;
	//    if (this->Phases.size() == 1)
	//    {
	//        int nidx1 = m_pcs->GetNodeValueIndex("PRESSURE1") + 1; //+1... new time level
	//        vec_string.push_back("Element; X; Y; Z; P-Geosys");
	//        for(unsigned long i = 0; i <  m_msh->ele_vector.size(); i++)
	//        {
	//            elem = m_msh->ele_vector[i];
	//            m_block = this->eclgrid[i];
	//            //Calculate the average of the pressure for all nodes of 1 element
	//            Val = 0;
	//            elem->GetNodes(Nodes);
	//            for (long j = 0; j < long(Nodes.Size()); j++)
	//                Val +=m_pcs->GetNodeValue(Nodes[j]->GetIndex(),nidx1) / Nodes.Size();

	//            temp.str("");
	//            temp.clear();
	//            temp << i;
	//            tempstring = temp.str();
	//            temp.str("");
	//            temp.clear();
	//            temp << m_block->x_barycentre;
	//            tempstring += "; " + temp.str();
	//            temp.str("");
	//            temp.clear();
	//            temp << m_block->y_barycentre;
	//            tempstring += "; " + temp.str();
	//            temp.str("");
	//            temp.clear();
	//            temp << m_block->z_barycentre;
	//            tempstring += "; " + temp.str();
	//            temp.str("");
	//            temp.clear();
	//            temp.precision(12);
	//            temp << Val;
	//            tempstring += "; " + temp.str();
	//            vec_string.push_back(tempstring);
	//        }
	//    }
	//    if (this->Phases.size() > 1)
	//    {
	//        int nidx1 = m_pcs->GetNodeValueIndex("PRESSURE1") + 1; //+1... new time level
	//        int nidx2 = m_pcs->GetNodeValueIndex("PRESSURE2") + 1;
	//        int nidx3 = m_pcs->GetNodeValueIndex("SATURATION1") + 1;
	//        vec_string.push_back(
	//                "Element; X; Y; Z; P-Geosys_phase1; P-Geosys_phase2; Saturation_Geosys_phase1");
	//        for(unsigned long i = 0; i <  m_msh->ele_vector.size(); i++)
	//        {
	//            elem = m_msh->ele_vector[i];
	//            m_block = this->eclgrid[i];
	//            temp.str("");
	//            temp.clear();
	//            temp << i;
	//            tempstring = temp.str();
	//            temp.str("");
	//            temp.clear();
	//            temp << m_block->x_barycentre;
	//            tempstring += "; " + temp.str();
	//            temp.str("");
	//            temp.clear();
	//            temp << m_block->y_barycentre;
	//            tempstring += "; " + temp.str();
	//            temp.str("");
	//            temp.clear();
	//            temp << m_block->z_barycentre;
	//            tempstring += "; " + temp.str();
	//            Val = 0;
	//            elem->GetNodes(Nodes);
	//            for (long j = 0; j < long(Nodes.Size()); j++)
	//                Val +=m_pcs->GetNodeValue(Nodes[j]->GetIndex(),nidx1) / Nodes.Size();
	//            temp.str("");
	//            temp.clear();
	//            temp.precision(12);
	//            temp << Val;
	//            tempstring += "; " + temp.str();

	//            Val = 0;
	//            for (long j = 0; j < long(Nodes.Size()); j++)
	//                Val +=m_pcs->GetNodeValue(Nodes[j]->GetIndex(),nidx2) / Nodes.Size();
	//            temp.str("");
	//            temp.clear();
	//            temp.precision(12);
	//            temp << Val;
	//            tempstring += "; " + temp.str();

	//            Val = 0;
	//            for (long j = 0; j < long(Nodes.Size()); j++)
	//                Val +=m_pcs->GetNodeValue(Nodes[j]->GetIndex(),nidx3) / Nodes.Size();
	//            temp.str("");
	//            temp.clear();
	//            temp.precision(12);
	//            temp << Val;
	//            tempstring += "; " + temp.str();

	//            vec_string.push_back(tempstring);
	//        }
	//    }
	//    //CWriteTextfiles *TextFile;
	//    //TextFile = new CWriteTextfiles;
	//    //int position = int(path.find_last_of("\\"));
	//    //std::string path_new;
	//    //path_new = path.substr(0,position);
	//    //position = int(path_new.find_last_of("\\"));
	//    //path_new = path_new.substr(0,position);
	//    //temp.str(""); temp.clear(); temp << m_pcs->Tim->step_current; tempstring = temp.str();
	//    //TextFile->Write_Text(path_new + "\\ElementValues_Geosys_" + tempstring + ".csv", vec_string);
	//} // End Test Output

	const clock_t finish = clock();
	const double time = (static_cast<double>(finish)-static_cast<double>(start)) / CLOCKS_PER_SEC;
	if (verbosity > 2)
		std::cout << "                                              Time: " << time << " seconds." << "\n";
}
/*-------------------------------------------------------------------------
   GeoSys - Function: ExecuteEclipse
   Task: Starts Eclipse
   Return: nothing
   Programming: 09/2009 BG
   Modification: 04/2014 WTP: streamlined the function, moved the initialization bits to other functions
   -------------------------------------------------------------------------*/
void CECLIPSEData::ExecuteEclipse(CReadTextfiles_ECL* eclDataFile, CReadTextfiles_ECL* eclFFile, long Timestep, CRFProcess* m_pcs)
{
	//bool Water_phase_exists = false, Oil_phase_exists = false, Gas_phase_exists = false;
	//int position;
	double timestep_length;
	std::string tempstring;
	std::string EclipseExe; // = "c:\\programme\\ecl\\2008.1\\bin\\pc\\eclipse.exe";
	/*WTP std::string projectname;
	std::string Filename;*/                 
	std::string Keyword;
	//WTP std::string Keyword_well;
	ostringstream temp; //WTP pretty time consuming constructor...
	vector <string> vecData;
	std::string root_folder, geosys_folder;
	//WTP CReadTextfiles_ECL* TextFile;
	std::string file_id = "_SAVE";
	
	clock_t start, finish;
	double time;

	EclipseExe = m_pcs->simulator_path; // path to eclipse executable

	//check if file exists
	if (!UsePrecalculatedFiles)
	if (!UseEclrun)
	if (CheckIfFileExists(EclipseExe) == false)
	{
		std::cout << " ERROR: The ECLIPSE executable could not be found! (" << EclipseExe << ")" << "\n";
        std::cout << flush;
		exit(0);
	};

	if (UsePrecalculatedFiles) // file id# for jenkins bm
	{
		int timestep = aktueller_zeitschritt + timestep_adjust_initial;
		temp.str(""); temp.clear(); temp << timestep;
		if (timestep < 10)
			file_id += "000" + temp.str();
		else if (timestep < 100)
			file_id += "00" + temp.str() ;
		else if (timestep < 1000)
			file_id += "0" + temp.str() ;
		else
			file_id += temp.str();
	}

	// Update the well data in the *.data file if an external *.well file is available
	//if(this->ecl_well.size() != 0)
	if (this->existWells == true)
	if (this->ReplaceWellRate(eclDataFile) == false)
	{
		std::cout << " ERROR: Replacing the well rates did not work correctly !" << "\n";
		exit(0);
	};

	//Read length of current timestep and recalculate it to days
	if (m_pcs->Tim->time_unit == "DAY")
		timestep_length = m_pcs->Tim->time_step_length;
	else
	{
		if (m_pcs->Tim->time_unit == "MINUTE")
			timestep_length = m_pcs->Tim->time_step_length / 60 / 24;
		else
		{
			if (m_pcs->Tim->time_unit == "SECOND")
				timestep_length = m_pcs->Tim->time_step_length / 60 / 60 / 24;
			else
			{
				std::cout << " ERROR: time unit " << m_pcs->Tim->time_unit << " is not yet considered." << "\n";
                std::cout << flush;
				exit(0);
			}
		}
	};
	//Write the timestep in the virtual eclipse file
	//if (this->UsePrecalculatedFiles == false || (this->UsePrecalculatedFiles == true && this->UseSaveEclipseDataFiles == true))
	//{
		Keyword = "TSTEP";
		tempstring = "1*";
		temp.str("");
		temp.clear();
		temp << timestep_length;
		tempstring = tempstring + temp.str();
		tempstring = tempstring + " /";
		vecData.clear();
		vecData.push_back(tempstring);
		if (ReplaceASectionInData(eclDataFile, Keyword, vecData, true) == false)
		{
			std::cout << " ERROR: Replacing a section in the file: " << "\n";
			std::cout << " " << pathECLProject + ".DATA" << "\n";
			std::cout << " didn't work for Keyword TSTEP!" << "\n";
            std::cout << flush;
			//system("Pause"); // SB redo wtp
			exit(0);
		}
		vecData.clear();
	//}

	if (Timestep == 1 && m_pcs->iter_outer_cpl == 0)
	{
		//Rewrite the .data file, *.Fxx file not needed since it is the first timestep
		CWriteTextfiles_ECL* OutputFile;
		OutputFile = new CWriteTextfiles_ECL;
		if (UsePrecalculatedFiles == false)
			OutputFile->Write_Text(pathECLProject + ".DATA", eclDataFile->Data);
		else
			OutputFile->Write_Text(pathECLProject + ".DATA" + file_id, eclDataFile->Data);

		if (this->UsePrecalculatedFiles == false)
		{
			//Run Eclipse the first time from the original *.data file
			// make string for external program call to ECLIPSE
			start = clock();
			int number_loops = 0;
			int maximum_loops = 10;
			do {
				if (verbosity > 1)
				{
					std::cout << "      " << number_loops + 1 << ". trial" << "\n";
					std::cout.flush();    // WTP
				}
				if (this->Windows_System == true)
				{
					if (verbosity < 1)
						tempstring = EclipseExe + " " + pathECLProject + " >NUL";
					else if (verbosity < 2)
						tempstring = EclipseExe + " " + pathECLProject + " > " + "ecl_t" + AddZero(Timestep, 4, true) + ".log";
					else
						tempstring = EclipseExe + " " + pathECLProject;
					// wtp debug
					//std::cout << tempstring << "\n";
					if (system(tempstring.c_str()))
					{
						DisplayMsgLn(" ERROR: Eclipse doesn't run properly!!! ");
                        std::cout << flush;
						exit(0);
					}
				}
				else
				{
					//tempstring = EclipseExe + " ." + pathECLProject.substr(root_folder.length(), root_folder.length()); // SB redo WTP
					tempstring = EclipseExe + " " + pathECLProject; 
					if (system(tempstring.c_str()))
					{
						DisplayMsgLn(" ERROR: Eclipse doesn't run properly!!! ");
                        std::cout << flush;
						exit(0);
					}
				}
				this->pathECLFFile = this->pathECLProject + ".F" + AddZero(Timestep + timestep_adjust_initial, 4, true);
				//std::cout << "path to F-File: " << pathECLFFile << "\n";
				number_loops += 1;
			} while ((CheckIfFileExists(pathECLFFile) == false) && (number_loops <= maximum_loops));

			if (number_loops > maximum_loops)
			{
				std::cout << " ERROR: The Eclipse execution does not work after " << number_loops <<
					" trials!" << "\n";
				exit(0);
			}

			finish = clock();
			time = (static_cast<double>(finish)-static_cast<double>(start)) / CLOCKS_PER_SEC;
			std::cout.flush();    // WTP: for clearer output
			if (verbosity > 2)
			{
				std::cout << "\n";
				//WTP std::cout << "  Timestep: " << Timestep << "\n";
				std::cout << "        ExecuteEclipse() called               Time: " << time << " seconds." << "\n";
			}
		}
	}
	else
	{
		// change Restart file
		// Changes are included in the restart File, because there it is possible to define new porosities ore permeabilties
		// every time at least the number of the time step which is used for restart has to be changed

		// increase the number of time steps used for restart
		Keyword = "RESTART";
		//SB redo WTP tempstring = folder + "TemporaryResults ";
		tempstring = "'" + this->pathECLFolder;
		// For the path to the restart file, ECLIPSE under LINUX wants "\" instead of "/"; replace for this section in folder string
		replace(tempstring.begin(), tempstring.end(), '/', '\\');
		tempstring += "TEMPORARYRESULTS";
		if (this->Windows_System == false) { //CB 1012
			replace(tempstring.begin(), tempstring.end(), '\\', '/');
			//WTP tempstring = "TEMPORARYRESULTS ";
		}
		tempstring += "' ";
		temp.str("");
		temp.clear();
		if (m_pcs->Iterative_Eclipse_coupling == true)
			temp << this->timestep_adjust_iteration_tot - 1 + timestep_adjust_initial;
		else 
			temp << Timestep - 1 + timestep_adjust_initial;
		tempstring = tempstring + temp.str();
		tempstring = tempstring + " /";
		//tempstring = tempstring + " SAVE FORMATTED /";
		vecData.clear();
		vecData.push_back(tempstring);
		//if (this->UsePrecalculatedFiles == false || (this->UsePrecalculatedFiles == true && this->UseSaveEclipseDataFiles == true))
		if (ReplaceASectionInData(eclDataFile, Keyword, vecData, true) == false)
		{
			std::cout << " ERROR: Replacing a section in the file: " << "\n";
			std::cout << " " << pathECLProject + ".DATA" << "\n";
			std::cout << " didn't work for Keyword RESTART!" << "\n";
            std::cout << flush;
			exit(0);
		}

		//Rewrite the .data file and the *FXX file
		CWriteTextfiles_ECL* OutputFile;
		OutputFile = new CWriteTextfiles_ECL;
		CWriteTextfiles_ECL* OutputFile2;
		OutputFile2 = new CWriteTextfiles_ECL;

		if (UsePrecalculatedFiles == false)
		{ 
			OutputFile->Write_Text(pathECLProject + ".DATA", eclDataFile->Data);
			OutputFile2->Write_Text(this->pathECLFFile, eclFFile->Data);
		}
		else
		{
			OutputFile->Write_Text(pathECLProject + ".DATA" + file_id, eclDataFile->Data);
			OutputFile2->Write_Text(this->pathECLFFile + file_id, eclFFile->Data);
		}

		////Release the memory, also helps to prevent access to old data
		delete eclFFile;
		delete eclDataFile;

		// Call ECLIPSE
		if (this->UsePrecalculatedFiles == false)
		{
			//Run Eclipse the first time from the original *.data file
			// make string for external program call to ECLIPSE
			start = clock();
			int number_loops = 0;
			int maximum_loops = 10;
			do {
				if (verbosity > 2)
				{
					std::cout << "      " << number_loops + 1 << ". trial" << "\n";
					std::cout.flush();    // WTP
				}
				if (this->Windows_System == true)
				{
					if (verbosity < 1)
						tempstring = EclipseExe + " " + pathECLProject + " >NUL";
					else if (verbosity < 2)
						tempstring = EclipseExe + " " + pathECLProject + " > " + "ecl_t" + AddZero(Timestep, 4, true) + ".log";
					else
						tempstring = EclipseExe + " " + pathECLProject;
					// wtp debug
					//std::cout << tempstring << "\n";
					if (system(tempstring.c_str()))
					{
						DisplayMsgLn(" ERROR: Eclipse doesn't run properly!!! ");
                        std::cout << flush;
						exit(0);
					}
				}
				else
				{
					//position = static_cast<int>(pathECLProject.find_last_of("/")); // WTP
					//root_folder = pathECLProject.substr(0, position);
					//position = static_cast<int>(root_folder.find_last_of("/"));
					//tempstring = EclipseExe + " ." + pathECLProject.substr(position, position);  //CB 1012
					tempstring = EclipseExe + " " + pathECLProject;
					if (system(tempstring.c_str()))
					{
						DisplayMsgLn(" ERROR: Eclipse doesn't run properly!!! ");
                        std::cout << flush;
						exit(0);
					}
				}
				if (m_pcs->Iterative_Eclipse_coupling == true)
					this->pathECLFFile = this->pathECLProject + ".F" + AddZero(this->timestep_adjust_iteration_tot + timestep_adjust_initial, 4, true);
				else
					this->pathECLFFile = this->pathECLProject + ".F" + AddZero(Timestep + timestep_adjust_initial, 4, true);
				//std::cout << "path to F-File: " << pathECLFFile << "\n";
				number_loops += 1;
			} while ((CheckIfFileExists(pathECLFFile) == false) && (number_loops <= maximum_loops));

			if (number_loops > maximum_loops)
			{
				std::cout << " ERROR: The Eclipse execution does not work after " << number_loops <<
					" trials!" << "\n";
				exit(0);
			}

			finish = clock();
			time = (static_cast<double>(finish)-static_cast<double>(start)) / CLOCKS_PER_SEC;
			std::cout.flush();    // WTP mitigates the overlapping output between OGS and ECL
			if (verbosity > 2)
			{
				std::cout << "\n";
				//WTP std::cout << "  Timestep: " << Timestep << "\n";
				std::cout << "        ExecuteEclipse() called               Time: " << time << " seconds." << "\n";
			}
		}
	} // end of Timestep > 1 condition
}

/*-------------------------------------------------------------------------
   GeoSys - Function: CleanUpEclipseFiles
   Task: Rename, Copy and Delete Files from the Eclipse Run
   Return: nothing
   Programming: 09/2009 BG
   Modification: 04/2014 WTP minor changes
   -------------------------------------------------------------------------*/
bool CECLIPSEData::CleanUpEclipseFiles(std::string folder, std::string projectname, long Timestep, CRFProcess* m_pcs)
{
	std::string systemcommand;
	std::string system_delete;
	std::string system_noquery;
	ostringstream temp;
	std::string extension;
	std::string laenge;

	double time;
	clock_t start, finish;
	start = clock();

	if (verbosity > 2)
		std::cout << "        CleanUpEclipseFiles()";

	if (this->Windows_System == true)
	{
		system_delete = "del ";
		system_noquery = " /Q";
	}
	else
	{// Linux system
		system_delete = "rm ";
		system_noquery = " ";
	}

	// Rename files for ECLIPSE restart
	//delete temporary result files
	systemcommand = system_delete + folder + "TEMPORARYRESULTS.*";
	if (system(systemcommand.c_str()))
			DisplayMsgLn(" WARNING: Could not delete temporary result files! ");

	/*  systemcommand = system_delete + projectname + ".FSAVE" + system_noquery;
	if (system(systemcommand.c_str())){
	// DisplayMsgLn("Could not delete result .FSAVE files! ");
	//SB ?? why return return 0;
	} */

	//if (UseSaveEclipseDataFiles == true)
	//{
		extension = ".F";
		if (m_pcs->Iterative_Eclipse_coupling == true) //KB0116
		{
			int timestep = this->timestep_adjust_iteration_tot + timestep_adjust_initial;
		temp.str(""); temp.clear(); temp << timestep;
		if  (timestep < 10)
			extension += "000" + temp.str() + " ";
		else if (timestep < 100)
			extension += "00" + temp.str() + " ";
		else if (timestep < 1000)
			extension += "0" + temp.str() + " ";
		else
			extension += temp.str() + " ";
		}
		else
		{
			int timestep = aktueller_zeitschritt + timestep_adjust_initial;
			temp.str(""); temp.clear(); temp << timestep;
			if (timestep < 10)
				extension += "000" + temp.str() + " ";
			else if (timestep < 100)
				extension += "00" + temp.str() + " ";
			else if (timestep < 1000)
				extension += "0" + temp.str() + " ";
			else
 				extension += temp.str() + " ";
		}
	//copy original result files
	if (this->Windows_System == true)
	{
		systemcommand = "ren " + projectname + extension + "TEMPORARYRESULTS.F*";
		//// wtp debug
		//std::cout << systemcommand << "\n";
		if (system(systemcommand.c_str()))
			DisplayMsgLn(" WARNING: Could not rename the temporary result files .F*! ");
		systemcommand = "ren " + projectname + ".FGRID " + "TEMPORARYRESULTS.F*";
		if (system(systemcommand.c_str()))
			DisplayMsgLn(" WARNING: Could not rename the temporary grid file! ");
	}
	else
	{
		systemcommand = "cp " + projectname + extension + folder + "TEMPORARYRESULTS" + extension;
		//// wtp debug
		//std::cout << systemcommand << "\n";
		if (system(systemcommand.c_str()))
			DisplayMsgLn(" WARNING: Could not rename the temporary result files .F*! ");
		systemcommand = "cp " + projectname + ".FGRID " + folder + "TEMPORARYRESULTS.FGRID";
		//// wtp debug
		//std::cout << systemcommand << "\n";
		if (system(systemcommand.c_str()))
			DisplayMsgLn(" WARNING: Could not rename the temporary grid file! ");
	}

	// delete all other F* files
	systemcommand = system_delete + projectname + ".F*" + system_noquery;
	if (verbosity < 2)
	{
		if (Windows_System)
			systemcommand += " 2>NUL";
		else
			systemcommand += " > /dev/null";
	}
	//// wtp debug
	//std::cout << systemcommand << "\n";
	if (system(systemcommand.c_str()))
		DisplayMsgLn(" WARNING: Could not delete result files .F*! ");

	systemcommand = system_delete + projectname + ".A*" + system_noquery;
	if (verbosity < 2)
	{
		if (Windows_System)
			systemcommand += " 2>NUL";
		else
			systemcommand += " > /dev/null";
	}
	//// wtp debug
	//std::cout << systemcommand << "\n";
	if (system(systemcommand.c_str()))
			DisplayMsgLn(" WARNING: Could not delete result files .A*! ");


	systemcommand = system_delete + projectname + ".P*" + system_noquery;
	if (verbosity < 2)
	{
		if (Windows_System)
			systemcommand += " 2>NUL";
		else
			systemcommand += " > /dev/null";
	}
	//// wtp debug
	//std::cout << systemcommand << "\n";
	if (system(systemcommand.c_str()))
			DisplayMsgLn(" WARNING: Could not delete result files .P*! ");
	
	systemcommand = system_delete + projectname + ".R*" + system_noquery;
	if (verbosity < 2)
	{
		if (Windows_System)
			systemcommand += " 2>NUL";
		else
			systemcommand += " > /dev/null";
	}
	//// wtp debug
	//std::cout << systemcommand << "\n";
	if (system(systemcommand.c_str()))
			DisplayMsgLn(" WARNING: Could not delete result files .R*! ");

	systemcommand = system_delete + projectname + ".M*" + system_noquery;
	if (verbosity < 2)
	{
		if (Windows_System)
			systemcommand += " 2>NUL";
		else
			systemcommand += " > /dev/null";
	}
	//// wtp debug
	//std::cout << systemcommand << "\n";
	if (system(systemcommand.c_str()))
			DisplayMsgLn(" WARNING: Could not delete result files .M*! ");

	systemcommand = system_delete + projectname + ".DB*" + system_noquery;
	if (verbosity < 2)
	{
		if (Windows_System)
			systemcommand += " 2>NUL";
		else
			systemcommand += " > /dev/null";
	}
	//// wtp debug
	//std::cout << systemcommand << "\n";
	if (system(systemcommand.c_str()))
			DisplayMsgLn(" WARNING: Could not delete result files .DB*! ");


	if (this->UseSaveEclipseDataFiles == true)
	{
		//std::string extension_new = extension.substr(0, extension.size() - 1) + "_SAVE";

		if (this->Windows_System == true)
			systemcommand = "copy " + folder + "TEMPORARYRESULTS" + extension + projectname + "_SAVE" + extension;
		else
			systemcommand = "cp " + folder + "TEMPORARYRESULTS" + extension + projectname + "_SAVE" + extension;
		////// wtp debug
		//std::cout << systemcommand << "\n";
		if (system(systemcommand.c_str()))
			DisplayMsgLn(" WARNING: Could not copy the temporary result files .F*! ");
	}

	finish = clock();
	time = (static_cast<double>(finish)-static_cast<double>(start)) / CLOCKS_PER_SEC;
	if (verbosity > 2)
		std::cout << "                 Time: " << time << " seconds." << "\n";

	return true;
}

/*-------------------------------------------------------------------------
GeoSys - Function: SaveEclipseInputFiles
Task: Save input files of the ecl run for file compare
Return: nothing
Programming: 08/2015 WTP
-------------------------------------------------------------------------*/
bool CECLIPSEData::SaveEclipseInputFiles(std::string folder, std::string projectname)
{
	std::string systemcommand;
	std::string system_delete;
	std::string system_noquery;
	ostringstream temp;
	std::string extension;
	std::string laenge;

	double time;
	clock_t start, finish;
	start = clock();

	if (verbosity > 2)
		std::cout << "        SaveEclipseInputFiles()";

	if (this->Windows_System == true)
	{
		system_delete = "del ";
		system_noquery = " /Q";
	}
	else
	{// Linux system
		system_delete = "rm ";
		system_noquery = " ";
	}

	extension = "_";
	int timestep = aktueller_zeitschritt + timestep_adjust_initial + this->timestep_adjust_iteration_tot;
	temp.str(""); temp.clear(); temp << timestep;
	if (timestep < 10)
		extension += "000" + temp.str() + " ";
	else if (timestep < 100)
		extension += "00" + temp.str() + " ";
	else if (timestep < 1000)
		extension += "0" + temp.str();

	//copy the *.data file
	if (this->Windows_System == true)
	{
		systemcommand = "copy " + folder + "TEMP.INC " + folder + "TEMP" + extension + ".INC";
		systemcommand = "copy " + projectname + ".DATA " + projectname + extension + ".DATA";
	}
	else
	{
		systemcommand = "cp " + folder + "TEMP.INC " + projectname + "_SAVE" + extension;
	}
	////// wtp debug
	//std::cout << systemcommand << "\n";
	if (system(systemcommand.c_str()))
		DisplayMsgLn(" WARNING: Could not copy the temporary result files .F*! ");

	if (this->Windows_System == true)
	{
		systemcommand = "ren " + projectname + extension + "TEMPORARYRESULTS.F*";
		//// wtp debug
		//std::cout << systemcommand << "\n";
		if (system(systemcommand.c_str()))
			DisplayMsgLn(" WARNING: Could not rename the temporary result files .F*! ");
		systemcommand = "ren " + projectname + ".FGRID " + "TEMPORARYRESULTS.F*";
		if (system(systemcommand.c_str()))
			DisplayMsgLn(" WARNING: Could not rename the temporary grid file! ");
	}
	else
	{
		systemcommand = "cp " + projectname + extension + folder + "TEMPORARYRESULTS" + extension;
		//// wtp debug
		//std::cout << systemcommand << "\n";
		if (system(systemcommand.c_str()))
			DisplayMsgLn(" WARNING: Could not rename the temporary result files .F*! ");
		systemcommand = "cp " + projectname + ".FGRID " + folder + "TEMPORARYRESULTS.FGRID";
		//// wtp debug
		//std::cout << systemcommand << "\n";
		if (system(systemcommand.c_str()))
			DisplayMsgLn(" WARNING: Could not rename the temporary grid file! ");
	}

	// delete all other F* files
	systemcommand = system_delete + projectname + ".F*" + system_noquery;
	if (verbosity < 2)
	{
		if (Windows_System)
			systemcommand += " 2>NUL";
		else
			systemcommand += " > /dev/null";
	}
	//// wtp debug
	//std::cout << systemcommand << "\n";
	if (system(systemcommand.c_str()))
		DisplayMsgLn(" WARNING: Could not delete result files .F*! ");
	return true;
}
/*-------------------------------------------------------------------------
   GeoSys - Function: CalculateDeltaGeoSysECL()
   Task: Calculates delta values for changes between ogs and eclipse
   Return: boolean
   Programming: 01/2014 WTP
   Modification:
   -------------------------------------------------------------------------*/
bool CECLIPSEData::CalculateDeltaGeoSysECL(CRFProcess* m_pcs)
{
	const double epsilon = 1e-7;

	CRFProcess* n_pcs = NULL;        // Process: dissovled comps in the liquid phase
	CRFProcess* g_pcs = NULL;        // Process: dissolved comps in the gas phase
	CRFProcess* o_pcs = NULL;        // Process: dissolved comps in the oil phase
	//int indexProcess;
	int indexConcentration_Gas = -99;
	int indexConcentration_Water = -99;
	int indexConcentration_Oil = -99;
	//get no of nodes in mesh
	long nnodes = fem_msh_vector[0]->GetNodesNumber(false);
	int phase1, phase2, phase3 = -1; // SB redo wtp

	// Determine the phase
	switch (static_cast<int>(this->Phases.size()))
	{
	case 1:
		phase1 = 0;		// --> Water (Liquid Flow)
		break;
	case 2:
		phase1 = 0;    // --> Water
		phase2 = 1;    // --> Gas
		break;
	case 3:
		// Assumption that there are 3 phases but water is only used for the boundaries
		// -> the oil and gas phase are the relevant one for the exchange with OGS
		// WTP: Is this still correct?
		//phase1 = 1;
		//phase2 = 2;
		if (this->E100 == true)
		{
			phase1 = 1;
			phase2 = 2;
		}
		else
		{
			phase1 = 0;
			phase2 = 1;
			phase3 = 2;
		}
		break;
	default:
		std::cout << " ERROR: There are not more than 3 phases possible!" << "\n";
		return 0;
		break;
	}

	int saturation1_index = -1;
	int saturation2_index = -1;
	int saturation3_index = -1;
	int pressure1_index = -1; // KB0714: If Only one phase
	int pressure2_index = -1;
	int pressure3_index = -1;

	if (static_cast<int>(this->Phases.size()) > 1)
	{
		// Depending on the no of phase get all the static indices
		saturation1_index = m_pcs->GetNodeValueIndex("SATURATION1") + 1; // +1: new timelevel
		saturation2_index = m_pcs->GetNodeValueIndex("SATURATION2") + 1; // +1: new timelevel
		// get the pressure indicies
		pressure2_index = m_pcs->GetNodeValueIndex("PRESSURE2") + 1; // +1: new timelevel
		if (static_cast<int>(this->Phases.size()) > 2)
		{
			saturation3_index = m_pcs->GetNodeValueIndex("SATURATION3") + 1; // +1: new timelevel
			pressure3_index = m_pcs->GetNodeValueIndex("PRESSURE3") + 1; // +1: new timelevel
		}
    }
	else // KB0714: only one phase
	{
		pressure1_index = m_pcs->GetNodeValueIndex("PRESSURE1") + 1;
	}

	//for (long j = 0; j < nnodes; j++)
	//{
	//	std::vector <double> phase_pressure1_ogs_test;
	//	phase_pressure1_ogs_test.clear();
	//	//double phase_pressure1_ecl = NodeData[j]->phase_pressure[phase1];
	//	phase_pressure1_ogs_test[j] = m_pcs->GetNodeValue(j, pressure1_index);
	//	//std::cout << "Ecl " << phase_pressure1_ecl << " OGS " << phase_pressure1_ogs << "\n";
	//}
	
	for (long j = 0; j < nnodes; j++)
	{
		/* wtp: changed to a uniform approach
		  use these variables:
		  std::vector <std::pair<std::string,int>> vec_OGS_process_index_comps_gas;
		  std::vector <std::pair<std::string,int>> vec_OGS_process_index_comps_water;
		  to get the process indicies. The string being the process name in OGS and the int value the corresponding process index
		  */
		// wtp debug
		//if (j == 5348 || j == 5350 || j == 5349 || j == 5351 || j == 5750 || j == 5752 || j == 5751 || j == 5753)
		//	std::cout << " found node" << "\n";
		// WTP: Run over all components and collect the component concentrations
		for (unsigned int i = 0; i < vec_components_ECL_OGS.size(); i++)
		{
			// store the current process data for the water phase
			std::string OGSProcessNameCompWater = "";
			int OGSProcessIndexCompWater = -99;
			int Current_Comp_Index_Water = -99;

			// store the current process data for the gas phase
			std::string OGSProcessNameCompGas = "";
			int OGSProcessIndexCompGas = -99;
			int Current_Comp_Index_Gas = -99;

			// store the current process data for the oil phase
			std::string OGSProcessNameCompOil = "";
			int OGSProcessIndexCompOil = -99;
			int Current_Comp_Index_Oil = -99;

			// Variables for the molar properties of the component
			//double MolarWeight_Comp = -1.;

			// assign the indices for the water phase
			if (vec_components_ECL_OGS[i][1] == 1)
			{
				for (unsigned int j = 0; j < vec_OGS_process_index_comps_water.size(); j++)
				{
					if (vec_components_ECL_OGS_pcs_names[i][1] == vec_OGS_process_index_comps_water[j].first)
					{
						// store the name and the index of the current (gas) process 
						OGSProcessNameCompWater = vec_OGS_process_index_comps_water[j].first;
						OGSProcessIndexCompWater = vec_OGS_process_index_comps_water[j].second;
						Current_Comp_Index_Water = j;
					}
				}
			}

			// assign the indices for the gas phase
			if (vec_components_ECL_OGS[i][2] == 1)
			{
				for (unsigned int j = 0; j < vec_OGS_process_index_comps_gas.size(); j++)
				{
					if (vec_components_ECL_OGS_pcs_names[i][2] == vec_OGS_process_index_comps_gas[j].first)
					{
						// store the name and the index of the current (gas) process 
						OGSProcessNameCompGas = vec_OGS_process_index_comps_gas[j].first;
						OGSProcessIndexCompGas = vec_OGS_process_index_comps_gas[j].second;
						Current_Comp_Index_Gas = j;
					}
				}
			}

			// assign the indices for the oil phase
			if (vec_components_ECL_OGS[i][3] == 1)
			{
				for (unsigned int j = 0; j < vec_OGS_process_index_comps_oil.size(); j++)
				{
					if (vec_components_ECL_OGS_pcs_names[i][3] == vec_OGS_process_index_comps_oil[j].first)
					{
						// store the name and the index of the current (gas) process 
						OGSProcessNameCompOil = vec_OGS_process_index_comps_oil[j].first;
						OGSProcessIndexCompOil = vec_OGS_process_index_comps_oil[j].second;
						Current_Comp_Index_Oil = j;
					}
				}
			}

			// get index of species concentration in nodevaluevector of this process
			if (OGSProcessIndexCompWater > -99)
			{
				n_pcs = pcs_vector[OGSProcessIndexCompWater];
				indexConcentration_Water = n_pcs->GetNodeValueIndex(n_pcs->pcs_primary_function_name[0]) + 1; // +1: new timelevel
			}

			if (OGSProcessIndexCompGas > -99)
			{
				g_pcs = pcs_vector[OGSProcessIndexCompGas];
				indexConcentration_Gas = g_pcs->GetNodeValueIndex(g_pcs->pcs_primary_function_name[0]) + 1; // +1: new timelevel
			}

			if (OGSProcessIndexCompOil > -99)
			{
				o_pcs = pcs_vector[OGSProcessIndexCompOil];
				indexConcentration_Oil = o_pcs->GetNodeValueIndex(o_pcs->pcs_primary_function_name[0]) + 1; // +1: new timelevel
			}

			// get the delta value for each component in each phase
			// first do everything for the dissolved stuff (in water)
			if (OGSProcessIndexCompWater == -99)
			{
				continue;
			}
			else
			{
				// first get the concentration value coming from OGS (e.g. changes due to kinetic reactions)
				double CompConc_OGS = n_pcs->GetNodeValue(j, indexConcentration_Water);
				// secondly get the original concentration value from ECL
				double CompConc_ECL = this->NodeData[j]->CompConcInWater_Nodes[Current_Comp_Index_Water];
				// calcucalte the difference in concentration
				double delta_CompConc = CompConc_OGS - CompConc_ECL;
				// display a warning if concentration in ogs is below zero
				if (CompConc_OGS < 0.0 && verbosity > 1)
				{
					std::cout << " WARNING: Concentration of component " << j << " is negative!" << "\n";
					std::cout << " Value: " << CompConc_OGS << "\n";
				}

				// if the abs. change is significant store it in the desig. vector structure.
				if (fabs(delta_CompConc) > 1.0E-25)
				{
					// if the old value was not zero evaluate if the change is less than epsilon %
					// true -> set delta to zero and neglect the difference in concentration
					if (this->NodeData[j]->CompConcInWater_Nodes[Current_Comp_Index_Water] > 0.0)
						if ((fabs(delta_CompConc) / this->NodeData[j]->CompConcInWater_Nodes[Current_Comp_Index_Water]) < epsilon)
							delta_CompConc = 0.;
				}
				else
					delta_CompConc = 0.0;
				// set delta 
				this->NodeData[j]->delta_CompConcInWater_Nodes[Current_Comp_Index_Water] = delta_CompConc;
			}

			// Now do the same for the gas phase...
			if (OGSProcessIndexCompGas == -99)
			{
				continue;
			}
			else
			{
				// first get the concentration value coming from OGS (e.g. changes due to kinetic reactions)
				double CompConc_OGS = g_pcs->GetNodeValue(j, indexConcentration_Gas);
				// secondly get the original concentration value from ECL
				double CompConc_ECL = this->NodeData[j]->CompConcInGas_Nodes[Current_Comp_Index_Gas];
				// calcucalte the difference in concentration
				double delta_CompConc = CompConc_OGS - CompConc_ECL;
				
				// display a warning if concentration in ogs is below zero
				if(CompConc_OGS < 0.0 && verbosity > 1)
				{
					std::cout << " WARNING: Concentration of component " << j << " is negative!" << "\n";
					std::cout << " Value: " << CompConc_OGS << "\n";	
				}

				// if the abs. change is significant store it in the desig. vector structure.
				if (fabs(delta_CompConc) > 1.0E-25)
				{
					// if the old value was not zero evaluate if the change is less than epsilon %
					// true -> set delta to zero and neglect the difference in concentration
					if (this->NodeData[j]->CompConcInGas_Nodes[Current_Comp_Index_Gas] > 0.0)
						if ((fabs(delta_CompConc) / this->NodeData[j]->CompConcInGas_Nodes[Current_Comp_Index_Gas]) < epsilon)
							delta_CompConc = 0.;
				}
				else
					delta_CompConc = 0.0;
				// set delta 
				this->NodeData[j]->delta_CompConcInGas_Nodes[Current_Comp_Index_Gas] = delta_CompConc;
				
			}
			// and here for the oil phase...
			if (OGSProcessIndexCompOil == -99)
			{
				continue;
			}
			else
			{
				// first get the concentration value coming from OGS (e.g. changes due to kinetic reactions)
				double CompConc_OGS = o_pcs->GetNodeValue(j, indexConcentration_Oil);
				// secondly get the original concentration value from ECL
				double CompConc_ECL = this->NodeData[j]->CompConcInOil_Nodes[Current_Comp_Index_Oil];
				// calcucalte the difference in concentration
				double delta_CompConc = CompConc_OGS - CompConc_ECL;
				// display a warning if concentration in ogs is below zero
				if (CompConc_OGS < 0.0 && verbosity > 1)
				{
					std::cout << " WARNING: Concentration of component " << j << " is negative!" << "\n";
					std::cout << " Value: " << CompConc_OGS << "\n";
				}

				// if the abs. change is significant store it in the desig. vector structure.
				if (fabs(delta_CompConc) > 1.0E-25)
				{
					// if the old value was not zero evaluate if the change is less than epsilon %
					// true -> set delta to zero and neglect the difference in concentration
					if (this->NodeData[j]->CompConcInOil_Nodes[Current_Comp_Index_Oil] > 0.0)
						if ((fabs(delta_CompConc) / this->NodeData[j]->CompConcInOil_Nodes[Current_Comp_Index_Oil]) < epsilon)
							delta_CompConc = 0.;
				}
				else
					delta_CompConc = 0.0;
				// set delta 
				this->NodeData[j]->delta_CompConcInOil_Nodes[Current_Comp_Index_Oil] = delta_CompConc;
			}
		}// end of component loop

		// Now get the gas saturation done
        if(static_cast<int>(this->Phases.size()) == 2 || (this->E100 == true && this->Phases.size() > 1))
		{
			// get the node values old and new
			double sat2_ecl = NodeData[j]->phase_saturation[phase2];
			double sat2_ogs = 1.0 - m_pcs->GetNodeValue(j, saturation1_index);
			// Check consistency
			if ((sat2_ogs < 0) || (sat2_ogs > 1))
			{
				std::cout << " ERROR: The gas saturation after the reactions is beyond physical limts: " << "\n";
				std::cout << " Node: " << j << " Sgas_old: " << sat2_ecl << " Sgas_new: " << sat2_ogs << "\n";
				return 0;
			}
			if ((sat2_ecl < 0) || (sat2_ecl > 1))
			{
				std::cout << " ERROR: The gas saturation from Eclipse is beyond physical limts: " << "\n";
				std::cout << " Node: " << j << " Sgas_old: " << sat2_ecl << " Sgas_new: " << sat2_ogs << "\n";
				return 0;
			}
			//calculate change of Sgas
			NodeData[j]->deltaSatGas = sat2_ogs - sat2_ecl;

			if (NodeData[j]->deltaSatGas != 0)
			{
				//neglect difference if it is too small
				if ((fabs(NodeData[j]->deltaSatGas) / sat2_ecl) < epsilon)
					NodeData[j]->deltaSatGas = 0;
				else // set the new Sgas sat2_ogs to NodeData
					NodeData[j]->phase_saturation[phase2] = sat2_ogs;
			}
			if ((NodeData[j]->phase_saturation[phase2] < 0) || (NodeData[j]->phase_saturation[phase2] > 1))
			{
				std::cout << " Sgas_new: " << NodeData[j]->phase_saturation[phase2] << " deltaSat: " << NodeData[j]->deltaSatGas << "\n";
				std::cout << "  Error in calculation of Sgas: " << NodeData[j]->phase_saturation[phase2] << "\n";
			}
		}
		else if (static_cast<int>(this->Phases.size()) == 3 && this->E100 != true)
		{
			// get the node values old and new
			double sat2_ecl = NodeData[j]->phase_saturation[phase2];
			double sat3_ecl = NodeData[j]->phase_saturation[phase3];
			double sat2_ogs = m_pcs->GetNodeValue(j, saturation2_index);
			double sat3_ogs = m_pcs->GetNodeValue(j, saturation3_index);
			// Check consistency of data
			if ((sat2_ogs < 0) || (sat2_ogs > 1))
			{
				std::cout << " ERROR: The gas saturation after the reactions is beyond physical limts: " << "\n";
				std::cout << " Node: " << j << " Sgas_old: " << sat2_ecl << " Sgas_new: " << sat2_ogs << "\n";
				return 0;
			}
			if ((sat2_ecl < 0) || (sat2_ecl > 1))
			{
				std::cout << " ERROR: The gas saturation from Eclipse is beyond physical limts: " << "\n";
				std::cout << " Node: " << j << " Sgas_old: " << sat2_ecl << " Sgas_new: " << sat2_ogs << "\n";
				return 0;
			}
			if ((sat3_ogs < 0) || (sat3_ogs > 1))
			{
				std::cout << " ERROR: The oil saturation after the reactions is beyond physical limts: " << "\n";
				std::cout << " Node: " << j << " Soil_old: " << sat3_ecl << " Soil_new: " << sat3_ogs << "\n";
				return 0;
			}
			if ((sat3_ecl < 0) || (sat3_ecl > 1))
			{
				std::cout << " ERROR: The oil saturation from Eclipse is beyond physical limts: " << "\n";
				std::cout << " Node: " << j << " Soil_old: " << sat3_ecl << " Soil_new: " << sat3_ogs << "\n";
				return 0;
			}
			//calculate change of Sgas
			NodeData[j]->deltaSatGas = sat2_ogs - sat2_ecl;
			if (NodeData[j]->deltaSatGas != 0)
			{
				//neglect difference if it is too small
				if ((fabs(NodeData[j]->deltaSatGas) / sat2_ecl) < epsilon)
					NodeData[j]->deltaSatGas = 0;
				else // set the new Sgas sat2_ogs to NodeData
					NodeData[j]->phase_saturation[phase2] = sat2_ogs;
			}
			if ((NodeData[j]->phase_saturation[phase2] < 0) || (NodeData[j]->phase_saturation[phase2] > 1))
			{
				std::cout << " Sgas_new: " << NodeData[j]->phase_saturation[phase2] << " deltaSat: " << NodeData[j]->deltaSatGas << "\n";
				std::cout << " Error in calculation of Sgas: " << NodeData[j]->phase_saturation[phase2] << "\n";
			}
			//calculate change of Soil
			NodeData[j]->deltaSatOil = sat3_ogs - sat3_ecl;
			if (NodeData[j]->deltaSatOil != 0)
			{
				//neglect difference if it is too small
				if ((fabs(NodeData[j]->deltaSatOil) / sat3_ecl) < epsilon)
					NodeData[j]->deltaSatOil = 0;
				else // set the new Sgas sat2_ogs to NodeData
					NodeData[j]->phase_saturation[phase3] = sat3_ogs;
			}
			if ((NodeData[j]->phase_saturation[phase3] < 0) || (NodeData[j]->phase_saturation[phase3] > 1))
			{
				std::cout << " Soil_new: " << NodeData[j]->phase_saturation[phase3] << " deltaSat: " << NodeData[j]->deltaSatOil << "\n";
				std::cout << " Error in calculation of Soil: " << NodeData[j]->phase_saturation[phase3] << "\n";
			}
		} // end of saturation calculations

		// Now calcualte changes in phase pressures
		if (static_cast<int>(this->Phases.size()) == 2 || (this->E100 == true && this->Phases.size() > 1))
		{
			// CB Return delta Pressure
			// get the node values old and new
			double phase_pressure2_ecl = NodeData[j]->phase_pressure[phase2];
			double phase_pressure2_ogs = m_pcs->GetNodeValue(j, pressure2_index);
			// Check consistency
			if (phase_pressure2_ogs < 0)
			{
				std::cout << " ERROR: The gas pressure after the reactions is beyond physical limts: " << "\n";
				std::cout << " Node: " << j << " PRESSURE2_old: " << phase_pressure2_ecl << " PRESSURE2_new: " << phase_pressure2_ogs << "\n";
				return 0;
			}
			if (phase_pressure2_ecl < 0)
			{
				std::cout << " ERROR: The gas pressure from Eclipse is beyond physical limts: " << "\n";
				std::cout << " Node: " << j << " PRESSURE2_old: " << phase_pressure2_ecl << " PRESSURE2_new: " << phase_pressure2_ogs << "\n";
				return 0;
			}

			//calculate change of PRESSURE2
			NodeData[j]->deltaPress2 = phase_pressure2_ogs - phase_pressure2_ecl;
			if (NodeData[j]->deltaPress2 != 0)
			{
				//neglect difference if it is too small
				if ((fabs(NodeData[j]->deltaPress2) / phase_pressure2_ecl) < epsilon)
				{
					NodeData[j]->deltaPress2 = 0;
				}
				else // set the new Sgas phase_pressure2_ogs to NodeData
					NodeData[j]->phase_pressure[phase2] = phase_pressure2_ogs;
			}
			if (NodeData[j]->phase_pressure[phase2] < 0)
			{
				std::cout << " PRESSURE2_new: " << NodeData[j]->phase_pressure[phase2] << " deltaPress2: " << NodeData[j]->deltaPress2 << "\n";
				std::cout << " Error in calculation of PRESSURE2: " << NodeData[j]->phase_pressure[phase2] << "\n";
			}
		}
		else if (static_cast<int>(this->Phases.size()) == 3 && this->E100 != true)
		{
			// get the node values old and new
			double phase_pressure2_ecl = NodeData[j]->phase_pressure[phase2];
			double phase_pressure2_ogs = m_pcs->GetNodeValue(j, pressure2_index);
			double phase_pressure3_ecl = NodeData[j]->phase_pressure[phase3];
			double phase_pressure3_ogs = m_pcs->GetNodeValue(j, pressure3_index);
			// Check consistency
			if (phase_pressure2_ogs < 0)
			{
				std::cout << " ERROR: The gas pressure after the reactions is beyond physical limts: " << "\n";
				std::cout << " Node: " << j << " PRESSURE2_old: " << phase_pressure2_ecl << " PRESSURE2_new: " << phase_pressure2_ogs << "\n";
				return 0;
			}
			if (phase_pressure2_ecl < 0)
			{
				std::cout << " ERROR: The gas pressure from Eclipse is beyond physical limts: " << "\n";
				std::cout << " Node: " << j << " PRESSURE2_old: " << phase_pressure2_ecl << " PRESSURE2_new: " << phase_pressure2_ogs << "\n";
				return 0;
			}
			if (phase_pressure3_ogs < 0)
			{
				std::cout << " ERROR: The  oil pressure after the reactions is beyond physical limts: " << "\n";
				std::cout << " Node: " << j << " PRESSURE3_old: " << phase_pressure3_ecl << " PRESSURE3_new: " << phase_pressure3_ogs << "\n";
				return 0;
			}
			if (phase_pressure3_ecl < 0)
			{
				std::cout << " ERROR: The oil pressure from Eclipse is beyond physical limts: " << "\n";
				std::cout << " Node: " << j << " PRESSURE3_old: " << phase_pressure3_ecl << " PRESSURE3_new: " << phase_pressure3_ogs << "\n";
				return 0;
			}

			//calculate change of PRESSURE2
			NodeData[j]->deltaPress2 = phase_pressure2_ogs - phase_pressure2_ecl;
			if (NodeData[j]->deltaPress2 != 0)
			{
				//neglect difference if it is too small
				if ((fabs(NodeData[j]->deltaPress2) / phase_pressure2_ecl) < epsilon)
				{
					NodeData[j]->deltaPress2 = 0;
				}
				else // set the new Sgas phase_pressure2_ogs to NodeData
					NodeData[j]->phase_pressure[phase2] = phase_pressure2_ogs;
			}
			if (NodeData[j]->phase_pressure[phase2] < 0)
			{
				std::cout << " PRESSURE2_new: " << NodeData[j]->phase_pressure[phase2] << " deltaPress2: " << NodeData[j]->deltaPress2 << "\n";
				std::cout << " Error in calculation of PRESSURE2: " << NodeData[j]->phase_pressure[phase2] << "\n";
			}
			//calculate change of PRESSURE3
			NodeData[j]->deltaPress3 = phase_pressure3_ogs - phase_pressure3_ecl;
			if (NodeData[j]->deltaPress3 != 0)
			{
				//neglect difference if it is too small
				if ((fabs(NodeData[j]->deltaPress3) / phase_pressure3_ecl) < epsilon)
				{
					NodeData[j]->deltaPress3 = 0;
				}
				else // set the new Sgas phase_pressure2_ogs to NodeData
					NodeData[j]->phase_pressure[phase3] = phase_pressure2_ogs;
			}
			if (NodeData[j]->phase_pressure[phase3] < 0)
			{
				std::cout << " PRESSURE3_new: " << NodeData[j]->phase_pressure[phase3] << " deltaPress3: " << NodeData[j]->deltaPress3 << "\n";
				std::cout << " Error in calculation of PRESSURE3: " << NodeData[j]->phase_pressure[phase3] << "\n";
			}
		}
		else if (static_cast<int>(this->Phases.size()) == 1)
		{ 
			// CB Return delta Pressure
			// get the node values old and new
			double phase_pressure1_ecl = NodeData[j]->phase_pressure[phase1];
			double phase_pressure1_ogs = m_pcs->GetNodeValue(j, pressure1_index);

			// Check consistency
			if (phase_pressure1_ogs < 0)
			{
				std::cout << " ERROR: The gas pressure after the reactions is beyond physical limts: " << "\n";
				std::cout << " Node: " << j << " PRESSURE1_old: " << phase_pressure1_ecl << " PRESSURE1_new: " << phase_pressure1_ogs << "\n";
				return 0;
			}
			if (phase_pressure1_ecl < 0)
			{
				std::cout << " ERROR: The gas pressure from Eclipse is beyond physical limts: " << "\n";
				std::cout << " Node: " << j << " PRESSURE1_old: " << phase_pressure1_ecl << " PRESSURE1_new: " << phase_pressure1_ogs << "\n";
				return 0;
			}

			//calculate change of PRESSURE2
			NodeData[j]->deltaPress1 = phase_pressure1_ogs - phase_pressure1_ecl;
			if (NodeData[j]->deltaPress1 != 0)
			{
				//neglect difference if it is too small
				if ((fabs(NodeData[j]->deltaPress1) / phase_pressure1_ecl) < epsilon)
				{
					NodeData[j]->deltaPress1 = 0;
				}
				else // set the new Sgas phase_pressure2_ogs to NodeData
					NodeData[j]->phase_pressure[phase2] = phase_pressure1_ogs;
			}
			if (NodeData[j]->phase_pressure[phase1] < 0)
			{
				std::cout << " PRESSURE1_new: " << NodeData[j]->phase_pressure[phase1] << " deltaPress1: " << NodeData[j]->deltaPress1 << "\n";
				std::cout << " Error in calculation of PRESSURE1: " << NodeData[j]->phase_pressure[phase1] << "\n";
			}
		}
	} // end of nnodes loop
	return true;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: InterpolateDeltaGeoSysECL()
   Task: Interpolates delta values for changes between ogs and eclipse
   Return: boolean
   Programming: 01/2014 WTP
   Modification:
   -------------------------------------------------------------------------*/
bool CECLIPSEData::InterpolateDeltaGeoSysECL(CRFProcess* m_pcs)
{
	//MeshLib::CElem* m_element = NULL;
	CFEMesh* m_msh = fem_msh_vector[0];
	bool test = true;
	double weight = -1.;
	const double epsilon = 1e-7;
	int variable_index_FVF_Liquid = -1;
	int variable_index_RS = -1;

	// Get all general keyword indices for the Eclipse data
	int variable_index_porevolume = this->GetVariableIndex("RPORV");
	int variable_index_water_saturation = this->GetVariableIndex("SWAT");
	int variable_index_gas_saturation = this->GetVariableIndex("SGAS");   // CB
	int variable_index_oil_saturation = this->GetVariableIndex("SOIL");
	int variable_index_oil_pressure = this->GetVariableIndex("POIL");
	int variable_index_gas_pressure = this->GetVariableIndex("PRESSURE");
	std::vector <int> vec_variable_indicies_Comp_MolarDensity;

	// get all MLSC indicies
	if (this->E100 != true)
	{
		for (unsigned int i = 0; i < this->vec_components_ECL_OGS.size(); i++)
		{
			//WTP: we need them anyway, even if the component is only transported in the gas phase
			if (this->vec_components_ECL_OGS[i][0] != -1) // Check if component is transported in the water phase 
			{
				std::stringstream ss;
				ss << vec_components_ECL_OGS[i][0];
				vec_variable_indicies_Comp_MolarDensity.push_back(GetVariableIndex("MLSC" + ss.str()));
			}
		}
	}
	else
	{
		// get the index for the liquid formation volume factor to convert RS back to surface values
		variable_index_FVF_Liquid = this->GetVariableIndex("BO");
		variable_index_RS = this->GetVariableIndex("RS");
	}

	// start the general element loop
	for (long i = 0; i < this->elements; i++)
	{
		//Modify the scheme to set vec_element_nodes:BW
		Math_Group::vec <MeshLib::CNode*> vec_element_nodes(this->eclgrid[i]->MshNodeindex.size());
		for (size_t tempi = 0; tempi < this->eclgrid[i]->MshNodeindex.size(); tempi++)
			vec_element_nodes[tempi] = m_msh->nod_vector[this->eclgrid[i]->MshNodeindex[tempi]];
		//m_element = m_msh->ele_vector[i];
		//Math_Group::vec <MeshLib::CNode*> vec_element_nodes(m_element->GetNodesNumber(false));	//Nodes number might be less than 8:BW
		//m_element->GetNodes(vec_element_nodes);

		std::vector<double> vec_dummy_mlsc;
		// First interpolate the change in component concentrations
		// WTP: Interpolate the change of component concentration from nodes to elements
		// A loop is needed for the interpolation of the component concentrations from the nodes to the elements
		// Now run over all the processes and collect the component concentrations in water
		int current_run = -1; // WTP: This has to go, local variable to select the mlsc value
		for (unsigned int j = 0; j < vec_components_ECL_OGS.size(); j++)
		{
			if (vec_components_ECL_OGS[j][0] != -1)
			{
				current_run++; // WTP this is nto final...just a variable to acceess the correct mlsc keyword 
				double MolarWeight_Comp = -1.;
				// store the current process data for the water phase
				std::string OGSProcessNameCompWater = "";
				int OGSProcessIndexCompWater = -99;
				int Current_Comp_Index_Water = -99;

				// store the current process data for the gas phase
				std::string OGSProcessNameCompGas = "";
				int OGSProcessIndexCompGas = -99;
				int Current_Comp_Index_Gas = -99;

				// store the current process data for the oil phase
				std::string OGSProcessNameCompOil = "";
				int OGSProcessIndexCompOil = -99;
				int Current_Comp_Index_Oil = -99;

				// first assign all indices
				// assign the indices for the water phase
				if (vec_components_ECL_OGS[j][1] == 1)
				{
					for (unsigned int k = 0; k < vec_OGS_process_index_comps_water.size(); k++)
					{
						if (vec_components_ECL_OGS_pcs_names[j][1] == vec_OGS_process_index_comps_water[k].first)
						{
							// store the name and the index of the current (gas) process 
							OGSProcessNameCompWater = vec_OGS_process_index_comps_water[k].first;
							OGSProcessIndexCompWater = vec_OGS_process_index_comps_water[k].second;
							Current_Comp_Index_Water = k;

							// get the molar properties of the component
							for (unsigned int l = 0; l < cp_vec.size(); l++)
							{
								if (cp_vec[l]->compname == vec_OGS_process_index_comps_water[k].first)
								{
									MolarWeight_Comp = cp_vec[l]->molar_weight;
									break;
								}
							}
						}
					}
				}
				if (this->E100 != true) // Only assign gas and oil indices if it is not E100
				{
					// assign the indices for the gas phase
					if (vec_components_ECL_OGS[j][2] == 1)
					{
						for (unsigned int k = 0; k < vec_OGS_process_index_comps_gas.size(); k++)
						{
							if (vec_components_ECL_OGS_pcs_names[j][2] == vec_OGS_process_index_comps_gas[k].first)
							{
								// store the name and the index of the current (gas) process 
								OGSProcessNameCompGas = vec_OGS_process_index_comps_gas[k].first;
								OGSProcessIndexCompGas = vec_OGS_process_index_comps_gas[k].second;
								Current_Comp_Index_Gas = k;
								if (MolarWeight_Comp == -1.)
								{
									// get the molar properties of the component but only if they are not set yet
									for (unsigned int l = 0; l < cp_vec.size(); l++)
									{
										if (cp_vec[l]->compname == vec_OGS_process_index_comps_gas[k].first)
										{
											MolarWeight_Comp = cp_vec[l]->molar_weight;
											break;
										}
									}
								}
							}
						}
					}
					// assign the indices for the oil phase
					if (vec_components_ECL_OGS[j][3] == 1)
					{
						for (unsigned int k = 0; k < vec_OGS_process_index_comps_oil.size(); k++)
						{
							if (vec_components_ECL_OGS_pcs_names[j][3] == vec_OGS_process_index_comps_oil[k].first)
							{
								// store the name and the index of the current (gas) process 
								OGSProcessNameCompOil = vec_OGS_process_index_comps_oil[k].first;
								OGSProcessIndexCompOil = vec_OGS_process_index_comps_oil[k].second;
								Current_Comp_Index_Oil = k;
								if (MolarWeight_Comp == -1.)
								{
									// get the molar properties of the component but only if they are not set yet
									for (unsigned int l = 0; l < cp_vec.size(); l++)
									{
										if (cp_vec[l]->compname == vec_OGS_process_index_comps_oil[k].first)
										{
											MolarWeight_Comp = cp_vec[l]->molar_weight;
											break;
										}
									}
								}
							}
						}
					}
				}
				if (this->E100 == true && OGSProcessIndexCompWater > -99)
				{
					// reset sdelta gas variable
					double delta_CompConc = 0.0;
					double newRS = 0.;
					// set the FVF for each element
					double FVF_Liquid = this->Data[i][variable_index_FVF_Liquid];
					weight = float(1.0 / vec_element_nodes.Size()); //arithmetic average SB redo wtp

					for (long k = 0; k < long(vec_element_nodes.Size()); k++)
					{
						// get the delta in dissolved gas components from the vector created previously
						double delta_comp_value = this->NodeData[vec_element_nodes[k]->GetIndex()]->delta_CompConcInWater_Nodes[Current_Comp_Index_Water];
						// add it up
						delta_CompConc += delta_comp_value * weight;
					}
					if (delta_CompConc< 0.)
					{
						//if (fabs(delta_RS) > (this->Data[j][this->variable_index_RS] - epsilon))
						if (fabs(delta_CompConc) >(this->vec_CompConc_Water_elements[i][Current_Comp_Index_Water] - epsilon))
						{
							if (verbosity > 1)
							{
								std::cout << " WARNING: The amount of gas removed from water is larger than the available amount of gas! The value is corrected!"
									<< "\n";
								for (long k = 0; k < long(vec_element_nodes.Size()); k++)
									std::cout << "     Node: " << vec_element_nodes[k]->GetIndex() << " Value: "
									<< this->NodeData[vec_element_nodes[k]->GetIndex()]->delta_CompConcInWater_Nodes[Current_Comp_Index_Water] << "\n";
								std::cout << "  Element: " << i << " delta conc: " << delta_CompConc << " Value_old: " <<
									this->vec_CompConc_Water_elements[i][Current_Comp_Index_Water];
							}
							delta_CompConc = -this->vec_CompConc_Water_elements[i][Current_Comp_Index_Water];
						}
					}
					if (delta_CompConc != 0.)
					{
						// if the change is insignificant just take the old value (remember the delta is interpolated...)
						if ((fabs(delta_CompConc) / this->vec_CompConc_Water_elements[i][Current_Comp_Index_Water]) < epsilon)
						{
							// WTP: vecRS needs to be a vector <vector <double>> field i think.....
							newRS = this->vec_CompConc_Water_elements[i][Current_Comp_Index_Water] * FVF_Liquid * MolarWeight_Comp / SurfaceDensity_Gas_E100;
						}
						else
						{
							// if the change is significant, take the old value and add the delta value. 
							newRS = (this->vec_CompConc_Water_elements[i][Current_Comp_Index_Water] + delta_CompConc)
								*  FVF_Liquid * MolarWeight_Comp
								/ SurfaceDensity_Gas_E100;
						}
					}
					else
					{
						// if the delta is 0 take the old value
						newRS = this->vec_CompConc_Water_elements[i][Current_Comp_Index_Water] * FVF_Liquid * MolarWeight_Comp / SurfaceDensity_Gas_E100;
					}
					vec_RS_elements.push_back(newRS);
					if (vec_RS_elements[vec_RS_elements.size() - 1] < 0)
					{
						std::cout << " ERROR: The new calculated amount of dissolved components in gas (RS) is negative!" << "\n";
						return 0; //exit(0) ??
					}
				}    // END of E100 case
				// Begin of the E300 case
				if (this->E100 != true) // && vec_components_ECL_OGS[0][j] != -1)
				{
					double deltaMolesWater = 0.;
					double deltaMolesGas = 0.;
					double deltaMolesOil = 0.;
					double RPORV = this->Data[i][variable_index_porevolume];

					if (OGSProcessIndexCompWater != -99)
					{
						// reset delta conc variable
						double delta_CompConc = 0.0;
						weight = float(1.0 / vec_element_nodes.Size()); //arithmetic average SB redo wtp
						for (long k = 0; k < long(vec_element_nodes.Size()); k++)
						{
							// get the delta in water components & add it up
							double delta_comp_value = this->NodeData[vec_element_nodes[k]->GetIndex()]->delta_CompConcInWater_Nodes[Current_Comp_Index_Water];
							delta_CompConc += delta_comp_value * weight;
						}
						if (delta_CompConc < 0.0 && fabs(delta_CompConc) > this->vec_CompConc_Water_elements[i][Current_Comp_Index_Water])
						{
							if (verbosity > 1)
							{
								std::cout << " WARNING: The amount of component removed from water is larger than the available amount! The value is corrected!"
									<< "\n";
								for (long k = 0; k < long(vec_element_nodes.Size()); k++)
								{
									std::cout << "     Node: " << vec_element_nodes[k]->GetIndex() << " Value: "
										<< this->NodeData[vec_element_nodes[k]->GetIndex()]->delta_CompConcInWater_Nodes[Current_Comp_Index_Water] << "\n";
								}
								std::cout << "  Element: " << i << " delta conc: " << delta_CompConc << " Value_old: " <<
									this->vec_CompConc_Water_elements[i][Current_Comp_Index_Water];
							}
							delta_CompConc = -this->vec_CompConc_Water_elements[i][Current_Comp_Index_Water];
						}
						if (fabs(delta_CompConc) > 0.0)
						{
							double SWAT = this->Data[i][variable_index_water_saturation];
							deltaMolesWater = delta_CompConc * SWAT * RPORV;
						}
						else
						{
							deltaMolesWater = 0.;
						}
					}; // END OF OGSProcessIndexCompWater != -99 condition
					if (OGSProcessIndexCompGas != -99)
					{
						double delta_CompConc = 0.0;
						weight = float(1.0 / vec_element_nodes.Size()); //arithmetic average SB redo wtp
						for (long k = 0; k < long(vec_element_nodes.Size()); k++)
						{
							// get the delta in gas components & add it up
							double delta_comp_value = this->NodeData[vec_element_nodes[k]->GetIndex()]->delta_CompConcInGas_Nodes[Current_Comp_Index_Gas];
							delta_CompConc += delta_comp_value * weight;
						}
						if (delta_CompConc < 0.0 && fabs(delta_CompConc) > this->vec_CompConc_Gas_elements[i][Current_Comp_Index_Gas])
						{
							if (verbosity > 1)
							{
								std::cout << " WARNING: The amount of component removed from gas is larger than the available amount! The value is corrected!"
									<< "\n";
								for (long k = 0; k < long(vec_element_nodes.Size()); k++)
								{
									std::cout << "     Node: " << vec_element_nodes[k]->GetIndex() << " Value: "
										<< this->NodeData[vec_element_nodes[k]->GetIndex()]->delta_CompConcInGas_Nodes[Current_Comp_Index_Gas] << "\n";
								}
								std::cout << "  Element: " << i << " delta conc: " << delta_CompConc << " Value_old: " <<
									this->vec_CompConc_Gas_elements[i][Current_Comp_Index_Gas];
							}
							delta_CompConc = -this->vec_CompConc_Gas_elements[i][Current_Comp_Index_Gas];
						}
						if (delta_CompConc > 0.0)
						{
							double SGAS = this->Data[i][variable_index_gas_saturation];
							deltaMolesGas = delta_CompConc * SGAS * RPORV;
						}
						else
						{
							deltaMolesGas = 0.;
						}
					};
					if (OGSProcessIndexCompOil != -99)
					{
						double delta_CompConc = 0.0;
						weight = float(1.0 / vec_element_nodes.Size()); //arithmetic average SB redo wtp
						for (long k = 0; k < long(vec_element_nodes.Size()); k++)
						{
							// get the delta in oil components & add it up
							double delta_comp_value = this->NodeData[vec_element_nodes[k]->GetIndex()]->delta_CompConcInOil_Nodes[Current_Comp_Index_Oil];
							delta_CompConc += delta_comp_value * weight;
						}
						if (delta_CompConc < 0.0 && fabs(delta_CompConc) > this->vec_CompConc_Oil_elements[i][Current_Comp_Index_Oil])
						{
							if (verbosity > 1)
							{
								std::cout << " WARNING: The amount of component removed from gas is larger than the available amount! The value is corrected!"
									<< "\n";
								for (long k = 0; k < long(vec_element_nodes.Size()); k++)
								{
									std::cout << "     Node: " << vec_element_nodes[k]->GetIndex() << " Value: "
										<< this->NodeData[vec_element_nodes[k]->GetIndex()]->delta_CompConcInOil_Nodes[Current_Comp_Index_Oil] << "\n";
								}
								std::cout << "  Element: " << i << " delta conc: " << delta_CompConc << " Value_old: " <<
									this->vec_CompConc_Oil_elements[i][Current_Comp_Index_Oil];
							}
							delta_CompConc = -this->vec_CompConc_Oil_elements[i][Current_Comp_Index_Oil];
						}
						if (delta_CompConc > 0.0)
						{
							double SOIL = this->Data[i][variable_index_oil_saturation];
							deltaMolesOil = delta_CompConc * SOIL * RPORV;
						}
						else
						{
							deltaMolesOil = 0.;
						}
					}
					// calculate new delta MLSC value based on the change in moles per phase
					double deltaMLSC = (deltaMolesWater + deltaMolesGas + deltaMolesOil) / (RPORV * 1000);
					double NewMLSC = this->Data[i][vec_variable_indicies_Comp_MolarDensity[current_run]] + deltaMLSC;
					//if (NewMLSC < 0.0)
					//	std::cout << " WARNING: New MLSC Value is negative!" << "\n";
					vec_dummy_mlsc.push_back(NewMLSC);
				} // END of if-condition: if(this->E100 != true)
			} // END of if-condition: if(vec_components_ECL_OGS[j][0] != -1)
		} // END of LOOP over ALL COMPONENTS

		// Update the MLSC vector if it is a E300 run
		if (this->E100 != true && this->dissolved_comps == true)
			vec_Comp_MLSC_elements.push_back(vec_dummy_mlsc);

		// Now interpolate the gas phase saturations and pressures from the nodes to the elements
		if (static_cast<int>(this->Phases.size()) > 1) // Always do it for the gas phase if we have more than one phase
		{
			double delta_Sgas = 0.0;
			weight = float(1.0 / vec_element_nodes.Size());        //arithmetic average
			// accumulate node contributions to element
			for (long j = 0; j < long(vec_element_nodes.Size()); j++)
			{
				double delta_value = NodeData[vec_element_nodes[j]->GetIndex()]->deltaSatGas;
				delta_Sgas += delta_value * weight;
			}
			if (delta_Sgas < 0.0) // decrease in Sgas
			{
				//check if delta_comp_valuesolved is not larger than the available volume of CO2 gas phase
				if (fabs(delta_Sgas) >(Data[i][variable_index_gas_saturation] - epsilon))
				{
					if (verbosity > 1)
					{
						std::cout << " WARNING: The volume of the gas removed is larger than the available volume of gas! The value is corrected!" << "\n";
						for (long j = 0; j < long(vec_element_nodes.Size()); j++)
						{
							std::cout << "     Node: " << vec_element_nodes[j]->GetIndex() << " delta_Sgas: "
								<< NodeData[vec_element_nodes[j]->GetIndex()]->deltaSatGas << "\n";
						};
						std::cout << "  Element: " << i << " deltaSatGas: " << delta_Sgas << " Sgas_old: " << Data[i][variable_index_gas_saturation] << "\n";
					}
					delta_Sgas = -Data[i][variable_index_gas_saturation]; // limit the delta
				}
			}
			// here set the delta to element
			if (delta_Sgas != 0.0)
			{
				if ((fabs(delta_Sgas) / Data[i][variable_index_gas_saturation]) < epsilon)
					vec_SGAS.push_back(Data[i][variable_index_gas_saturation]); // old saturation
				else
					vec_SGAS.push_back(Data[i][variable_index_gas_saturation] + delta_Sgas);
			}
			else
				vec_SGAS.push_back(Data[i][variable_index_gas_saturation]);

			if (vec_SGAS[vec_SGAS.size() - 1] < 0)
			{
				std::cout << " ERROR: The newly calculated gas saturation is negative!" << "\n";
				return 0;
			}

			// Now calculate interpolate the pressure change from the nodes to the elements
			double delta_pressure2 = 0.0;
			weight = float(1.0 / vec_element_nodes.Size());        //arithmetic average
			// accumulate node contributions to element
			for (long j = 0; j < long(vec_element_nodes.Size()); j++)
			{
				double delta_value = NodeData[vec_element_nodes[j]->GetIndex()]->deltaPress2;
				delta_pressure2 += delta_value * weight;
			}
			if (delta_pressure2 < 0.0)  // decrease in Pressure2
			{
				//check if delta_comp_valuesolved is not larger than the old Pressure of CO2 gas phase
				if (fabs(delta_pressure2) >(Data[i][variable_index_gas_pressure] - epsilon)) 
				{
					if (verbosity > 1)
					{
						std::cout << " WARNING: The pressure reduction in the gas phase is larger than the old gas pressure value! The value is corrected!" << "\n";
						for (long j = 0; j < long(vec_element_nodes.Size()); j++)
							std::cout << "     Node: " << vec_element_nodes[j]->GetIndex() << " Press: "
							<< NodeData[vec_element_nodes[j]->GetIndex()]->deltaPress2 << "\n";
						std::cout << "  Element: " << i << " deltaPress: " << delta_pressure2 << " Press_old: "
							<< Data[i][variable_index_gas_pressure];
					}
					delta_pressure2 = -Data[i][variable_index_gas_pressure]; // limit the delta
				}
			}
			// here set the delta ro element
			if (delta_pressure2 != 0.0)
			{
				if ((fabs(delta_pressure2) / Data[i][variable_index_gas_pressure]) < epsilon)
					vec_PRESS2.push_back(Data[i][variable_index_gas_pressure]); // old spressure
				else
					vec_PRESS2.push_back(Data[i][variable_index_gas_pressure] + delta_pressure2);
			}
			else
				vec_PRESS2.push_back(Data[i][variable_index_gas_pressure]);
			if (vec_PRESS2[vec_PRESS2.size() - 1] < 0)
			{
				std::cout << " ERROR: The newly calculated gas pressure is negative!" << "\n";
				return 0;
			}
			// Pa -> bar: /10000
			vec_PRESS2[i] /= 1.0e+5;
			
            if (M_Process == false)//&& m_pcs->M_feedback == false) // KB0815
                continue;
			else if (m_pcs->M_feedback != false) // KB0116 mechanical feedback allowed
				vec_PRESS2[i] += m_pcs->d_strain_2[i];

				// Now interpolate the saturation and pressure data for the oil phase
			// if we have three phases present
			if (static_cast<int>(this->Phases.size()) > 2 && this->E100 != true)
			{
				// Now interpolate the oil phase saturations from the nodes to the elements
				double delta_Soil = 0.0;
				weight = float(1.0 / vec_element_nodes.Size());        //arithmetic average
				// accumulate node contributions to element
				for (long j = 0; j < long(vec_element_nodes.Size()); j++)
				{
					double delta_value = NodeData[vec_element_nodes[j]->GetIndex()]->deltaSatOil;
					delta_Soil += delta_value * weight;
				}
				if (delta_Soil < 0.0) // decrease in soil
				{
					//check if delta is not larger than the available volume of oil phase
					if (fabs(delta_Soil) >(Data[i][variable_index_oil_saturation] - epsilon))
					{
						if (verbosity > 1)
						{
							std::cout << " WARNING: The volume of the oil phase removed is larger than the available volume of oil! The value is corrected!" << "\n";
							for (long j = 0; j < long(vec_element_nodes.Size()); j++)
							{
								std::cout << "     Node: " << vec_element_nodes[j]->GetIndex() << " Soil: "
									<< NodeData[vec_element_nodes[j]->GetIndex()]->deltaSatOil << "\n";
							};
							std::cout << "  Element: " << i << " deltaSatOil: " << delta_Soil << " Soil_old: " << Data[i][variable_index_oil_saturation] << "\n";
						}
						delta_Soil = -Data[i][variable_index_oil_saturation]; // limit the delta
					}
				}
				// here set the delta ro element
				if (delta_Soil != 0.0)
				{
					if ((fabs(delta_Soil) / Data[i][variable_index_oil_saturation]) < epsilon)
						vec_SOIL.push_back(Data[i][variable_index_oil_saturation]); // old satu
					else
						vec_SOIL.push_back(Data[i][variable_index_oil_saturation] + delta_Soil);
				}
				else
					vec_SOIL.push_back(Data[i][variable_index_oil_saturation]);

				if (vec_SOIL[vec_SOIL.size() - 1] < 0)
				{
					std::cout << " ERROR: The newly calculated oil saturation is negative!" << "\n";
					return 0;
				}
				// Now calculate interpolate the pressure change from the nodes to the elements
				double delta_pressure3 = 0.0;
				weight = float(1.0 / vec_element_nodes.Size());        //arithmetic average
				// accumulate node contributions to element
				for (long j = 0; j < long(vec_element_nodes.Size()); j++)
				{
					double delta_value = NodeData[vec_element_nodes[j]->GetIndex()]->deltaPress3;
					delta_pressure3 += delta_value * weight;
				}
				if (delta_pressure3 < 0.0)  // decrease in Pressure2
				{
					//check if delta_comp_valuesolved is not larger than the old Pressure of CO2 gas phase
					if (fabs(delta_pressure3) >(Data[i][variable_index_oil_pressure] - epsilon)) 
					{
						if (verbosity > 1)
						{
							std::cout << " WARNING: The pressure reduction in the oil phase is larger than the oil pressure value! The value is corrected!" << "\n";
							for (long j = 0; j < long(vec_element_nodes.Size()); j++)
								std::cout << "     Node: " << vec_element_nodes[j]->GetIndex() << " Press: "
								<< NodeData[vec_element_nodes[j]->GetIndex()]->deltaPress3 << "\n";
							std::cout << "  Element: " << i << " deltaPress: " << delta_pressure3 << " Press_old: "
								<< Data[i][variable_index_oil_pressure];
						}
						delta_pressure3 = -Data[i][variable_index_oil_pressure]; // limit the delta
					}
				}
				// here set the delta ro element
				if (delta_pressure3 != 0.0)
				{
					if ((fabs(delta_pressure3) / Data[i][variable_index_oil_pressure]) < epsilon)
						vec_PRESS3.push_back(Data[i][variable_index_oil_pressure]); // old pressure
					else
						vec_PRESS3.push_back(Data[i][variable_index_oil_pressure] + delta_pressure3);
				}
				else
					vec_PRESS3.push_back(Data[i][variable_index_oil_pressure]);
				if (vec_PRESS3[vec_PRESS3.size() - 1] < 0)
				{
					std::cout << " ERROR: The newly calculated oil pressure is negative!" << "\n";
					return 0;
				}
				// Pa -> bar: /10000
				vec_PRESS3[i] /= 1.0e+5;
			} // end of phases > 2 condition
		} // end of phases > 1 condition
		else if (static_cast<int>(this->Phases.size()) == 1)   // KB0714: One Phase only
		{
			int variable_index_water_pressure = this->GetVariableIndex("PRESSURE");
			// Now calculate interpolate the pressure change from the nodes to the elements
			double delta_pressure1 = 0.0;
			weight = float(1.0 / vec_element_nodes.Size());        //arithmetic average
			// accumulate node contributions to element
			for (long j = 0; j < long(vec_element_nodes.Size()); j++)
			{
				double delta_value = NodeData[vec_element_nodes[j]->GetIndex()]->deltaPress1;
				delta_pressure1 += delta_value * weight;
			}
			if (delta_pressure1 < 0.0)  // decrease in Pressure1
			{
				//check if delta_comp_valuesolved is not larger than the old Pressure of CO2 gas phase
				if (fabs(delta_pressure1) >(Data[i][variable_index_water_pressure] - epsilon))
				{
					if (verbosity > 1)
					{
						std::cout << " WARNING: The pressure reduction in the gas phase is larger than the old gas pressure value! The value is corrected!" << "\n";
						for (long j = 0; j < long(vec_element_nodes.Size()); j++)
							std::cout << "     Node: " << vec_element_nodes[j]->GetIndex() << " Press: "
							<< NodeData[vec_element_nodes[j]->GetIndex()]->deltaPress1 << "\n";
						std::cout << "  Element: " << i << " deltaPress: " << delta_pressure1 << " Press_old: "
							<< Data[i][variable_index_water_pressure];
					}
					delta_pressure1 = -Data[i][variable_index_water_pressure]; // limit the delta
				}
			}
			// here set the delta ro element
			if (delta_pressure1 != 0.0)
			{
				if ((fabs(delta_pressure1) / Data[i][variable_index_water_pressure]) < epsilon)
					vec_PRESS1.push_back(Data[i][variable_index_water_pressure]); // old spressure
				else
					vec_PRESS1.push_back(Data[i][variable_index_water_pressure] + delta_pressure1);
			}
			else
				vec_PRESS1.push_back(Data[i][variable_index_water_pressure]);
			if (vec_PRESS1[vec_PRESS1.size() - 1] < 0)
			{
				std::cout << " ERROR: The newly calculated gas pressure is negative!" << "\n";
				return 0;
			}
			// Pa -> bar: /10000
			vec_PRESS1[i] /= 1.0e+5;

			//Testoutput KB
			//cout << vec_PRESS1[i] << "\n";

			if (M_Process == false)// KB0815
				continue;
			else if (m_pcs->M_feedback)
				//if (m_pcs->d_strain_2[i] > 0) vec_PRESS1[i] -= m_pcs->d_strain_2[i]; // increase in strain --> decrease in pressure
				vec_PRESS1[i] += m_pcs->d_strain_2[i]; // increase in strain --> decrease in pressure
				//else vec_PRESS1[i] += m_pcs->d_strain_2[i]; // decrease in strain --> increase in pressure 

		}
	} // end of element loop
	return test;
}


/*-------------------------------------------------------------------------
   GeoSys - Function: WriteDataBackToEclipse
   Task: Gather Data from OGS, write it back to Eclipse (.data, include files)
   Return: int
   Programming: xx/xx
   Modification: 04/2013 WTP added temperature coupling
   11/2013 WTP added multiphase-multicomponent support
   03/2014 WTP reduced function length by consolidation & outsourcing of functions
   04/2014 WTP added support for 3phases
   07/2014 KB added support for 1 phase and deformation feedback
   -------------------------------------------------------------------------*/
int CECLIPSEData::WriteDataBackToEclipse(CReadTextfiles_ECL* eclFFile, CReadTextfiles_ECL* eclDataFile, CRFProcess* m_pcs, std::string folder)
{
	std::string Filename;
	//MeshLib::CElem* m_element = NULL;
	//CFEMesh* m_msh = fem_msh_vector[0];
	clock_t start, finish;
	double time;
	vector <std::string> vecString;
	std::string Keyword;
	ostringstream temp;
	std::string tempstring;
	//int idx = -1;

	typeExponentialNumber tempNumber;
	//WTP const double epsilon = 1e-7;
	int j_max;

	//get no of nodes in mesh
	//long nnodes = fem_msh_vector[0]->GetNodesNumber(false);
	//CB PoroPermData exchange
	bool poroflag = false;
	bool permxxflag = false;
	//bool permyyflag = false;
	//bool permzzflag = false;
	vector <double> poroperm;
	vector <string> vecString2;
	vector <string> vecString3;
	std::string tempstring2;
	std::string tempstring3;
	double *tensor = NULL;
	long count;
	//CB_merge_0513 
	//int idx_kxx, idx_kyy, idx_kzz; //WTP
	int idx_n = 0;
	std::ostringstream sstream;
	sstream.precision(8);

	//WTP variables for temperautre coupling to ecl
	//CRFProcess* t_pcs = NULL;
	//int idx_T = -1;	// WTP Index of temperature data in
	//double temp_temp = -99.;
	//double ele_temp = -99.;

	if (verbosity > 2)
		std::cout << "      WriteDataBackToEclipse" << "\n";

	start = clock();

	//WTP 04/2014: FUNCTION CALL FOR DELTA CALCULATION
	CalculateDeltaGeoSysECL(m_pcs);

	//WTP 04/2014: FUNCTION CALL FOR THE INTERPOLATION OF THE DELTA VALUES
    InterpolateDeltaGeoSysECL(m_pcs);

    // NOW WRITE THE DATA INTO THE FILES
	//write dissolved components (RS/MLSC), SGAS, PRESSURE, POROPERM  and temperature data into the Eclipse restart file
    if (this->E100 == true && this->Phases.size() > 1)
	{
		//write RS into the Eclipse restart file
		Keyword = " 'RS";
		vecString.clear(); // SB redo wtp

		//check consitency of data
		if (static_cast<int>(vec_RS_elements.size()) != this->elements)
		{
			std::cout <<
				" ERROR: Number of RS data entries does not match element count!"
				<< "\n";
			//system("Pause"); // SB redo wtp
			return 0;
		}
		for (long i = 0; i < this->elements; i = i + 4)
		{
			tempstring = "";
			j_max = 4;
			if ((this->elements - i) < 4)
				j_max = (this->elements - i);
			for (int j = 0; j < j_max; j++)
			{
				tempstring += "   ";
				tempNumber = this->RoundEXP(vec_RS_elements[i + j], 8);    // HC FOR 1 COMPONENT ONLY!
				tempstring += this->AddZero(tempNumber.Number, 8, false);
				tempstring += "E";
				if (tempNumber.Exponent >= 0)
				{
					tempstring += "+";
					tempstring += this->AddZero(tempNumber.Exponent, 2, true);
				}
				else
					tempstring += this->AddZero(tempNumber.Exponent, 3, true);
			}
			vecString.push_back(tempstring);
		}
		if (ReplaceASectionInData(eclFFile, Keyword, vecString, false) == false)
		{
			std::cout << " ERROR: Replacing a section in the virtual file for Keyword " << Keyword <<
				" didn't work!" << "\n";
			//system("Pause"); // SB redo wtp
			return 0;
		}
	}
    else if (this->E100 == false) // this is the E300 run
	{
		// local test value
		int current_run = -1;
		//Filename = folder + "TemporaryResults.F" + AddZero(m_pcs->Tim->step_current - 1,4,true);
		for (unsigned int n = 0; n < this->vec_components_ECL_OGS.size(); n++)
		{

			if (vec_components_ECL_OGS[n][0] != -1)
			{
				current_run++;
				std::stringstream ss;
				ss << vec_components_ECL_OGS[n][0];
				std::string dummy_str = "MLSC" + ss.str();
				Keyword = " '" + dummy_str;
				ss.clear();

				vecString.clear();
				//check consitency of data
				if (static_cast<int>(vec_Comp_MLSC_elements.size()) != this->elements)
					//if (int(vec_dummy_mlsc[current_run].size()) != this->elements)
				{
					std::cout <<
						" ERROR: Number of MLSC data entries does not match element count!"
						<< "\n";
					system("Pause");
					return 0;
				}
				for (long i = 0; i < this->elements; i = i + 4)
				{
					tempstring = "";
					j_max = 4;
					if ((this->elements - i) < 4)
						j_max = (this->elements - i);
					for (int j = 0; j < j_max; j++)
					{
						tempstring += "   ";
						tempNumber = this->RoundEXP(vec_Comp_MLSC_elements[i + j][current_run], 8);
						tempstring += this->AddZero(tempNumber.Number, 8, false);
						tempstring += "E";
						if (tempNumber.Exponent >= 0)
						{
							tempstring += "+";
							tempstring += this->AddZero(tempNumber.Exponent, 2, true);
						}
						else
							tempstring += this->AddZero(tempNumber.Exponent, 3, true);
					}
					vecString.push_back(tempstring);
				}

				if (ReplaceASectionInData(eclFFile, Keyword, vecString, false) == false)
				{
					std::cout << " ERROR: Replacing a section in the virtual file for Keyword " << Keyword <<
						" didn't work!" << "\n";
					//system("Pause"); // SB redo wtp
					return 0;
				}
			}
		}
	}

	// SB redo wtp
	// CB SGAS independent of e100 / e300
	Keyword = " 'SGAS";
	vecString.clear();
	//check consitency of data
	if (this->Phases.size() > 1)
	{
		if (static_cast<int>(vec_SGAS.size()) != this->elements)
		{
			std::cout << " ERROR: Number of SGAS data entries does not match element count!" << "\n";
			return 0;
		}
		for (long i = 0; i < this->elements; i = i + 4)
		{
			tempstring = "";
			j_max = 4;
			if ((this->elements - i) < 4)
				j_max = (this->elements - i);
			for (int j = 0; j < j_max; j++)
			{
				tempstring += "   "; // 1lz
				tempNumber = this->RoundEXP(vec_SGAS[i + j], 8);
				tempstring += this->AddZero(tempNumber.Number, 8, false); // not for poroperm
				tempstring += "E";
				if (tempNumber.Exponent >= 0)
				{
					tempstring += "+";
					tempstring += this->AddZero(tempNumber.Exponent, 2, true);
				}
				else
					tempstring += this->AddZero(tempNumber.Exponent, 3, true);
			}
			vecString.push_back(tempstring);
		}
		if (ReplaceASectionInData(eclFFile, Keyword, vecString, false) == false)
		{
			std::cout << " ERROR: Replacing a section in the virtual file for Keyword " << Keyword <<
				" didn't work!" << "\n";
			return 0;
		}
	}

	// CB PRESS independent of e100 / e300
	Keyword = " 'PRESSURE";
	vecString.clear();
	if (m_pcs->therzagi == 1)// KB:Manipulating the fluid pressure for special BC in Therzagi Benchmark
	{
		if (vec_PRESS2.size() == 0)
		{
			//for (long i = 0; i < vec_PRESS1.size(); i++)
			//{
			for (long l = 0; l < this->eclgrid.size(); l++)
			{
				if (this->eclgrid[l]->row == 1)

					vec_PRESS1[l] = 1.00000000E+00;
			}
			//}
		}
	}
	if (this->Phases.size() > 1)
	{
		//check consitency of data
		if (static_cast<int>(vec_PRESS2.size()) != this->elements)
		{
			std::cout << " ERROR: Number of PRESSURE2 data entries does not match element count!" << "\n";
			return 0;
		}
		for (long i = 0; i < this->elements; i = i + 4)
		{
			tempstring = "";
			j_max = 4;
			if ((this->elements - i) < 4)
				j_max = (this->elements - i);
			for (int j = 0; j < j_max; j++)
			{
				tempstring += "   "; // 1lz
				tempNumber = this->RoundEXP(vec_PRESS2[i + j], 8);
				tempstring += this->AddZero(tempNumber.Number, 8, false); // not for poroperm
				tempstring += "E";
				if (tempNumber.Exponent >= 0)
				{
					tempstring += "+";
					tempstring += this->AddZero(tempNumber.Exponent, 2, true);
				}
				else
					tempstring += this->AddZero(tempNumber.Exponent, 3, true);
			}
			vecString.push_back(tempstring);
		}
		if (ReplaceASectionInData(eclFFile, Keyword, vecString, false) == false)
		{
			std::cout << " ERROR: Replacing a section in the virtual file for Keyword " << Keyword <<
				" didn't work!" << "\n";
			return 0;
		}
	}
	if (this->Phases.size() == 1)
	{
		//check consitency of data
		if (static_cast<int>(vec_PRESS1.size()) != this->elements)
		{
			std::cout << " ERROR: Number of PRESSURE1 data entries does not match element count!" << "\n";
			return 0;
		}
		for (long i = 0; i < this->elements; i = i + 4)
		{
			tempstring = "";
			j_max = 4;
			if ((this->elements - i) < 4)
				j_max = (this->elements - i);
			for (int j = 0; j < j_max; j++)
			{
				tempstring += "   "; // 1lz
				tempNumber = this->RoundEXP(vec_PRESS1[i + j], 8);
				tempstring += this->AddZero(tempNumber.Number, 8, false); // not for poroperm
				tempstring += "E";
				if (tempNumber.Exponent >= 0)
				{
					tempstring += "+";
					tempstring += this->AddZero(tempNumber.Exponent, 2, true);
				}
				else
					tempstring += this->AddZero(tempNumber.Exponent, 3, true);
			}
			vecString.push_back(tempstring);
		}
		if (ReplaceASectionInData(eclFFile, Keyword, vecString, false) == false)
		{
			std::cout << " ERROR: Replacing a section in the virtual file for Keyword " << Keyword <<
				" didn't work!" << "\n";
			return 0;
		}
	}

	// CB Poropermdata independent of e100 / e300
	//for (long i = 0; i < mmp_vector.size(); i++) {
	for (std::size_t i = 0; i < mmp_vector.size(); i++) // WTP
	{
		if (mmp_vector[i]->porosity_model == 12)
		{
			poroflag = true; // poro update from geochemistry  
			idx_n = m_pcs->GetElementValueIndex("POROSITY") + 1;
			if (idx_n < 0)
				std::cout << " WARNING: No POROSITY ele index found with PCS " << m_pcs->getProcessType() << "." << "\n";
		}
		if (mmp_vector[i]->permeability_porosity_model == 8)    // perm update from porochange
			permxxflag = true;
	}
	if (permxxflag)        //write perm into the Eclipse restart file
	{
		j_max = 4;
		count = 0;
		vecString.clear();
		vecString2.clear();
		vecString3.clear();
		for (long i = 0; i < long(eclipse_ele_active_flag.size()); i = i + j_max)
		{
			tempstring = "";
			tempstring2 = "";
			tempstring3 = "";
			if (long(eclipse_ele_active_flag.size()) - i < j_max)
				j_max = (long(eclipse_ele_active_flag.size()) - i);
			for (int j = 0; j < j_max; j++)
			{
				tempstring += " ";  // 1lz
				tempstring2 += " "; // 1lz
				tempstring3 += " "; // 1lz
				if (eclipse_ele_active_flag[i])
				{
					tensor = mmp_vector[0]->PermeabilityTensor(count);
					//x
					sstream << fixed << scientific << tensor[0] / 9.869233e-16; // m²->mD
					tempstring += sstream.str();
					sstream.str(""); sstream.clear();
					//y
					sstream << fixed << scientific << tensor[4] / 9.869233e-16;
					tempstring2 += sstream.str();
					sstream.str(""); sstream.clear();
					//z
					sstream << fixed << scientific << tensor[8] / 9.869233e-16;
					tempstring3 += sstream.str();
					sstream.str(""); sstream.clear();
					count++;
				}
				else
				{
					tempstring += "0.123";
					tempstring2 += "0.123";
					tempstring3 += "0.123";
				}
			}
			vecString.push_back(tempstring);
			vecString2.push_back(tempstring2);
			vecString3.push_back(tempstring3);
		}

		if (this->Radialmodell)
			Keyword = "PERMR";
		else
			Keyword = "PERMX";

		if (PoroPermIncludeFile)
			WriteIncludeFile(pathECLFolder + "GRID_PROPS.INC", Keyword, vecString, false);
		else if (ReplaceASectionInData(eclDataFile, Keyword, vecString, false) == false)
		{
			std::cout << " ERROR: Replacing a section in the virtual file for Keyword " << Keyword <<
				" didn't work!" << "\n";
			return 0;
		}

		if (this->Radialmodell)
			Keyword = "PERMTHT";
		else
			Keyword = "PERMY";

		if (PoroPermIncludeFile)
			WriteIncludeFile(pathECLFolder + "GRID_PROPS.INC", Keyword, vecString2, true);
		else if (ReplaceASectionInData(eclDataFile, Keyword, vecString, false) == false)
		{
			std::cout << " ERROR: Replacing a section in the virtual file for Keyword " << Keyword <<
				" didn't work!" << "\n";
			return 0;
		}
		Keyword = "PERMZ";
		if (PoroPermIncludeFile)
			WriteIncludeFile(pathECLFolder + "GRID_PROPS.INC", Keyword, vecString3, true);
		else if (ReplaceASectionInData(eclDataFile, Keyword, vecString, false) == false)
		{
			std::cout << " ERROR: Replacing a section in the virtual file for Keyword " << Keyword <<
				" didn't work!" << "\n";
			return 0;
		}
	}
	if (poroflag)
	{        //write poro into the Eclipse restart file
		j_max = 4;
		count = 0;
		vecString.clear();
		for (long i = 0; i < long(eclipse_ele_active_flag.size()); i = i + j_max)
		{
			tempstring = "";
			if (long(eclipse_ele_active_flag.size()) - i < j_max)
				j_max = (long(eclipse_ele_active_flag.size()) - i);
			for (int j = 0; j < j_max; j++)
			{
				tempstring += " "; // 1lz
				if (eclipse_ele_active_flag[i])
				{
					sstream << fixed << scientific << m_pcs->GetElementValue(count, idx_n);
					tempstring += sstream.str();
					count++;
					sstream.str(""); sstream.clear();
				}
				else
					tempstring += "0.123";
			}
			vecString.push_back(tempstring);
		}

		Keyword = "PORO";

		if (PoroPermIncludeFile)
			WriteIncludeFile(pathECLFolder + "GRID_PROPS.INC", Keyword, vecString, true);
		else if (ReplaceASectionInData(eclDataFile, Keyword, vecString, false) == false)
		{
			std::cout << " ERROR: Replacing a section in the virtual file for Keyword " << Keyword <<
				" didn't work!" << "\n";
			return 0;
		}
	}
    if (T_Process == true)
        GetHeatDataFromOGS(1);
    	

	// clearing the vectors
	vecString.clear();
	vecString2.clear();
	vecString3.clear();
	vec_Comp_MLSC_elements.clear();
	vec_RS_elements.clear();
	vec_SWAT.clear();
	vec_SGAS.clear();
	vec_SOIL.clear();
	vec_PRESS1.clear(); //KB
	vec_PRESS2.clear();
	vec_PRESS3.clear();

	finish = clock();

	if (verbosity > 2)
		std::cout << "        done.";
	time = (static_cast<double>(finish)-static_cast<double>(start)) / CLOCKS_PER_SEC;
	if (verbosity > 2)
		std::cout << "                    Time: " << time << " seconds." << "\n";

	return 1;
}


/*-------------------------------------------------------------------------
GeoSys - Function: GetHeatDataFromOGS
Task: Function to get heat data fomr OGS and hand it over to ECL via an include file
Return: nothing
Programming: 01/2015 WTP
Modification:
-------------------------------------------------------------------------*/
    void CECLIPSEData::GetHeatDataFromOGS(int idx_adjust)
{
    // setting bool flag for save function (not implemented yet)
    this->TempIncludeFile = true;
    CRFProcess* t_pcs = NULL;
    t_pcs = PCSGet("HEAT_TRANSPORT");
    int idx_T = t_pcs->GetNodeValueIndex("TEMPERATURE1") + idx_adjust;

    //MeshLib::CElem* m_element = NULL;
    CFEMesh* m_msh = fem_msh_vector[0];

    //clock_t start, finish;
    //double time;

    std::string Filename;
    vector <std::string> vecString;
    std::string Keyword = "TEMPI";
    std::string tempstring;
    int j_max = 4;
    std::ostringstream sstream;
    sstream.precision(8);
    int index_adjust = 0;   //BW, correct the index of active cells
    // running a loop over all, active and non-active elements of Eclipse
    for (long i = 0; i < long(eclipse_ele_active_flag.size()); i = i + j_max)
    {
        tempstring = "";
        if (long(eclipse_ele_active_flag.size()) - i < j_max)
            j_max = (long(eclipse_ele_active_flag.size()) - i);

        for (int j = 0; j < j_max; j++)
        {
            tempstring += " ";
            // if the current element is active e.g. also used in Geosys get the temperature
            if (eclipse_ele_active_flag[i + j])
            {
                // here get the index of every node in the element 
                Math_Group::vec <MeshLib::CNode*> vec_element_nodes(this->eclgrid[i + j - index_adjust]->MshNodeindex.size()); //BW:might not less than 8 nodes data need to write back
                for (int tempi = 0; tempi < static_cast<int>(this->eclgrid[i + j - index_adjust]->MshNodeindex.size()); tempi++)
                    vec_element_nodes[tempi] = m_msh->nod_vector[this->eclgrid[i + j - index_adjust]->MshNodeindex[tempi]];
                double ele_temp = 0.0;
                double weight = float(1.0 / vec_element_nodes.Size()); // getting the correct weight for each nodal value
                // run a loop over the number of nodes of the element
                for (long k = 0; k < long(vec_element_nodes.Size()); k++)
                {
                    // save the node index of the nodes
                    long loc_nd_idx = vec_element_nodes[k]->GetIndex();
                    // read the corresponging temperature
                    double temp_temp = t_pcs->GetNodeValue(loc_nd_idx, idx_T);
                    // add it up 
                    ele_temp += temp_temp * weight;
                }
                // Since the temperature units used in Geosys and ECL differ they have to be converted
                sstream << fixed << scientific << ele_temp - 273.15;
                tempstring += sstream.str();
                sstream.str(""); sstream.clear();
            }
            else
            {
                // all inactive cells get some random number
                tempstring += "999.9";
                // BW: inactive cell count to adjust the index
                index_adjust ++;
            }
        }
        vecString.push_back(tempstring);
    }
    // write the data into the corresponding include file
	if (UsePrecalculatedFiles == false)
		WriteIncludeFile(pathECLFolder + "TEMP.INC", Keyword, vecString, false);
	else
	{
		ostringstream temp; //WTP pretty time consuming constructor...
		std::string file_id = "_SAVE";

		int timestep = aktueller_zeitschritt + timestep_adjust_initial + this->timestep_adjust_iteration_tot;
		temp.str(""); temp.clear(); temp << timestep;
		if (timestep < 10)
			file_id += "000" + temp.str();
		else if (timestep < 100)
			file_id += "00" + temp.str();
		else if (timestep < 1000)
			file_id += "0" + temp.str();
		else
			file_id += temp.str();

		WriteIncludeFile(pathECLFolder + "TEMP.INC" + file_id, Keyword, vecString, false);
	}
}
/*-------------------------------------------------------------------------
   GeoSys - Function: InitializeProject
   Task: Function to initialize the OGS ECL run. Consists of functions which are only
   executed in the first time step.
   Return: nothing
   Programming: 02/2014 WTP
   Modification:
   -------------------------------------------------------------------------*/
void CECLIPSEData::InitializeProject(CRFProcess* m_pcs, CReadTextfiles_ECL* eclDataFile)
{
	// Create a int vector field of the process names that are going to be exchanged between ECL and OGS
	// format being: Component number in ECL - flag for transported dissolved in water - flag for transportet in gas
	// with 1=true and 0=false
	// no_comps_ECL = m_pcs->vec_component_pcs_names.size(); // get the total number of components
	std::vector <int> vec_dummy_int;
	typedef std::vector<std::vector<std::string> >::const_iterator CI;        // defintion of a constant iterator to loop over a vector<vector<string>> field
	for (CI pcs_name = m_pcs->vec_component_pcs_names.begin();                // instead of using normal counting indicies the defined iterator can be used
		pcs_name != m_pcs->vec_component_pcs_names.end(); ++pcs_name)    // resulting in a more direct data access
	{
		vec_dummy_int.push_back(-1);
		for (int j = 1; j < 4; j++)
		{
			string tempstring = (*pcs_name)[j];
			if (tempstring == "void" || tempstring == "0")
				vec_dummy_int.push_back(0);
			else
			{
				vec_dummy_int.push_back(1);
				if (j == 1)
					dissolved_comps = true;        // bool flag for queries
			}
		}
		vec_components_ECL_OGS.push_back(vec_dummy_int);
		vec_dummy_int.clear();
	}

	// Copy flags from pcs input to eclipse data structure
	this->UsePrecalculatedFiles = m_pcs->PrecalculatedFiles;
	this->UseSaveEclipseDataFiles = m_pcs->SaveEclipseDataFiles;
	this->reservoir_conditions = m_pcs->ECLunits_reservoir_conditions;

	// Analyse Eclipse Data-File
	this->AnalyzeDataFromInputFile(eclDataFile, m_pcs);

	if (this->SurfaceDensity_Gas_E100 <= 0 && this->E100 == true)
	{
		std::cout <<
			" ERROR: The gas density at surface conditions was not read properly (rho <= 0!)." << "\n";
		exit(0);
	}

	// Read the external .well-File if needed
	if (this->existWells == true)
		this->ReadWellData(pathECLWellFile);

	// Update path to eclipse executable if needed
	if (UseEclrun)
	{
		if (E100)
			m_pcs->simulator_path += " eclipse";
		else
			m_pcs->simulator_path += " e300";
	}

	//WTP: Consistency check
	// if multi phase flow is selected in ogs  but only one phase velocity is available in eclipse display a warning
	if (this->Phases.size() == 1 &&
		m_pcs->getProcessType() == FiniteElement::MULTI_PHASE_FLOW)
	{
		std::cout << "  ERROR: Only one Phase is definded in ECL but OGS is set to MULTI_PHASE_FLOW!";
		std::cout << "           -> Not all phase velocities needed by OGS will be available!" << "\n";
		exit(0);
	}
}

/*-------------------------------------------------------------------------
   GeoSys - Function: GeneralBookkeeping
   Task: Function to execute several functions in the first time step
   Return: nothing
   Programming: 02/2014 WTP
   Modification: Add information about collapsed cell 04/2014 WB

   -------------------------------------------------------------------------*/
void CECLIPSEData::GeneralBookkeeping(void)
{
	// Read ECLIPSE model grid, check consistency with GeoSys mesh and construct faces
	this->ReadEclipseGrid(pathECLProject + ".FGRID");
	// Read a List containing nodeindex and corresponding OGS Elements
	this->ReadCorrespondingList(pathECLProject + ".list");

	//Check if Eclipse grid is identical with Geosys grid
	if (CompareElementsGeosysEclipse() == 0)
	{
		std::cout <<
			" ERROR: The Eclipse grid is not identical to the Geosys grid!"
			<< "\n";
        std::cout << flush;
		//system("Pause");
		exit(0);
	}

    //Read boundary conditions
	//this->ReadPositionBoundaryCondition(pathECLProject + ".DATA");
	//Create Faces
	//std::cout << " creating faces " << "\n";
	this->CreateFaces(pathECLProject);
	//Determine neighbouring elements (only once)
	//std::cout << " determine Neighbours " << "\n";
	this->DetermineNeighbourElements(pathECLProject);
	// Connect Faces to each block
	//std::cout << " connect faces to elements " << "\n";
	this->ConnectFacesToElements();
	// check for radial model and one column
	int iindex = 0;
	if (this->Radialmodell == true)
	{
		if (this->rows > 1) // check if more than one column is active in radial model
		{
			int count = 0;
			for (long ii = 0; ii < long(this->eclgrid.size()); ii++)
			{
				int found = 0;
				if ((this->eclgrid[ii]->layer == 1) && (this->eclgrid[ii]->column == 1))
				{
					found++;
					if (this->eclgrid[ii]->active == 1)
					{
						count++;
						iindex = ii;
					}
				}
				if (found == this->rows)
					break;  // foun all cells in layer ==1 and column == 1
			}
			if (count > 1)
				std::cout << "\n" << "\n" <<
				" ERROR:  Definition of radial flow model has more than one active column."
				<< "\n";
		}
		// test if faces are perpenducular to I coordinate axis
		// WTP use static casts
		//for(int ii = 0; ii < int(this->eclgrid[iindex]->connected_faces.size());
		for (int ii = 0; ii < static_cast<int>(this->eclgrid[iindex]->connected_faces.size()); ii++)
		{
			if (this->faces[ii]->model_axis.find("I") == 0)
			{
				double* nvec = this->faces[ii]->PlaneEquation->GetNormalVector();
				if (nvec[0] == 1.0)
					this->Radial_I = true;
			}
		}
		if (Radial_I != true)
		{
			std::cout << "\n" <<
				" ERROR: I-face is not perpendicular to coordinate axis " << "\n";
		}
	}
}

/*-------------------------------------------------------------------------
   GeoSys - Function: RunEclipse
   Task: Preprocessing, running Eclipse, Postprocessing
   Return: nothing
   Programming: 09/2009 BG / SB
   Modification: 04/2013 WTP: added routine to save .data files from Ecl
   03/2014 WTP: tried to improve code readability, changed data/file handling
   07/2015 WTP: included verbosity flags
   -------------------------------------------------------------------------*/
int CECLIPSEData::RunEclipse(long Timestep, CRFProcess* m_pcs)
{
	const clock_t start = clock();
	std::cout << "\n Starting OGS-ECLIPSE coupling routine... " << "\n";

	// get the actual timestep
	if (m_pcs->Tim->step_current == 1)
	{
		actual_time = 0;
	}
	else
	{
		actual_time += m_pcs->Tim->time_step_vector[m_pcs->Tim->step_current - 2];
	}

	if (m_pcs->Iterative_Eclipse_coupling == true)
		this->timestep_adjust_iteration_tot = this->timestep_adjust_iteration_tot + 1;
	// set the filenames & paths
	if (this->SetFilenamesAndPaths(m_pcs, Timestep) == false)
	{
		std::cout << " ERROR: Filenames or paths are not correct! " << "\n";
		exit(0);
	}
	// Get geomechanics process
	if (M_Process)
		for (size_t i = 0; i < pcs_vector.size(); i++)
			if (isDeformationProcess(pcs_vector[i]->getProcessType()))
			{
				dm_pcs = (CRFProcessDeformation*)pcs_vector[i];
				break;
			};

	// declare reading function for the *.F***  and the *.Data textfile first
	CReadTextfiles_ECL* eclFFile;
	eclFFile = new CReadTextfiles_ECL;
	CReadTextfiles_ECL* eclDataFile;
	eclDataFile = new CReadTextfiles_ECL;

	// Now read the .Data file --> Maybe add some conditions later?
	if (eclDataFile->Read_Text(this->pathECLProject + ".DATA") == true)
	{
		std::cout << " ERROR: Could not read the textfile: " << this->pathECLProject + ".DATA" << "! \n";
		exit(0);
	}

	// Read the *.F** file --> further conditions should be added later
	if (Timestep == 1 && m_pcs->iter_outer_cpl == 0)
    {
        InitializeProject(m_pcs, eclDataFile);
    }
    if (Timestep > 1 || (m_pcs->iter_outer_cpl > 0))
    {
		// WTP DEBUG std::cout << " PATH ECL FILE: " << pathECLFFile << "\n";
        if (eclFFile->Read_Text(pathECLFFile) == true)
        {
            std::cout << " ERROR: Could not read the textfile: " << pathECLFFile << "! \n";
            std::cout << flush;
            exit(0);
        }
        //write amount of dissolved gas and later maybe also water density back to Eclipse; ToDo: Write new Water density back to Eclipse
        //if ((this->dissolved_comps == true) || (T_Process == true) || dm_pcs)
		if ((this->dissolved_comps == true) || (T_Process == true) || (M_Process == true))
		{
			if (this->UsePrecalculatedFiles == true)
			{
					std::cout << " ATTENTION: Running simulation using precalculated *.F*** files! " << "\n";
					std::cout << "            No feedback on ECLIPSE simulation possible!"
					<< "\n";
			}
            //else
			if (this->WriteDataBackToEclipse(eclFFile, eclDataFile, m_pcs, pathECLFolder) == 0)
            {
                std::cout <<
					" ERROR: WriteDataBackToEclipse() was not finished properly!"
                    << "\n";
                std::cout << flush;
                exit(0);
            }
        }
	}
	else
	{
		if (T_Process == true && timestep_adjust_initial > 0) // TEMPI only works with restarts anyway
        {
            this->ReadEclipseGrid(pathECLProject + ".FGRID");
            this->ReadCorrespondingList(pathECLProject + ".list");
            GetHeatDataFromOGS(0);
        }
	}

	//Execute Eclipse, try several times in case of problems finding the license
	if (verbosity > 2)
		std::cout << "        Calling Eclipse " << "\n";
	ExecuteEclipse(eclDataFile, eclFFile, Timestep, m_pcs);

	if (Timestep == 1 && m_pcs->iter_outer_cpl == 0)
	{
		GeneralBookkeeping();
	}

	
	//Read the ECLIPSE model output data
	if ( m_pcs->Iterative_Eclipse_coupling == true ) //KB0116
		this->ReadEclipseData(pathECLProject + ".F" + AddZero(this->timestep_adjust_iteration_tot + timestep_adjust_initial, 4, true));
	else
		this->ReadEclipseData(pathECLProject + ".F" + AddZero(Timestep + timestep_adjust_initial, 4, true));
	
	// WTP if any component is transported, convert units to 
	// DEBUG: Consistency Check
	if (this->dissolved_comps == true)
	{
        this->ConvertEclipseDataToUniformUnits(m_pcs, Timestep);
	}
	
	this->MakeNodeVector();
	
	// WTP changed to c++ style casts
	//for (int i = 0; i < int(static_cast<int>(this->Phases.size()); i++)
	for (int i = 0; i < static_cast<int>(this->Phases.size()); i++)
	{
		//Get the flow for each face
		this->GetFlowForFaces(i);
		//this->CalcBlockBudget(i);
		//this->GetVelForFaces();
		//interpolate face flows to nodes
		this->InterpolateDataFromBlocksToNodes(m_pcs, pathECLFolder, i, Timestep);
	}
	
	WriteDataToGeoSys(m_pcs, pathECLFolder);
	
	// WTP Save old .data files from eclipse for benchmarking if wanted 
	// we should also add saving of the *.FXXX files since OGS writes stuff back to ECLIPSE using them
	// --> might want to add saving the temperature files
	/*if (UseSaveEclipseDataFiles == true)
	{
		SaveEclipseDataFile(Timestep, m_pcs);
	}*/

	if (UsePrecalculatedFiles == false || UseSaveEclipseDataFiles == true)
	{
			CleanUpEclipseFiles(pathECLFolder, pathECLProject, Timestep, m_pcs);
	}

	const clock_t finish = clock();
	//WTP changed to c++ casts
	const double time = (static_cast<double>(finish)-static_cast<double>(start)) / CLOCKS_PER_SEC;
	std::cout << " -> Total time for Eclipse coupling routine: " << time << " seconds." << "\n";
	std::cout << "\n";
	std::cout << flush;
	return 1;
}

//TODO:

// interpolation of pressure to nodes-> volume is just devided by the number of nodes -> not exact if shape is not a quader

/*-------------------------------------------------------------------------
GeoSys - Function: WriteOutput_2DSection
Task: writes some output
Return:
Programming: XX/XXXX KB / CB
Modification:
-------------------------------------------------------------------------*/
void CECLIPSEData::WriteOutput_2DSection(long Timestep) {

    vector<string> pressure;
    vector<string> soil;
    vector<string> sgas;
    vector<string> swat;
    vector<string> oildens;
    vector<string> gasdens;
    vector<string> watdens;

    std::cout << " WriteOutput_2DSection()";

    std::string Filename, Filename_base, time;
    ofstream textfile;
    Filename_base = "Nusse_Out_2D_";

    // WTP to avoid getting warning due to conversion from long to string one could use a sstream
    //time.push_back(Timestep);
    std::stringstream strstream;
    strstream << Timestep;
    strstream >> time;

    Filename.clear();
    Filename = Filename_base + time + ".txt";

    textfile.open(Filename.data(), ios::out);

	textfile << "I" << "\t" << "J" << "\t" << "K" << "\t" << "x-barycentre" << "\t" << "y-barycentre" << "\t" << "z-barycentre";
	textfile << "\t" << "x1" << "\t" << "y1" << "\t" << "z1" << "\t" << "x2" << "\t" << "y2" << "\t" << "z2";
	textfile << "\t" << "x3" << "\t" << "y3" << "\t" << "z3" << "\t" << "x4" << "\t" << "y4" << "\t" << "z4";
	textfile << "\t" << "x5" << "\t" << "y5" << "\t" << "z5" << "\t" << "x6" << "\t" << "y6" << "\t" << "z6";
	textfile << "\t" << "x7" << "\t" << "y7" << "\t" << "z7" << "\t" << "x8" << "\t" << "y8" << "\t" << "z8";
    textfile << "\n";

	for (std::size_t i = 0; i < this->eclgrid.size(); i++)
	{
		//if (this->eclgrid[i]->row == 127)
		if (this->eclgrid[i]->column == 45)
        {
            textfile << this->eclgrid[i]->column << "\t" << this->eclgrid[i]->row << "\t" << this->eclgrid[i]->layer << "\t";
            textfile << scientific;
            textfile.precision(10);
            textfile << eclgrid[i]->x_barycentre << "\t" << eclgrid[i]->y_barycentre << "\t" << eclgrid[i]->z_barycentre << "\t";
			textfile << eclgrid[i]->x_coordinates[0] << "\t" << eclgrid[i]->y_coordinates[0] << "\t" << eclgrid[i]->z_coordinates[0] << "\t";
			textfile << eclgrid[i]->x_coordinates[1] << "\t" << eclgrid[i]->y_coordinates[1] << "\t" << eclgrid[i]->z_coordinates[1] << "\t";
			textfile << eclgrid[i]->x_coordinates[2] << "\t" << eclgrid[i]->y_coordinates[2] << "\t" << eclgrid[i]->z_coordinates[2] << "\t";
			textfile << eclgrid[i]->x_coordinates[3] << "\t" << eclgrid[i]->y_coordinates[3] << "\t" << eclgrid[i]->z_coordinates[3] << "\t";
			textfile << eclgrid[i]->x_coordinates[4] << "\t" << eclgrid[i]->y_coordinates[4] << "\t" << eclgrid[i]->z_coordinates[4] << "\t";
			textfile << eclgrid[i]->x_coordinates[5] << "\t" << eclgrid[i]->y_coordinates[5] << "\t" << eclgrid[i]->z_coordinates[5] << "\t";
			textfile << eclgrid[i]->x_coordinates[6] << "\t" << eclgrid[i]->y_coordinates[6] << "\t" << eclgrid[i]->z_coordinates[6] << "\t";
			textfile << eclgrid[i]->x_coordinates[7] << "\t" << eclgrid[i]->y_coordinates[7] << "\t" << eclgrid[i]->z_coordinates[7] << "\t";
            textfile << "\n";
        }
    }

    //Write general ouput file
    textfile << endl;
    textfile << "\n";
    textfile.close();
};

/**-------------------------------------------------------------------------
GeoSys - Function: SaveEclipseDataFile
Task: Saves the current .data file from Eclipse for benchmarking
The names of the .data files are enumerated with the oldest
having the lowest number (analog to the stored result files of Ecl).
Return: void
Programming: 04/2013 WTP
Modification:
-------------------------------------------------------------------------*/
void CECLIPSEData::SaveEclipseDataFile(long Timestep, CRFProcess* m_pcs) {
	std::string DataFilename_save;
	std::string file_extension;
	std::string projectname;
	std::string Filename;
	std::string systemcommand;
	std::string systemcommand2;
	std::string system_delete;
	std::string system_noquery;
	ostringstream temp;
	std::string extension;
	std::string laenge;
	CReadTextfiles_ECL *TextFile_in;
	CWriteTextfiles_ECL *TextFile_out;
	clock_t start, finish;
	double time;

	if (verbosity > 2)
		std::cout << "        SaveEclipseDataFile()";

	start = clock();

	// Getting the eclipse filename
	Filename = m_pcs->simulator_model_path;

	// check if filename is given with or without extension
	file_extension = Filename.substr(Filename.length() - 5, Filename.length());
	if ((file_extension.compare(".data") == 0) || (file_extension.compare(".DATA") == 0))
		projectname = Filename.substr(0, Filename.length() - 5);
	else
		projectname = Filename;
	if (Timestep > 1 || (m_pcs->iter_outer_cpl > 0))
		projectname = projectname + "_RESTART";

	//Set Eclipse model data filename
	DataFilename_save = projectname + "_" + AddZero(Timestep, 4, true) + file_extension;

	// Reading the original .data File
	TextFile_in = new CReadTextfiles_ECL;
	if (!TextFile_in->Read_Text(projectname + file_extension)) {
		// Copying the content of the original .data file to a new file
		TextFile_out = new CWriteTextfiles_ECL;
		TextFile_out->Write_Text(DataFilename_save, TextFile_in->Data);
	}
	else if (verbosity > 1)
		std::cout << " WARNING: Copying eclipse *.data files did not work properly" << "\n";

	finish = clock();
	time = (static_cast<double>(finish)-static_cast<double>(start)) / CLOCKS_PER_SEC;
	if (verbosity > 2)
		std::cout << "                 Time: " << time << " seconds." << "\n";

	return;
};

std::vector<std::vector<double> > CECLIPSEData::GetMolComponentsSorted(std::vector<std::vector<double> > vec_mComps) {
	//bool error = false;
	std::string keys[] = { "H2", "N2", "O2", "CH4", "CO2" };
	int idx[] = { -1, -1, -1, -1, -1 };
	// run thorugh vec_cnames_ecl and get the correct indicies,save them in idx
	// use the indicies in idx to sort vec_mComps_nodes
	for (int i = 0; i < 5; i++) // currently 5 components are supported
	{
		std::string current_comp = keys[i];
		// now loop over all available comps 
		for (unsigned int j = 0; j < Components.size(); j++)
		{
			std::string ecl_comps = Components[j];
			if (ecl_comps.compare(current_comp) == 0)
			{
				idx[i] = j;
				break;
			}
		}
	}
	// write mol fraction data into variable structure 
	for (unsigned int i = 0; i < this->NodeData.size(); i++)
	{
		std::vector <double> vec_dummy_dbl;
		double dummy_sum = 0.0;
		for (int j = 0; j < 5; j++)
		{
			const int current_idx = idx[j];
			if (current_idx == -1)
				vec_dummy_dbl.push_back(0.0);
			else
			{
				vec_dummy_dbl.push_back(this->NodeData[i]->mCompInGas_Nodes[current_idx]);
				dummy_sum += this->NodeData[i]->mCompInGas_Nodes[current_idx];
			}
		}
		for (int j = 0; j < 5; j++)
		{
			if (dummy_sum < 1.0E-10)
				vec_dummy_dbl[j] = 0.;
			else
				vec_dummy_dbl[j] = vec_dummy_dbl[j] / dummy_sum;
		}
		vec_mComps[i] = vec_dummy_dbl;
	}

	return vec_mComps;
};