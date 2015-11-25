/*
   rf_react_int.cpp
   Reaction interface
 */

#include "rf_react_int.h"
#include "display.h"
#include "timer.h"

#include "files0.h"
#include "makros.h"
#include "mathlib.h"
#include "rf_ic_new.h"
#include "rf_kinreact.h"
#include "rf_mmp_new.h"
#include "rf_pcs.h"
//#include "rf_react.h"
#include "rf_tim_new.h"
#include "rfmat_cp.h"
#include "stdio.h"
#include "tools.h"
#include <cfloat>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <signal.h>
#include <vector>
#include <sstream>
#include <IAPWS-IF97.h>


// Elem object
#include "fem_ele_std.h"

#include "Density.h"
#include "VLE.h"



#include <vector>
using namespace std;
using namespace Math_Group;
extern double gravity_constant;

vector <REACTINT*> REACTINT_vec;


static double dmaxarg1, dmaxarg2;
#define DMAX(a,b) (dmaxarg1 = (a),dmaxarg2 = (b),(dmaxarg1) > (dmaxarg2) ? \
  (dmaxarg1) : (dmaxarg2))
#define DMIN(a,b) (dmaxarg1 = (a),dmaxarg2 = (b),(dmaxarg1) < (dmaxarg2) ? \
  (dmaxarg1) : (dmaxarg2))




/**************************************************************************/
// REACINT class
/**************************************************************************/



/**************************************************************************
 ReacInt-Method:
 Task: Constructor 
 Programing:
 //CB 01.2011 CB First implementation
 **************************************************************************/
REACTINT::REACTINT(void){

  unitconversion = false;
  constantdensity = true;
  constanttemperature = true;
  constantpressure = true;
  c_TT = 298.15;  // std conditions 25°C
  c_PP = 1.0;     // std conditions 1 bar
  ssp_outstep = -1;
  pcs_outstep = -1;
  dump_min_moles = false;
  dump_all_pcs = false;
  dump_mass_integrals = false;
  s_water_limit = false;
  WaterSatLimit = 0.0;
  WaterSpeciesName = "NULL";
  NeutralCO2name = "NULL";
  SodiumSpeciesName = "NULL";
  icOutput = false;
  icSolidUpdate = false;
  readNodePoro = false;
  vle_flag = false;
  pcs_rename_init_flag = false;	
  pcs_rename_pre_flag = false;	
  pcs_rename_post_flag = false;	
  poroupdate_flag = false;
  heatpump_2DhTO2Dv = false;
  heatpump_Z = -9999;
    //0-Li, 1-Na, 2-K, 3-Mg, 4-Ca, 5-Cl, 6-SO4, 7-CO3
  //double mv[8];
  //mv[0]=0; mv[1]=2.71; mv[2]=0; mv[3]=0; mv[4]=0; mv[5]=2.71; mv[6]=0; mv[7]=0; 
  //double test = Density_CO2brine(393.15, 202, 2.71, 0.0);
  //test = Density_CO2_MultiBrine(393.15, 202, mv, 0);
  this->node_porosity.clear();
  this->node_ini_porosity.clear();
  this->water_conc.clear();
  this->water_conc_species.clear();
  this->sp_varind.clear();
  //this->sp_pcsind.clear();
  this->Temp_store.clear();
  Temp_GHP_mapidx.clear();

  //CB ReactDeact
  ReactDeactMode = -1;
  ReactDeactEpsilon = -1.0;                      
  ReactDeactCThresh = -1.0;
  ReactDeactRelative = false;
  ReactDeactFlag = false;
  ReactDeactPlotFlag = 0;
  ReactDeact.clear();                   // flags for individual nodes
  React_dCdT.clear();                   // Sum of reaction rates for individual nodes
  ReactNeighborhood.clear();            // node indices of local neighborhood around individual nodes

  // copy conc
  copy_concentrations = false; //SB 01.2011
  copy_nodes.clear();   //SB 01.2011
  radial = false;       // default case
  batch = false;
  TwoDinThreeD = false;

}
/**************************************************************************
 ReacInt-Method:
 Task: Destructor 
 Programing:
 //CB 01.2011 CB First implementation
 **************************************************************************/
REACTINT::~REACTINT(void){
}

/**************************************************************************
 ReacInt-Method:
 Task: general global Read function 
 Programing:
 //CB 01.2011 CB First implementation
 **************************************************************************/
bool REACINTRead(string file_base_name, const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name){

  char line[MAX_ZEILE];
  string sub_line;
  string line_string;
  string rei_file_name;
  ios::pos_type position;
  string m_name, sp_name;

  //========================================================================
  // File handling
  rei_file_name = file_base_name + REI_FILE_EXTENSION;
  ifstream rei_file (rei_file_name.data(),ios::in);
  if (!rei_file.good()) 
    return false;

  REACTINT *m_rei = NULL; 

  rei_file.seekg(0L,ios::beg); // rewind?
  //========================================================================
  // Keyword loop
  std::cout << "REIRead" << "\n" << std::flush;
  while (!rei_file.eof()) {
    rei_file.getline(line,MAX_ZEILE);
    line_string = line;
    if(line_string.find("#STOP")!=string::npos)
      break;
    //----------------------------------------------------------------------
    if(line_string.find("#REACTION_INTERFACE")!=string::npos) { // keyword found  // set up object
      m_rei = new REACTINT();
      m_rei->Read(&rei_file, geo_obj, unique_name);

      REACTINT_vec.push_back(m_rei);
    } // keyword found
  } // eof
  // Close input file
  rei_file.close();

  if(!m_rei){
	  cout << " No keyword #REACTION_INTERFACE specified - setting defaults" << "\n";
	  m_rei = new REACTINT();
      REACTINT_vec.push_back(m_rei);
  }
  return true;
}

/**************************************************************************
 ReacInt-Method:
 Task: ReacInt Read function 
 Programing:
 //CB 01.2011 CB First implementation
 **************************************************************************/
bool REACTINT::Read(ifstream* rfd_file, const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name)
{
  char line[MAX_ZEILE];
  string dummystring, speciesname;
  string line_string, line_str1;
  bool new_keyword = false, new_subkeyword = false;
  string hash("#"), dollar("$");
  std::stringstream in;
  long index, index1;
  double temp = 0;
  string dummy;

  //for pcs rename
  std::vector<std::string> pcs_name;
  std::vector<double>      pcs_stoi;
  std::string line_str2, ncomp_x;
  int j;
  VLE_type condition, condition_p;
  VLE_conditions.clear();
  VLE_pressure.clear();
		

  //========================================================================
  while (!new_keyword) {
	  index = rfd_file->tellg();
    if(!GetLineFromFile(line,rfd_file)) 
      break;
    line_string = line;
    if(line_string.find(hash)!=string::npos) {
      new_keyword = true;
  	  rfd_file->seekg(index); //Dateipointer zurücksetzen, sonst ist das nächste keyword weg
      break;
    }
	  /* Keywords nacheinander durchsuchen */
    //....................................................................
    if(line_string.find("$MOL_PER")!=string::npos) { // subkeyword found
	  in.str(GetLineFromFile1(rfd_file));
      in >> dummystring ;
      if(dummystring.compare("VOLUME")==0)
        unitconversion = true;
      else if(dummystring.compare("WEIGHT")==0)
        unitconversion = false;
      else
        cout << " Warning in REACTINT::Read - No valid keyword for $MOL_PER (unit conversions)." << "\n";
      in.clear();
    }
    if(line_string.find("$WATER_CONCENTRATION")!=string::npos) { // subkeyword found
	    in.str(GetLineFromFile1(rfd_file));
      in >> dummystring ;
      in.clear();
      if(dummystring.compare("CONSTANT")==0)
        constantdensity = true;
      else if(dummystring.compare("VARIABLE")==0){
        constantdensity = false;
        // Now read water concentration species list

        while(!new_subkeyword){
		      index1 = rfd_file->tellg();
		      line_str1 = GetLineFromFile1(rfd_file);
		      // Check for end of data block
	        if((line_str1.find(hash)!=string::npos) || (line_str1.find(dollar)!=string::npos)) {
   		      if(line_str1.find(hash)!=string::npos) 
              new_keyword = true;
   			    rfd_file->seekg(index1); //Dateipointer zurücksetzen, sonst ist das nächste subkeyword weg
	          break;
		      }
          // Here, read individual species data 
          in.str(line_str1);
          in >> speciesname ; 
          if(speciesname.compare("NIX")!=0){ //check for read in
            water_conc_species.push_back(speciesname);
            // store in osme vector 
            speciesname = "NIX";
          }
          else
			      DisplayMsgLn(" ERROR reading Water concentration relevant Species  - skipping");
		      in.clear();
	      }

      }
      else
        cout << " Warning in REACTINT::Read - No valid keyword for $WATER_CONCENTRATION." << "\n";
    }
    // water species name
    if(line_string.find("$WATER_SPECIES_NAME")!=string::npos) { // subkeyword found
	    in.str(GetLineFromFile1(rfd_file));
      in >> WaterSpeciesName ;
      in.clear();
    }
    // co2 species name
    if(line_string.find("$DISSOLVED_NEUTRAL_CO2_SPECIES_NAME")!=string::npos) { // subkeyword found
	    in.str(GetLineFromFile1(rfd_file));
      in >> NeutralCO2name ;
      in.clear();
    }
    if(line_string.find("$SODIUM_SPECIES_NAME")!=string::npos) { // subkeyword found
	    in.str(GetLineFromFile1(rfd_file));
      in >> SodiumSpeciesName ;
      in.clear();
    }
    if(line_string.find("$PRESSURE")!=string::npos) { // subkeyword found
	    in.str(GetLineFromFile1(rfd_file));
      in >> dummystring >> temp;
      if(dummystring.compare("CONSTANT")==0){
        constantpressure  = true;
        c_PP = temp;
      }
      else if(dummystring.compare("VARIABLE")==0)
        constantpressure  = false;
      else
        cout << " Warning in REACTINT::Read - No valid keyword for $PRESSURE." << "\n";
      in.clear();     
    }
    if(line_string.find("$TEMPERATURE")!=string::npos) { // subkeyword found
	    in.str(GetLineFromFile1(rfd_file));
      in >> dummystring >> temp;
      if(dummystring.compare("CONSTANT")==0){
        constanttemperature  = true;
        c_TT = temp;
      }
      else if(dummystring.compare("VARIABLE")==0)
        constanttemperature  = false;
      else
        cout << " Warning in REACTINT::Read - No valid keyword for $TEMPERATURE." << "\n";
	    in.clear();
    }
	// AB 12.05.2015: ... for geting salt concentration in the liquid 
	// hence CONCENTRATION1 refers only for delta salt in density driven study
	if (line_string.find("$INITIAL_SALT_CONCENTRATION") != string::npos) { // subkeyword found
		in.str(GetLineFromFile1(rfd_file));
		in >> c_salt;
		in.clear();
	}
    if(line_string.find("$WATER_SATURATION_LIMIT")!=string::npos) { // subkeyword found
	    in.str(GetLineFromFile1(rfd_file));
      in >> this->WaterSatLimit ;
      s_water_limit = true;
      in.clear();
    }
    if(line_string.find("$RESIDUAL")!=string::npos) { // subkeyword found
	    in.str(GetLineFromFile1(rfd_file));
      in >> dummystring ;
      if(dummystring.compare("RECORD")==0)
        residual_back = true;
      else if(dummystring.compare("REMOVE")==0)
        residual_back = false;
      else
        cout << " Warning in REACTINT::Read - No valid keyword for $RESIDUAL (for equilibrium reaction input value)." << "\n";
      in.clear();
    }
    if(line_string.find("$SOLID_SPECIES_DUMP_MOLES")!=string::npos) { // subkeyword found
	    in.str(GetLineFromFile1(rfd_file));
      in >> this->ssp_outstep ;
      dump_min_moles = true;
      in.clear();
    }
    if(line_string.find("$ALL_PCS_DUMP")!=string::npos) { // subkeyword found
	    in.str(GetLineFromFile1(rfd_file));
      in >> this->pcs_outstep ;
      dump_all_pcs = true;
      in.clear();
    }
    if(line_string.find("$INITIAL_CONDITION_OUTPUT")!=string::npos) { // subkeyword found
      icOutput = true;
    }
    if(line_string.find("$UPDATE_INITIAL_SOLID_COMPOSITION")!=string::npos) { // subkeyword found
      icSolidUpdate = true;
    }
	if(line_string.find("$VLE")!=string::npos){ // subkeyword found
		vle_flag = true;
        cout << " VLE load..." << "\n";		
		in.str(GetLineFromFile1(rfd_file));
		in >> ncomp_x;  condition.vp_name=ncomp_x;
		in >> ncomp_x; 	condition.aq_name=ncomp_x;
		in.clear();
		this->VLE_conditions.push_back(condition);
	}

	if(line_string.find("$P_VLE")!=string::npos){ // subkeyword found
		vle_p_flag = true;
        cout << " P_VLE load..." << "\n";		
		in.str(GetLineFromFile1(rfd_file));
		in >> ncomp_x;  condition_p.vp_name=ncomp_x;
		in >> ncomp_x; 	condition_p.aq_name=ncomp_x;
		in.clear();
		this->VLE_pressure.push_back(condition_p);
	}

    if(line_string.find("$POROSITY_RESTART")!=string::npos) { // subkeyword found
	    in.str(GetLineFromFile1(rfd_file));
      in >> this->porofile ;
      readNodePoro = true;
      in.clear();
    }
    if(line_string.find("$HEATPUMP_2DH_TO_2DV")!=string::npos) { // subkeyword found
	  in.str(GetLineFromFile1(rfd_file));
      in >> heatpump_Z;
      in.clear();
      heatpump_2DhTO2Dv = true;
    }

    //....................................................................
    // subkeyword found

    bool allow = false;
    bool allow_not = false;
    bool reallow = false;
    bool OK = true;
    string s_geo_type, s_geo_name;
    string dummy;




    allow = allow_not = reallow = false;
    if (line_string.find("$ALLOW_REACTIONS") != string::npos)
      allow = true;       // 1) AllowReactions: This allows reactions only on specific geometries
    if (line_string.find("$NO_REACTIONS") != string::npos)
      allow_not = true;   // 2) NoReactions, this suppresses reactions on specific geometries, can be combined with 1) 
    if (line_string.find("$RE_ALLOW_REACTIONS") != string::npos)
      reallow = true;       // 1) AllowReactions: This allows reactions only on specific geometries
    if ((allow == true) || (allow_not == true) || (reallow == true))
    {
      while (OK)
      {
        index1 = rfd_file->tellg();
        if (!GetLineFromFile(line, rfd_file))
          break;
        line_str1 = line;
        if ((line_str1.find(hash) != string::npos) || (line_str1.find(dollar) != string::npos))
        {
          OK = false;
          rfd_file->seekg(index1);           //Dateipointer zur�cksetzen, sonst ist das n�chste keyword weg
          break;
        }
        in.str(line_str1);
        in >> s_geo_type >> s_geo_name;
        if (allow_not == true){
          NoReactGeoType.push_back(s_geo_type);
          NoReactGeoName.push_back(s_geo_name);
        }
        else if (allow == true){
          AllowReactGeoType.push_back(s_geo_type);
          AllowReactGeoName.push_back(s_geo_name);
        }
        else if (reallow == true){
          ReAllowReactGeoType.push_back(s_geo_type);
          ReAllowReactGeoName.push_back(s_geo_name);
        }

        // TF 06/2010 - for the change string-ID to size_t-ID
        size_t geo_obj_idx(std::numeric_limits<size_t>::max());

        if (s_geo_type.find("POINT") != std::string::npos)
        {
          // get the point vector and set the geo_obj_idx
          if (!((geo_obj.getPointVecObj(unique_name))->getElementIDByName(s_geo_name, geo_obj_idx)))
          {
            std::cerr << "error in CKinReactData::Read: (type=" << s_geo_type << "): " << s_geo_name << " point name not found!" << "\n";
            exit(1);
          }
        }
        if (s_geo_type.find("POLYLINE") != std::string::npos)
        {
          // get the point vector and set the geo_obj_idx
          if (!((geo_obj.getPolylineVecObj(unique_name))->getElementIDByName(s_geo_name, geo_obj_idx)))
          {
            std::cerr << "error in CKinReactData::Read: polyline name " << s_geo_name << " not found!" << "\n";
            exit(1);
          }
        }

        //
        if (allow_not == true)
          NoReactGeoID.push_back(geo_obj_idx);
        else if (allow == true)
          AllowReactGeoID.push_back(geo_obj_idx);
        else if (reallow == true)
          ReAllowReactGeoID.push_back(geo_obj_idx);

        in.clear();
      }
    }


    if (line_string.find("$REACTION_DEACTIVATION") != string::npos) {
      in.str(GetLineFromFile1(rfd_file));
      in >> ReactDeactMode >> ReactDeactEpsilon >> ReactDeactPlotFlag >> dummy >> ReactDeactCThresh;
      in.clear();
      if (dummy.compare("RELATIVE") == 0)
        ReactDeactRelative = true;
      else
        ReactDeactRelative = false;
    }

    //....................................................................
    if (line_string.find("$COPY_CONCENTRATIONS") != string::npos) { //SB 01.2011
      std::cout << "COPY_CONCENTRATIONS allows copying concentrations after reactions to neighbouring"
        << " nodes for symmetric models:" << "\n";
      std::cout << " - BATCH: 2 node quasi 0D models." << "\n";
      std::cout << " - RADIAL symmetric quasi 2D hex models with y = 0 as symmetry axis." << "\n";
      std::cout << " - 2DIN3D: vertical x-z cross sections in full 3D hex models with arbitrary y as symmetry axis." << "\n";

      copy_concentrations = true;
      in.str(GetLineFromFile1(rfd_file));
      in >> dummy;
      in.clear();
      if (dummy.compare("RADIAL") == 0)
        radial = true; // default case
      else if (dummy.compare("BATCH") == 0)
        batch = true;
      else if (dummy.compare("2DIN3D") == 0)
        TwoDinThreeD = true;
      else {
        cout << "Error in CKinReactData::Read: $COPY_CONCENTRATIONS: No copy mode defined."
          << "\n" << " Define either RADIAL or BATCH or 2DIN3D. Terminating now!" << "\n";
        exit(1);
      }
    }


  //------------------------------------------------------

    //for pcs rename init, pre, post position
	if(line_string.find("$PCS_RENAME_INIT")!=string::npos){ // subkeyword found
  pcs_rename_init_flag = true;	
		this->pcs_rename0_init.clear();
		this->pcs_rename1_init.clear();
		this->pcs_rename1_stoi_init.clear();
		this->pow_stoi_init.clear();
		while(1){
			pcs_name.clear();
			pcs_stoi.clear();
			in.str(GetLineFromFile1(rfd_file));
			line_str2=in.str();
			in >> ncomp_x;
			if(!ncomp_x.find("END")) {
				in.clear();
				break;
			}			
			this->pcs_rename0_init.push_back(ncomp_x);
			for(j=1;j<(int)string2vector(line_str2).size();j+=2){
				in >> ncomp_x;
				if(ncomp_x=="*")
					pcs_stoi.push_back(999999);
				else if(ncomp_x=="/")
					pcs_stoi.push_back(-999999);
				else if(ncomp_x=="exp")
					pcs_stoi.push_back(-999998);
				else if(ncomp_x=="log")
					pcs_stoi.push_back(-999997);
				else if(ncomp_x=="log10")
					pcs_stoi.push_back(-999996);
				else if(ncomp_x=="delta")
					pcs_stoi.push_back(-999990);				
				else if(ncomp_x.substr(0,3)=="pow"){
					pcs_stoi.push_back(-999995);
					this->pow_stoi_init.push_back((double)atof(ncomp_x.substr(3).c_str()));
				}
				else
					pcs_stoi.push_back((double)atof(ncomp_x.c_str()));
				in >> ncomp_x;
				pcs_name.push_back(ncomp_x);
			}
			this->pcs_rename1_init.push_back(pcs_name);
			this->pcs_rename1_stoi_init.push_back(pcs_stoi);
			in.clear();				
		}
	}
	if(line_string.find("$PCS_RENAME_PRE")!=string::npos){ // subkeyword found
		pcs_rename_pre_flag = true;	
  this->pcs_rename0_pre.clear();
		this->pcs_rename1_pre.clear();
		this->pcs_rename1_stoi_pre.clear();
		while(1){
			pcs_name.clear();
			pcs_stoi.clear();
			in.str(GetLineFromFile1(rfd_file));
			line_str2=in.str();
			in >> ncomp_x;
			if(!ncomp_x.find("END")) {
				in.clear();
				break;
			}			
			this->pcs_rename0_pre.push_back(ncomp_x);
			for(j=1;j<(int)string2vector(line_str2).size();j+=2){
				in >> ncomp_x;
				if(ncomp_x=="*")
					pcs_stoi.push_back(999999);
				else if(ncomp_x=="/")
					pcs_stoi.push_back(-999999);
				else if(ncomp_x=="exp")
					pcs_stoi.push_back(-999998);
				else if(ncomp_x=="log")
					pcs_stoi.push_back(-999997);
				else if(ncomp_x=="log10")
					pcs_stoi.push_back(-999996);
				else if(ncomp_x=="delta")
					pcs_stoi.push_back(-999990);
				else if(ncomp_x.substr(0,3)=="pow"){
					pcs_stoi.push_back(-999995);
					this->pow_stoi_pre.push_back((double)atof(ncomp_x.substr(3).c_str()));
				}
				else
					pcs_stoi.push_back((double)atof(ncomp_x.c_str()));
				in >> ncomp_x;
				pcs_name.push_back(ncomp_x);
			}
			this->pcs_rename1_pre.push_back(pcs_name);
			this->pcs_rename1_stoi_pre.push_back(pcs_stoi);
			in.clear();				
		}
	}
	if(line_string.find("$PCS_RENAME_POST")!=string::npos){ // subkeyword found
  pcs_rename_post_flag = true;	
  this->pcs_rename0_post.clear();
		this->pcs_rename1_post.clear();
		this->pcs_rename1_stoi_post.clear();
		while(1){
			pcs_name.clear();
			pcs_stoi.clear();
			in.str(GetLineFromFile1(rfd_file));
			line_str2=in.str();
			in >> ncomp_x;
			if(!ncomp_x.find("END")) {
				in.clear();
				break;
			}			
			this->pcs_rename0_post.push_back(ncomp_x);
			for(j=1;j<(int)string2vector(line_str2).size();j+=2){
				in >> ncomp_x;
				if(ncomp_x=="*")
					pcs_stoi.push_back(999999);
				else if(ncomp_x=="/")
					pcs_stoi.push_back(-999999);
				else if(ncomp_x=="exp")
					pcs_stoi.push_back(-999998);
				else if(ncomp_x=="log")
					pcs_stoi.push_back(-999997);
				else if(ncomp_x=="log10")
					pcs_stoi.push_back(-999996);
				else if(ncomp_x=="delta")
					pcs_stoi.push_back(-999990);
				else if(ncomp_x.substr(0,3)=="pow"){
					pcs_stoi.push_back(-999995);
					this->pow_stoi_post.push_back((double)atof(ncomp_x.substr(3).c_str()));
				}
				else
					pcs_stoi.push_back((double)atof(ncomp_x.c_str()));
				in >> ncomp_x;
				pcs_name.push_back(ncomp_x);
			}
			this->pcs_rename1_post.push_back(pcs_name);
			this->pcs_rename1_stoi_post.push_back(pcs_stoi);
			in.clear();				
		}
	}
  }
  return true;
}



/**************************************************************************
 ReacInt-Method:
 Task: return the ReacInt object 
 Programing:
 //CB 01.2011 CB First implementation
 **************************************************************************/
REACTINT* REACTINT::GetREACTINT(void)
{
  REACTINT *rei = NULL;
  /* Tests */
  if(REACTINT_vec.capacity() > 0)
	rei = REACTINT_vec[0];
  return rei;
}

/**************************************************************************
 ReacInt-Method:
 Task: Initialize the ReacInt object
   - set initial poro and perm as element values of the flow pcs
   - prepare nodal poro vector
   - prepare water concentration vector
   - 
 Programing:
 //CB 01.2011 CB First implementation
 **************************************************************************/
void REACTINT::InitREACTINT(const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name){

  long i,j;
  double poro = 0;  
  int nc = (int)cp_vec.size();
  vector<int> id;
  string sp_name;
  int idx;//,idummy

  PrintGCCdata();



  if(icOutput)
    OUTData(0.0,aktueller_zeitschritt, true);

  //CompProperties *m_cp = NULL;
  CRFProcess *m_pcs = NULL;
  CFEMesh* m_msh = fem_msh_vector[0]; 
  if (m_msh == NULL)
  {
    cout << "No mesh in InitREACTINT" << "\n";
    exit(1);
  }
  const size_t nnodes = m_msh->nod_vector.size();

  m_pcs = PCSGetFlow();
  flowtype = m_pcs->type;

  /********************************************************/
  //Set up vector is_a_CCBC
  long ll, lm;
  std::vector<long> nodes_vector;
  std::string s_geo_name, s_geo_type;
  size_t s_geo_id;

  // Initialize vector is_a_CCBC
  if (is_a_CCBC.size() == 0)
    for (size_t l = 0; l < nnodes; l++)
      is_a_CCBC.push_back(false);

  // Go through specified geometry elements

  // 1) AllowReactions: This allows reactions only on specific geometries
  // First, set all nodes in vector is_a_CCBC inactive, if $ALLOW_REACTIONS was set
  if (AllowReactGeoName.size()>0)
    for (size_t l = 0; l < m_msh->nod_vector.size(); l++)
      is_a_CCBC[l] = true;

  // then, activate the nodes on geometries allowed for reactions
  for (size_t j = 0; j < AllowReactGeoName.size(); j++)
  {
    s_geo_name = AllowReactGeoName[j];
    s_geo_type = AllowReactGeoType[j];
    s_geo_id = AllowReactGeoID[j];
    //------------------------------------------------------------------
    if (s_geo_type.compare("POINT") == 0)
    {
      const std::vector<GEOLIB::Point*>* pnt_vec(geo_obj.getPointVec(unique_name));
      lm = m_msh->GetNODOnPNT((*pnt_vec)[AllowReactGeoID[j]]);
      is_a_CCBC[lm] = false;
    }
    //------------------------------------------------------------------
    if (s_geo_type.compare("POLYLINE") == 0)
    {

      //CGLPolyline *ply (polyline_vector[s_geo_id]);
      std::vector<GEOLIB::Polyline*> const* const ply_vec(geo_obj.getPolylineVec(unique_name));
      GEOLIB::Polyline const* const ply((*ply_vec)[s_geo_id]);
      if (ply)
      {
        m_msh->GetNODOnPLY(ply, nodes_vector);
        for (size_t i = 0; i < nodes_vector.size(); i++)
        {
          ll = nodes_vector[i];
          lm = ll;                         //+ShiftInNodeVector;
          is_a_CCBC[lm] = false;
        }
      }
    }                                        // if(POLYLINE)
    //------------------------------------------------------------------
    if (s_geo_type.compare("SURFACE") == 0)
    {
      Surface *m_surface = NULL;
      m_surface = GEOGetSFCByName(s_geo_name);
      if (m_surface)
      {
        m_msh->GetNODOnSFC(m_surface, nodes_vector); //kg44 if used together with PETSC one maybe needs to set the for_ic flag...see rf_id_new.cpp
        for (size_t i = 0; i < nodes_vector.size(); i++)
        {
          ll = nodes_vector[i];
          lm = ll;                         //+ShiftInNodeVector;
          is_a_CCBC[lm] = false;
        }
      }
    }
  }                                           // end of for(j=0;j<m_krd->NoReactGeoName.size()...
  //test output
  /*  cout << " Vector KinReactData::is_a_CCBC: " << "\n";
  for(l=0; l< (long)m_msh->nod_vector.size();l++)
  if(m_krd->is_a_CCBC[l] == true) cout << l <<", ";
  cout << "\n";
  */
  nodes_vector.clear();


  // 2) NoReactions, this suppresses reactions on specific geometries
  //    can be used in combination with 1)
  for (size_t j = 0; j < NoReactGeoName.size(); j++)
  {
    s_geo_name = NoReactGeoName[j];
    s_geo_type = NoReactGeoType[j];
    s_geo_id =   NoReactGeoID[j];
    //------------------------------------------------------------------
    if (s_geo_type.compare("POINT") == 0)
    {
      const std::vector<GEOLIB::Point*>* pnt_vec(geo_obj.getPointVec(unique_name));
      if (pnt_vec)
        is_a_CCBC[m_msh->GetNODOnPNT((*pnt_vec)[NoReactGeoID[j]])] = true;
    }
    //------------------------------------------------------------------
    if (s_geo_type.compare("POLYLINE") == 0)
    {
      CGLPolyline* ply(polyline_vector[s_geo_id]);

      std::vector<GEOLIB::Polyline*> const* const ply_vec(geo_obj.getPolylineVec(unique_name));
      GEOLIB::Polyline const* const polyline((*ply_vec)[s_geo_id]);
      double msh_min_edge_length(m_msh->getMinEdgeLength());
      m_msh->setMinEdgeLength(ply->epsilon);
      m_msh->GetNODOnPLY(polyline, nodes_vector);
      m_msh->setMinEdgeLength(msh_min_edge_length);
      for (size_t i = 0; i < nodes_vector.size(); i++)
        is_a_CCBC[nodes_vector[i]] = true;
      nodes_vector.clear();
    } // if(POLYLINE)
    //------------------------------------------------------------------
    if (s_geo_type.compare("SURFACE") == 0)
    {
      Surface* m_surface = NULL;
      m_surface = GEOGetSFCByName(s_geo_name);
      if (m_surface)
      {
        m_msh->GetNODOnSFC(m_surface, nodes_vector); //kg44 if used together with PETSC one maybe needs to set the for_ic flag...see rf_id_new.cpp
        for (size_t i = 0; i < nodes_vector.size(); i++)
        {
          long ll = nodes_vector[i]; //+ShiftInNodeVector;
          is_a_CCBC[ll] = true;
        }
      }
    }
  } // end of for(j=0;j<m_krd->NoReactGeoName.size()...
  //test output
  /*  cout << " Vector KinReactData::is_a_CCBC: " << "\n";
  for(l=0; l< (long)m_msh->nod_vector.size();l++)
  if(m_krd->is_a_CCBC[l] == true) cout << l <<", ";
  cout << "\n";
  */
  nodes_vector.clear();


  // 3) ReAllowReactions: This will switch on some of previously deactivated nodes on specific geometries
  for (size_t j = 0; j < ReAllowReactGeoName.size(); j++)
  {
    s_geo_name = ReAllowReactGeoName[j];
    s_geo_type = ReAllowReactGeoType[j];
    s_geo_id = ReAllowReactGeoID[j];
    //------------------------------------------------------------------
    if (s_geo_type.compare("POINT") == 0)
    {
      const std::vector<GEOLIB::Point*>* pnt_vec(geo_obj.getPointVec(unique_name));
      lm = m_msh->GetNODOnPNT((*pnt_vec)[ReAllowReactGeoID[j]]);
      is_a_CCBC[lm] = false;
    }
    //------------------------------------------------------------------
    if (s_geo_type.compare("POLYLINE") == 0)
    {

      //CGLPolyline *ply (polyline_vector[s_geo_id]);
      std::vector<GEOLIB::Polyline*> const* const ply_vec(geo_obj.getPolylineVec(unique_name));
      GEOLIB::Polyline const* const ply((*ply_vec)[s_geo_id]);
      if (ply)
      {
        m_msh->GetNODOnPLY(ply, nodes_vector);
        for (size_t i = 0; i < nodes_vector.size(); i++)
        {
          ll = nodes_vector[i];
          lm = ll;                         //+ShiftInNodeVector;
          is_a_CCBC[lm] = false;
        }
      }
    }                                        // if(POLYLINE)
    //------------------------------------------------------------------
    if (s_geo_type.compare("SURFACE") == 0)
    {
      Surface *m_surface = NULL;
      m_surface = GEOGetSFCByName(s_geo_name);
      if (m_surface)
      {
        m_msh->GetNODOnSFC(m_surface, nodes_vector); //kg44 if used together with PETSC one maybe needs to set the for_ic flag...see rf_id_new.cpp
        for (size_t i = 0; i < nodes_vector.size(); i++)
        {
          ll = nodes_vector[i];
          lm = ll;                         //+ShiftInNodeVector;
          is_a_CCBC[lm] = false;
        }
      }
    }
  }                                           // end of for(j=0;j<m_krd->NoReactGeoName.size()...
  //test output
  /*  cout << " Vector KinReactData::is_a_CCBC: " << "\n";
  for(l=0; l< (long)m_msh->nod_vector.size();l++)
  if(m_krd->is_a_CCBC[l] == true) cout << l <<", ";
  cout << "\n";
  */
  nodes_vector.clear();

//----------





  // for all mass transport processes, false: copy oldTL->newTL
  CopyAllConcentrationsToOtherTimeLevel(false);

  // AB & CB
  // Porosity & permeability preprocessing function
  SetInitialPorosityAndPermToElementValue();
  
  // set up vector for nodal porosities
  for(i=0;i<nnodes;i++){
    poro = GetNodePhaseVolume(i, 0, 0);
    node_porosity.push_back(poro);
    node_ini_porosity.push_back(poro);
  }
  
  if(this->readNodePoro)
    ReadRestartNodePoro(nnodes);

  // set up and initialize vector for water concentration
  CFluidProperties* m_mfp = NULL;
	 m_mfp = MFPGet( "LIQUID" );
  double variables[3];
  variables[0]=variables[1]=variables[2]=0;
  switch(m_mfp->density_model){
         case 0:                                  // rho = f(x)
            break;
         case 1:                                  // rho = const
            break;
         case 2:                                  // rho(p) = rho_0*(1+beta_p*(p-p_0))
           variables[0] = GetPressure(0);
            break;
         case 3:                                  // rho(C) = rho_0*(1+beta_C*(C-C_0))
           //variables[2] = GetConcentration(0);
           break;
         case 4:                                  // rho(T) = rho_0*(1+beta_T*(T-T_0))
            variables[1] = GetTemperature(0);            
            break;
         case 5:                                  // rho(C,T) = rho_0*(1+beta_C*(C-C_0)+beta_T*(T-T_0))
           variables[1] = GetTemperature(0);            
           //variables[2] = GetConcentration(0);
            break;
         case 6:                                  // rho(p,T) = rho_0*(1+beta_p*(p-p_0)+beta_T*(T-T_0))
           variables[0] = GetPressure(0);
           variables[1] = GetTemperature(0);            
            break;
         case 7:                                  // Pefect gas. WW
            break;
         case 8:                                  // M14 von JdJ // 25.1.12 Added by CB for density output AB-model
           variables[0] = GetPressure(0);
           variables[1] = GetTemperature(0);            
           //variables[2] = GetConcentration(0);
            break;
         default:
            std::cout << "Density in ReactInt: no valid model" << "\n";
            break;
  }
  double dens = m_mfp->Density(variables);
  if(unitconversion)
    for(i=0;i<nnodes;i++)
      water_conc.push_back(MOLH2OPERKG*dens);

  // set up vectors sp_pcs and sp_varind, taken from KinReactData
  //if(unitconversion){
    for(i=0;i<nc;i++) {
      sp_name = cp_vec[i]->compname;
      //idummy = -1;
      //m_pcs = PCSGet("MASS_TRANSPORT",sp_name);// CB HS update
      m_pcs = cp_vec[cp_name_2_idx[sp_name]]->getProcess();
      if(m_pcs != NULL){
        //idummy = PCSGetPCSIndex("MASS_TRANSPORT",sp_name);
        //sp_pcsind.push_back(idummy);
        idx = m_pcs->GetNodeValueIndex(sp_name)+1; // new timelevel
        sp_varind.push_back(idx);
      }
    }
  //}

  if(unitconversion ){
    if(constantdensity == false){
      for(j=0;j<10;j++) id.push_back(0);
      // set up and initialize matrix for water concentration computation
      for(i=0;i<nc;i++){ // set the species relevant for water concentration
        for(j=0;j<(int)id.size();j++) id[j]=0;
        // search in list, if species is relevant
        for(j=0;j<int(water_conc_species.size());j++) {
          if(cp_vec[i]->compname.compare(water_conc_species[j])==0 && cp_vec[i]->transport_phase==0 ) {// species found
            if(cp_vec[i]->iupac_formula[0]) // a formula has been defined
              id=formula2index(cp_vec[i]->iupac_formula);
            else // otherwise use species name directly. attention, this could result in wrong concentration of basis species!!!
              id=formula2index(cp_vec[i]->compname);
            break; // skip the rest of water_conc_species vector
          }
        }

        // Test
        if(cp_vec[i]->iupac_formula[0]) cout << cp_vec[i]->iupac_formula << ": " << "\n"; 
        else  cout << cp_vec[i]->compname << ": " << "\n"; 
        for(j=0;j<(int)id.size();j++) cout << " " << j ;    cout << "\n"; 
        for(j=0;j<(int)id.size();j++) cout << " " << id[j]; cout << "\n";

        // now store the id vector
        ElementInSpecies.push_back(id);
      }
    }
    // this should be called as well after transport and before any reaction calculations
    // necessary here as well for calculation of initial equilibrium
    CalcWaterConc();
  } // if(unitconversion)

  // prepare vector for reaction deactivation due to dry out from eclipse
  if(s_water_limit && WaterSatLimit>0)
     for(i=0;i<nnodes;i++) dried_out_nodes.push_back(false); // initialize for eclipse coupling
  
  // clean up
  id.clear();

  //DL 2012.2.17 for VLE conditions, set aq_idx and vp_idx
  if(vle_flag == true){
    cout << " VLE " << "\n";
    int cond;
    for(i=0;i<long(this->VLE_conditions.size());i++){
	    cond=0;
	    for(j=0;j<(int)pcs_vector.size();j++){
		    if(pcs_vector[j]->pcs_primary_function_name[0] == VLE_conditions[i].aq_name){
			    VLE_conditions[i].aq_idx=j;
			    cout << " " << VLE_conditions[i].aq_name;
			    cond=1;
		    }
	    }
	    if(cond==0){
		    cout << " VLE aq_name is not a pcs name " << "\n";
		    exit(0);
	    }		  
	    cond=0;
	    for(j=0;j<(int)pcs_vector.size();j++){
		    if(pcs_vector[j]->pcs_primary_function_name[0] == VLE_conditions[i].vp_name){
			    VLE_conditions[i].vp_idx=j;
			    cout << " " << VLE_conditions[i].vp_name;
			    cond=1;
		    }
	    }
	    if(cond==0){
		    cout << " VLE vp_name is not a pcs name " << "\n";
		    exit(0);
	    }
	    cout << "\n";
    }
  }


//DL 2013.5.10 for VLE pressure conditions, set aq_idx and vp_idx
  if(vle_p_flag == true){
  cout << " VLE_P " << "\n";
  int cond;
  for(i=0;i<(long)this->VLE_pressure.size();i++){
	  cond=0;
	  for(j=0;j<(int)pcs_vector.size();j++){
		  if(pcs_vector[j]->pcs_primary_function_name[0] == VLE_pressure[i].aq_name){
			  VLE_pressure[i].aq_idx=j;
			  cout << " " << VLE_pressure[i].aq_name;
			  cond=1;
		  }
	  }
	  if(cond==0){
		  cout << " VLE_P aq_name is not a pcs name " << "\n";
		  exit(0);
	  }		  
	  cond=0;
	  for(j=0;j<(int)pcs_vector.size();j++){
		  if(pcs_vector[j]->pcs_primary_function_name[0] == VLE_pressure[i].vp_name){
			  VLE_pressure[i].vp_idx=j;
			  cout << " " << VLE_pressure[i].vp_name;
			  cond=1;
		  }
	  }
	  if(cond==0){
		  cout << " VLE_P vp_name is not a pcs name " << "\n";
		  exit(0);
	  }
	  cout << "\n";
  }
  }
//====================================




//DL 2011.10.05  for pcs rename operation
  	int ii, ix, jx, is_idx, no_pcs;
	double value;
	vector<int> idx_tmp;
 no_pcs = (int)pcs_vector.size();
 
 if(pcs_rename_init_flag == true) {	
		  cout << "\n" << " PCS_RENAME_INIT " << "\n";	
	  this->pcs_rename0_idx_init.clear();
	  this->pcs_rename1_idx_init.clear();
	  for(j=0;j<(int)pcs_rename0_init.size();j++){
		  idx_tmp.clear();
		  for(jx=0;jx<(int)pcs_rename1_init[j].size();jx++)
			  idx_tmp.push_back(-1);
		  pcs_rename1_idx_init.push_back(idx_tmp);
		  pcs_rename0_idx_init.push_back(-1);
	  }
	  for(i=0;i<no_pcs;i++){
		  m_pcs = pcs_vector[i];
		  for(j=0;j<(int)pcs_rename0_init.size();j++){
			  if(m_pcs->pcs_primary_function_name[0] == pcs_rename0_init[j])
				  pcs_rename0_idx_init[j]=i;
			  for(jx=0;jx<(int)pcs_rename1_init[j].size();jx++)
				  if(m_pcs->pcs_primary_function_name[0] == pcs_rename1_init[j][jx])
			 		  pcs_rename1_idx_init[j][jx]=i;
		  }
	  }
	  for(i=0;i<(int)pcs_rename0_init.size();i++){
		  is_idx=1;
		  if(pcs_rename0_idx_init[i]==-1)
			  is_idx=0;
		  for(jx=0;jx<(int)pcs_rename1_init[i].size();jx++)
			  if(pcs_rename1_idx_init[i][jx]==-1)
				  is_idx=0;
		  if(is_idx==0){
			  cout << " Warning!!!, line " << i+1 << " can not be found in PCS list ! " << "\n";
			  exit(0);
		  }	
		  cout << " " << setw(12) << pcs_rename0_init[i] << " <--- ";
		  for(jx=0;jx<(int)pcs_rename1_init[i].size();jx++)
			  cout << "    "  << pcs_rename1_stoi_init[i][jx] << "  " << pcs_rename1_init[i][jx];
		  cout << "\n";
	  }
 }


 if(pcs_rename_pre_flag == true) {	
	  cout << "\n" << " PCS_RENAME_PRE " << "\n";	
	  this->pcs_rename0_idx_pre.clear();
	  this->pcs_rename1_idx_pre.clear();
	  for(j=0;j<(int)pcs_rename0_pre.size();j++){
		  idx_tmp.clear();
		  for(jx=0;jx<(int)pcs_rename1_pre[j].size();jx++)
			  idx_tmp.push_back(-1);
		  pcs_rename1_idx_pre.push_back(idx_tmp);
		  pcs_rename0_idx_pre.push_back(-1);
	  }
	  for(i=0;i<no_pcs;i++){
		  m_pcs = pcs_vector[i];
		  for(j=0;j<(int)pcs_rename0_pre.size();j++){
			  if(m_pcs->pcs_primary_function_name[0] == pcs_rename0_pre[j])
				  pcs_rename0_idx_pre[j]=i;
			  for(jx=0;jx<(int)pcs_rename1_pre[j].size();jx++)
				  if(m_pcs->pcs_primary_function_name[0] == pcs_rename1_pre[j][jx])
			 		  pcs_rename1_idx_pre[j][jx]=i;
		  }
	  }
	  for(i=0;i<(int)pcs_rename0_pre.size();i++){
		  is_idx=1;
		  if(pcs_rename0_idx_pre[i]==-1)
			  is_idx=0;
		  for(jx=0;jx<(int)pcs_rename1_pre[i].size();jx++)
			  if(pcs_rename1_idx_pre[i][jx]==-1)
				  is_idx=0;
		  if(is_idx==0){
			  cout << " Warning!!!, line " << i+1 << " can not be found in PCS list ! " << "\n";
			  exit(0);
		  }	
		  cout << " " << setw(12) << pcs_rename0_pre[i] << " <--- ";
		  for(jx=0;jx<(int)pcs_rename1_pre[i].size();jx++)
			  cout << "    "  << pcs_rename1_stoi_pre[i][jx] << "  " << pcs_rename1_pre[i][jx];
		  cout << "\n";
	  }
 }

 
 if(pcs_rename_post_flag == true) {	
   cout << "\n" << " PCS_RENAME_POST " << "\n";	
	  this->pcs_rename0_idx_post.clear();
	  this->pcs_rename1_idx_post.clear();
	  for(j=0;j<(int)pcs_rename0_post.size();j++){
		  idx_tmp.clear();
		  for(jx=0;jx<(int)pcs_rename1_post[j].size();jx++)
			  idx_tmp.push_back(-1);
		  pcs_rename1_idx_post.push_back(idx_tmp);
		  pcs_rename0_idx_post.push_back(-1);
	  }
	  for(i=0;i<no_pcs;i++){
		  m_pcs = pcs_vector[i];
		  for(j=0;j<(int)pcs_rename0_post.size();j++){
			  if(m_pcs->pcs_primary_function_name[0] == pcs_rename0_post[j])
				  pcs_rename0_idx_post[j]=i;
			  for(jx=0;jx<(int)pcs_rename1_post[j].size();jx++)
				  if(m_pcs->pcs_primary_function_name[0] == pcs_rename1_post[j][jx])
			 		  pcs_rename1_idx_post[j][jx]=i;
		  }
	  }
	  for(i=0;i<(int)pcs_rename0_post.size();i++){
		  is_idx=1;
		  if(pcs_rename0_idx_post[i]==-1)
			  is_idx=0;
		  for(jx=0;jx<(int)pcs_rename1_post[i].size();jx++)
			  if(pcs_rename1_idx_post[i][jx]==-1)
				  is_idx=0;
		  if(is_idx==0){
			  cout << " Warning!!!, line " << i+1 << " can not be found in PCS list ! " << "\n";
			  exit(0);
		  }	
		  cout << " " << setw(12) << pcs_rename0_post[i] << " <--- ";
		  for(jx=0;jx<(int)pcs_rename1_post[i].size();jx++)
			  cout << "    "  << pcs_rename1_stoi_post[i][jx] << "  " << pcs_rename1_post[i][jx];
		  cout << "\n";
	  }
 }

	for(i=0;i<no_pcs;i++){
		m_pcs = pcs_vector[i];
		if(m_pcs->getProcessType()==FiniteElement::MASS_TRANSPORT){
			this->nodenumber = (long) m_pcs->m_msh->GetNodesNumber(false);
			// elenumber = (long) m_pcs->m_msh->ele_vector.size();
		}
	}
 
 if(pcs_rename_init_flag == true){
	  int f=1, ic;
	  //return pcs rename init, do rename operation
	  for(ii=0;ii<this->nodenumber;ii++){ //node loop
		  ic=0;
		  for(i=0;i<(int)pcs_rename0_idx_init.size();i++){
			  value=0.0;
			  for(ix=0;ix<(int)pcs_rename1_idx_init[i].size();ix++){
				  m_pcs= pcs_vector[pcs_rename1_idx_init[i][ix]];
				  if(pcs_rename1_stoi_init[i][ix]==999999)
					  value *= m_pcs->GetNodeValue(ii,f);
				  else if(pcs_rename1_stoi_init[i][ix]==-999999)
					  value /= m_pcs->GetNodeValue(ii,f);
				  else if(pcs_rename1_stoi_init[i][ix]==-999998)
					  value += exp(m_pcs->GetNodeValue(ii,f));
				  else if(pcs_rename1_stoi_init[i][ix]==-999997)
					  value += log(m_pcs->GetNodeValue(ii,f));
				  else if(pcs_rename1_stoi_init[i][ix]==-999996)
					  value += log10(m_pcs->GetNodeValue(ii,f));
				  else if(pcs_rename1_stoi_init[i][ix]==-999990)
					  value += m_pcs->GetNodeValue(ii,f)-m_pcs->GetNodeValue(ii,0);
				  else if(pcs_rename1_stoi_init[i][ix]==-999995){
					  value += pow(m_pcs->GetNodeValue(ii,f), this->pow_stoi_init[ic]);
					  ic ++;
                  }
				  else
					  value += m_pcs->GetNodeValue(ii,f) *pcs_rename1_stoi_init[i][ix];
              }
			  m_pcs= pcs_vector[pcs_rename0_idx_init[i]];
			  m_pcs->SetNodeValue(ii,f,value);
		  }
	  }
 }

 this->t_step=0;


 /********************************************************/
 // CB Set neighborhood relations for reaction deactivation switch
 // CB reaction deactivation
 //int annode_idx, lll, llll, lllll, nnodpneigh, duplicate,nn, nnodsneigh 
 //int dnele_idx; // ,  pneighnod_idx, snele_idx,
 //  int dnnode, snnode, duplicate2;
 MeshLib::CElem* m_dnele = NULL;
 MeshLib::CElem* m_snele = NULL;
 MeshLib::CNode* m_dnnod = NULL;
 MeshLib::CNode* m_nod = NULL;
 MeshLib::CElem* m_ele = NULL;
 std::vector<int> ReactNeighborNodes;
 vec<long> secnnodesindices(8);
 vec<long> primnnodesindices(8);

 if ((ReactDeactEpsilon) < 0)
   ReactDeactFlag = false;
 else
   ReactDeactFlag = true;
 if (ReactDeactFlag)
 {
   for (size_t l = 0; l < nodenumber; l++)
   {
     // initialize all nodes as active
     ReactDeact.push_back(false);
     // initialize all reaction terms as zero
     React_dCdT.push_back(0);
     m_nod = m_msh->nod_vector[l]; // get current node object

     // loop over PRIMARY NEIGHBOR elements
     // this routine takes longer, but gets all the nodes
     for (size_t ll = 0; ll < m_nod->getConnectedElementIDs().size();ll++)
     {
       // get index of direct neighbour element
       size_t dnele_idx = m_nod->getConnectedElementIDs()[ll];
       // get direct neighbour element object
       m_dnele = m_msh->ele_vector[dnele_idx];
       // get the neighbor element node indices
       m_dnele->GetNodeIndeces(primnnodesindices);
       // get the neighbor element number of nodes
       int nnodpneigh = m_dnele->GetNodesNumber(false);
       // loop over primary neighbour element number of nodes
       for (int lll = 0; lll < nnodpneigh; lll++)
       {
         // get the current node index
         long pneighnod_idx = primnnodesindices[lll];
         // get the node object
         m_dnnod = m_msh->nod_vector[pneighnod_idx];
         // loop over the connected elements, this now includes the SECONDARY NEIGHBORS
         for (size_t llll = 0;llll < m_dnnod->getConnectedElementIDs().size();llll++)
         {
           // get the current connected element indices
           size_t snele_idx = m_dnnod->getConnectedElementIDs()[llll];
               // get secondary neighbour element object
               m_snele = m_msh->ele_vector[snele_idx];
               // get the neighbor element node indices
               m_snele->GetNodeIndeces(secnnodesindices);
               // get the neighbor element number of nodes
               int nnodsneigh = m_snele->GetNodesNumber(false);
               // loop over secondary neighbour element number of nodes
               for (int lllll = 0; lllll < nnodsneigh; lllll++)
               {
                 bool duplicate(false); // flag, initialize
                 int annode_idx = secnnodesindices[lllll];
                 // check vs. all previously identified nodes to prevent duplicate entries
                 for (size_t nn = 0; nn < ReactNeighborNodes.size(); nn++) // node has already been found
                   if (annode_idx == ReactNeighborNodes[nn])
                   {
                     duplicate = true; // set flag to "not fresh"
                     break; // skip rest of comparisons
                   }
                 if (!duplicate) // fresh node
                   // push back node index
                   ReactNeighborNodes.push_back(annode_idx);
                 else
                   duplicate = false;  // was no fresh node, reinitialize flag
               }
         }
       }
     }
     // Add local node neighborhood indices vector to global vector
     // This is the same for both above neighborhood search routines
     ReactNeighborhood.push_back(ReactNeighborNodes);
     ReactNeighborNodes.clear();
   }
   // debug
   //for (nn=0;nn<ReactNeighborhood.size();nn++){
   //  cout << nn << " " << ReactNeighborhood[nn].size() << " " ;
   //  for (ll=0;ll<ReactNeighborhood[nn].size();ll++)
   //    cout << ReactNeighborhood[nn][ll] << " " ;
   //  cout << "\n";
   //}
   concentrationmatrix = (double**)malloc(((long)nodenumber) * sizeof(double*));
   for (size_t ix = 0; ix < nodenumber; ix++)
     concentrationmatrix[ix] = (double*)malloc((nc)* sizeof(double));
   for (size_t ix = 0; ix < nodenumber; ix++)
     for (size_t ixx = 0; ixx < nc; ixx++)
       concentrationmatrix[ix][ixx] = 0;
 } // if (ReactDeactFlag)

}

/**************************************************************************
 ReacInt-Method:
 Task: General Preprocessing for reactions
   - calculate water concentration at all nodes
 Programing:
 //CB 01.2011 CB First implementation
 **************************************************************************/
void REACTINT::ReactionPreProcessing(void){
	if(unitconversion)
	  	CalcWaterConc();
 if(s_water_limit && WaterSatLimit>0)
    CheckForDriedOutNodes();

	
 // CB Reaction deactivation for this time step
 if (ReactDeactFlag){
   if (aktueller_zeitschritt > 2)
     ReactionDeactivation(nodenumber);      // Check if nodes should be deactivated or activated for this time step
   if (ClockTimeVec.size()>0){
     ClockTimeVec[0]->StopTime("ReactDeact");  // CB time
     ClockTimeVec[0]->StartTime();  // CB time
   }
 }

 
 
 //return pcs rename pre, do rename operation
	int i, ii, ix, f=1, ic;
	double value;
	CRFProcess *m_pcs = NULL;
 if (pcs_rename_pre_flag == true){
   for(ii=0;ii<this->nodenumber;ii++){ //node loop
		   ic=0;
		   for(i=0;i<(int)pcs_rename0_idx_pre.size();i++){
			    value=0.0;
			    for(ix=0;ix<(int)pcs_rename1_idx_pre[i].size();ix++){
				    m_pcs= pcs_vector[pcs_rename1_idx_pre[i][ix]];
				    if(pcs_rename1_stoi_pre[i][ix]==999999)
					    value *= m_pcs->GetNodeValue(ii,f);
				    else if(pcs_rename1_stoi_pre[i][ix]==-999999)
					    value /= m_pcs->GetNodeValue(ii,f);
				    else if(pcs_rename1_stoi_pre[i][ix]==-999998)
					    value += exp(m_pcs->GetNodeValue(ii,f));
				    else if(pcs_rename1_stoi_pre[i][ix]==-999997)
					    value += log(m_pcs->GetNodeValue(ii,f));
				    else if(pcs_rename1_stoi_pre[i][ix]==-999996)
					    value += log10(m_pcs->GetNodeValue(ii,f));
				    else if(pcs_rename1_stoi_pre[i][ix]==-999990)
					    value += m_pcs->GetNodeValue(ii,f)-m_pcs->GetNodeValue(ii,0);
				    else if(pcs_rename1_stoi_pre[i][ix]==-999995){
					    value += pow(m_pcs->GetNodeValue(ii,f), this->pow_stoi_pre[ic]);
					    ic ++;
				    }
				    else
					    value += m_pcs->GetNodeValue(ii,f) *pcs_rename1_stoi_pre[i][ix];
			    }
			    m_pcs= pcs_vector[pcs_rename0_idx_pre[i]];
			    m_pcs->SetNodeValue(ii,f,value);
		   }
	  }
 }
	return;
}

/**************************************************************************
 ReacInt-Method:
 Task: General Postprocessing for reactions
   - calculate porosity update from mineral reactions
   - calculate permeability update from porosity change
 Programing:
 //CB 01.2011 CB First implementation
 **************************************************************************/
void REACTINT::ReactionPostProcessing(bool initial){
	 
  int i, ii, ix, f=1, ic;
  double value;
  CRFProcess *m_pcs = NULL;

  cout << " ReactionPostProcessing" << "\n";

  if(initial) // only after initial reactions, i.e. before first real time step
    this->CopyAllConcentrationsToOtherTimeLevel(true); 
  else {
  
    // CB Reaction deactivation for next time step
    if (ReactDeactFlag) {
      //cout << "    Kinetic reactions executed at " << count << " of " << nnodes << " nodes." << "\n";
      if (ReactDeactMode == 2) // For mode 2 the C_new must be updated by C_old (C after eraction of last time step)
        ReactDeactSetOldReactionTerms(nodenumber);
      if (ReactDeactPlotFlag == 1)
        ReactDeactPlotFlagsToTec();
    }
    

    
    //  CopySymmetricConcentrationsInRadialModel();
    if(poroupdate_flag==true){
      PorosityVolumetricReactionUpdate(); // porosity change
      PermeabilityPorosityUpdate();       // resulting permeability change 
    }       
    if(dump_min_moles==true)
      DumpSolidSpeciesMoles();
    if(dump_all_pcs==true)
      DumpAllVariables();
    if(dump_mass_integrals==true)
      DumpMassIntegrals();

    //return pcs rename post, do rename operation
    if(this->pcs_rename_post_flag==true){
      for(ii=0;ii<this->nodenumber;ii++){ //node loop
		      ic=0;
		      for(i=0;i<(int)pcs_rename0_idx_post.size();i++){
			       value=0.0;
			       for(ix=0;ix<(int)pcs_rename1_idx_post[i].size();ix++){
				        m_pcs= pcs_vector[pcs_rename1_idx_post[i][ix]];
				        if(pcs_rename1_stoi_post[i][ix]==999999)
					         value *= m_pcs->GetNodeValue(ii,f);
				        else if(pcs_rename1_stoi_post[i][ix]==-999999)
					         value /= m_pcs->GetNodeValue(ii,f);
				        else if(pcs_rename1_stoi_post[i][ix]==-999998)
					         value += exp(m_pcs->GetNodeValue(ii,f));
				        else if(pcs_rename1_stoi_post[i][ix]==-999997)
					         value += log(m_pcs->GetNodeValue(ii,f));
				        else if(pcs_rename1_stoi_post[i][ix]==-999996)
					         value += log10(m_pcs->GetNodeValue(ii,f));
	      			  else if(pcs_rename1_stoi_post[i][ix]==-999990)
					         value += m_pcs->GetNodeValue(ii,f)-m_pcs->GetNodeValue(ii,0);
				        else if(pcs_rename1_stoi_post[i][ix]==-999995){
					         value += pow(m_pcs->GetNodeValue(ii,f), this->pow_stoi_post[ic]);
					         ic ++;
            }
				        else
					         value += m_pcs->GetNodeValue(ii,f) *pcs_rename1_stoi_post[i][ix];
          }
			       m_pcs= pcs_vector[pcs_rename0_idx_post[i]];
			       m_pcs->SetNodeValue(ii,f,value);
        }
      }
    }
  }

  return;
}

/**************************************************************************
 ReacInt-Method:
 Task: not used at the moment
 Programing:
 //CB 01.2011 CB First implementation
 **************************************************************************/
double REACTINT::CalcDensityOfWater(void){

  double dens = 0;
  /*
  // Test density function water
  // Li Na K Mg Ca Cl SO4 CO3
  mv[0]=0;
  mv[1]=species_value[6];
  mv[2]=0;
  mv[3]=species_value[9];
  mv[4]=species_value[8];
  mv[5]=species_value[7];
  mv[6]=0;
  mv[7]=species_value[5];
  mCO2 = 0;

  t0=double(clock());
  dens = Density_CO2_MultiBrine(TT,PP,mv,mCO2);
  t1+=double(clock())-t0;
  cout << "\n" << dens ;
  //cout <<  double(t1/CLOCKS_PER_SEC) << "\n";
  */
  return dens;
}

/**************************************************************************
 ReacInt-Method:
 Task: Return Pressure in bar
 Programing:
 //CB 01.2011 CB First implementation
 **************************************************************************/
double REACTINT::GetPressure(long node){

  double  PP = 1;
  double factor = 1, z;
  int	idxp=0, idxp2=0, timelevel = 1;
  bool ppmodel = false;
  
  CRFProcess *m_pcs_p = NULL;
  CFEMesh* m_msh = fem_msh_vector[0]; 
  CFluidProperties* m_mfp = NULL;

  if(constantpressure) 
    return c_PP;

	m_mfp = MFPGet( "LIQUID" );


  if(aktueller_zeitschritt==0) timelevel = 0;
  // get the Pressure info Todo: check this function for correct pressure conversion
  m_pcs_p = PCSGetFlow();

  if (flowtype==1)
    if(m_pcs_p->getProcessType()==FiniteElement::LIQUID_FLOW)
      flowtype=0;

  switch(flowtype){
    case 0: //Liquid_Flow
      idxp = m_pcs_p->GetNodeValueIndex ("PRESSURE1") + timelevel;
      factor = 1 / 1.0e+5 ;
      break;
    case 1: //Groundwater Flow
      idxp = m_pcs_p->GetNodeValueIndex ("HEAD") + timelevel;
      factor = gravity_constant * m_mfp->Density() / 1.0e+5 ;
      break;
    case 66: //Overland Flow
      break;
    case 5: //Air Flow
      break;
    case 11: //Componental Flow
      break;
    case 1212: //Multiphase Flow
      idxp = m_pcs_p->GetNodeValueIndex ("PRESSURE1") + timelevel; // p_cap
      idxp2 = m_pcs_p->GetNodeValueIndex ("PRESSURE2") + timelevel;// p_nw
      factor = 1 / 1.0e+5 ;
      ppmodel = true;
      break;
    case 12: // Two_phase_Flow
      idxp = m_pcs_p->GetNodeValueIndex ("PRESSURE1") + timelevel;
      factor = 1 / 1.0e+5 ;
      break;
    case 1313: // PS_GLOBAL
      idxp = m_pcs_p->GetNodeValueIndex ("PRESSURE1") + timelevel;
      factor = 1 / 1.0e+5 ;
      break;    
    case 22: // Richards flow
      idxp = m_pcs_p->GetNodeValueIndex ("PRESSURE1") + timelevel;
      factor = 1 / 1.0e+5 ;
      break;
    default:
      break;
  }
  
  // get pressure    todo CB
  if(flowtype==1){ //Groundwater flow head --> Pressure
    z = m_msh->nod_vector[node]->getData()[2]; // CB_merge_0513 ?? check
    PP = (m_pcs_p->GetNodeValue (node, idxp) - z) * factor;
  }
  else
    PP = m_pcs_p->GetNodeValue (node, idxp) * factor;
  if(ppmodel) // add nw pressure to capillary pressure
    PP += m_pcs_p->GetNodeValue (node, idxp2) * factor;
  if(PP<=0 || PP > 1000){ // limit pressure
    cout << "Warning in REACTINT::GetPressure(): unreasonable P: " << PP ;
    cout << " bar; limiting to P = 1 bar" << "\n";
    PP = 1;
  }
  return PP; // bar
}

/**************************************************************************
 ReacInt-Method:
 Task: Calculate water concentration at all nodes
 Programing:
 //CB 01.2011 CB First implementation
 **************************************************************************/
double REACTINT::GetTemperature(long node){

  double TT=298.15;
  int	idxt=0, timelevel = 1;
//  bool heattransport = false;
  
  if(constanttemperature) 
    return c_TT;

  if(aktueller_zeitschritt<1) timelevel = 0;
//  if(aktueller_zeitschritt<=1) timelevel = 0; // CB-DL??
  CRFProcess *m_pcs_t = NULL;

  // get the Temperature info
  if((m_pcs_t = PCSGet("HEAT_TRANSPORT"))) {
	idxt = m_pcs_t->GetNodeValueIndex ("TEMPERATURE1") + timelevel;
    TT = m_pcs_t->GetNodeValue (node, idxt);
  }
  else 
    TT = c_TT;

  return TT;
}

/**************************************************************************
 ReacInt-Method:
 Task: Check all nodes, if water sat is lower than WaterSatLimit
 Programing:
 //CB 11.2011 CB First implementation
 **************************************************************************/
void REACTINT::CheckForDriedOutNodes(void){
  
  long i, nnode;

  CFEMesh* m_msh = fem_msh_vector[0]; 
  nnode = (long) m_msh->nod_vector.size();

  if(aktueller_zeitschritt>1){
    for (i=0;i<nnode;i++){
      if(GetWaterSaturation(i) < WaterSatLimit) 
        dried_out_nodes[i] = true;
      else 
        dried_out_nodes[i] = false;
    }
  }    

  return;
}


/**************************************************************************
 ReacInt-Method:
 Task: For all species, copy node values to another time level
 old->new : forward = false
 new->old : forward = true
 Programing:
 //CB 10.2011 CB First implementation
 **************************************************************************/
void REACTINT::CopyAllConcentrationsToOtherTimeLevel(bool forward){

  int i;
  //CRFProcess *m_pcs = NULL;
  
  for (i = 0; i < int(cp_vec.size()); i++){
    cp_vec[i]->getProcess()->CopyTimestepNODValues(forward); // true: new->old

  // test
  //for (int j = 0; j < 5; j++){
  //  std::cout << cp_vec[i]->getProcess()->GetNodeValue(j, 0) << "     ";
  //  std::cout << cp_vec[i]->getProcess()->GetNodeValue(j, 1) << "\n";
  //}
  //std::cout << "\n";
}


  return;
}

/**************************************************************************
 ReacInt-Method:

 Programing:
 //CB 10.2011 CB First implementation
 **************************************************************************/
void REACTINT::Heatpump_2DhTO2Dv_Mapping(bool forward){

  int	idxt=0, timelevel = 1;
  long j=0;
  long i=0;
  double const *coord; //[3]={0,0,0};

   
  CRFProcess *m_pcs_t = NULL;
  long nnode = long(fem_msh_vector[0]->nod_vector.size());

  // get the Temperature info
  if((m_pcs_t = PCSGet("HEAT_TRANSPORT")))
	  idxt = m_pcs_t->GetNodeValueIndex("TEMPERATURE1") + timelevel;
  else
    return;
 
  // backward mapping 
  if(!forward){ 
    if(aktueller_zeitschritt==1)
      return;
    else 
      for(i=0;i<nnode;i++){
        m_pcs_t->SetNodeValue(i, idxt, Temp_store[i]);   // new TL
        m_pcs_t->SetNodeValue(i, idxt-1, Temp_store[i]); // old TL
      }
  }

  // forward mapping
  
  //check output
  //if(aktueller_zeitschritt==25) return;
  
  std::vector<std::vector<double> > Temp_GHP_xyz;
  std::vector<std::vector<double> > XYZ;
  std::vector <long> Temp_GHP_nod_idx;
  std::vector <double> xyz;  

  if(aktueller_zeitschritt==1){
    for(i=0;i<3;i++) xyz.push_back(0);
    for(i=0;i<nnode;i++){
      // first save all Temperatures
      Temp_store.push_back(m_pcs_t->GetNodeValue(i, idxt));
      // now store mapping relations
      coord = fem_msh_vector[0]->nod_vector[i]->getData();
      for(j=0;j<3;j++) xyz[j]=coord[j];
      // find nodes on plume center line and store xyz
      if(coord[2]==heatpump_Z) {
        Temp_GHP_xyz.push_back(xyz);
        Temp_GHP_nod_idx.push_back(i);
      }
      // store all other node corrdinates
      XYZ.push_back(xyz);
    }
    // now check for each node the corresponding map node on center line
    for(i=0;i<nnode;i++){
      for(j=0;j<long(Temp_GHP_xyz.size());j++)
        if((XYZ[i][0]==Temp_GHP_xyz[j][0])) // x cordinate matches -> only regular mesh so far
          Temp_GHP_mapidx.push_back(Temp_GHP_nod_idx[j]);
    }
  }
  else{ // after 1st time step
    for(i=0;i<nnode;i++)
      Temp_store[i]= m_pcs_t->GetNodeValue(i, idxt);  
  }

  // now map center line to xz plane
  for(i=0;i<nnode;i++)
    m_pcs_t->SetNodeValue(i, idxt, Temp_store[Temp_GHP_mapidx[i]]);

  // clean up
  xyz.clear(); 
  Temp_GHP_xyz.clear();
  Temp_GHP_nod_idx.clear();
  XYZ.clear();

  return;

}



void REACTINT::ReadRestartNodePoro(long nnodes){

  long i,j;
  char zeile[1000];
  double porosity;
  // get medium property
//  CMediumProperties *m_mat_mp = NULL;
  MeshLib::CElem* m_ele = NULL;
	CFEMesh* m_msh = fem_msh_vector[0];
  long n_ele = long(m_msh->ele_vector.size());
  CRFProcess* m_pcs_flow = NULL;
	  m_pcs_flow = PCSGetFlow();
  int idx_n = m_pcs_flow->GetElementValueIndex("POROSITY");
int node_idx;


  ifstream ein (this->porofile.data(),ios::in);
  if (!ein.good()){
	   cout << " Input file not found in ReadRestartNodePoro" << "\n";
	   exit(0);
  }
  
//  ein.str(GetLineFromFile1(rfd_file));
  for (i=0;i< nnodes;i++){
    ein >> this->node_porosity[i];
    ein.getline(zeile,1000);
  }

  // now set the new element porosities from node_porosity for mass transport and other processes
  // loop over all elements
  for(long n = 0; n < n_ele; n++){
    m_ele = m_msh->ele_vector[n];
    long nn = m_ele->GetNodesNumber(false); // nodes_index.Size();
    porosity  = 0;
    for(j = 0; j < nn; j++){
        node_idx = m_ele->GetNodeIndex(j);
        porosity  += node_porosity[node_idx];
    }
    //arithmetic average ... inverse distance?
    porosity  /= (double) nn;
    //set the element value
    m_pcs_flow->SetElementValue( n, idx_n+1, porosity );
  }

  return;
}



/**************************************************************************
 ReacInt-Method:
 Task: Calculate water concentration at all nodes
 Programing:
 //CB 01.2011 CB First implementation
 **************************************************************************/
void REACTINT::CalcWaterConc(void){

  long i,j,k=0;
  double cv[8], cv_CO2, TT=298.15, PP = 1;
  double conc;//factor = 1, 
  int	idxCO2=0, widx=-1; // timelevel = 1, idxp=0, idxp2=0, idxt=0, 
  //bool heattransport = false;
  //bool ppmodel = false;
  bool CO2 = false;
  //bool timeflag  = false;
  int nc = (int)cp_vec.size();
  
  CompProperties *m_cp = NULL;
  //CRFProcess *m_pcs_t = NULL;
  //CRFProcess *m_pcs_p = NULL;
  CFEMesh* m_msh = fem_msh_vector[0]; 
  //MeshLib::CNode* m_nod = NULL;
  long nnodes = (long) m_msh->nod_vector.size();
  CFluidProperties* m_mfp = NULL;
 	m_mfp = MFPGet( "LIQUID" );
  //double dens = 0;

  //if(aktueller_zeitschritt==0) timelevel = 0;

  // Get water concentration index
  for(j=0;j<nc;j++){
    m_cp = cp_vec[j];    
    if(m_cp->compname.compare(WaterSpeciesName)==0 || m_cp->compname.compare("H2O")==0 || m_cp->compname.compare("H2O_liquid")==0  || m_cp->compname.compare("water_liquid")==0){
      widx = j;
      break;
    }
  }

  if(constantdensity) {   
    for(i=0;i<nnodes;i++){
      water_conc[i] =   m_mfp->Density() * MOLH2OPERKG;  // here density of pure H2O phase must be used
      //water_conc[i] =   4.54912323E+04;
      if(widx > -1) 
        //pcs_vector[sp_pcsind[widx]]->SetNodeValue(i,1,water_conc[i]);
        cp_vec[widx]->getProcess()->SetNodeValue(i,1,water_conc[i]); 
    }
    return;
  }

  // Get CO2_in_water concentration index
  // this should be CO2_dissolved, not total CO2, 
  // so get secondary species CO2, not C(4) if you are using Phreeqc
  for(j=0;j<nc;j++){
    m_cp = cp_vec[j];    
    if(m_cp->compname.compare(this->NeutralCO2name)==0 ||  m_cp->compname.compare("CO2_w")==0 || m_cp->compname.compare("CO2")==0 ){ 
      CO2 = true;
      idxCO2 = j;
      break;
    }
  }
  
  //if(aktueller_zeitschritt==0 || aktueller_zeitschritt==1) timeflag = true;

  // node loop to calculate water concentration
  for(i=0;i<nnodes;i++){
    TT = this->GetTemperature(i);
    PP = this->GetPressure(i);
    //Weight function of water
    // Li Na K Mg Ca Cl SO4 CO3 (only charged CO2 species)
    for(k=0;k<8;k++) cv[k]= 0;
    // loop over all species for current node
    for(j=0;j<nc;j++){
      // get concentration; this requires existence of vectors sp_pcsind and sp_var_ind
      // todo: what about double species, e.g. primary and secondary species??
      //conc = pcs_vector[sp_pcsind[j]]->GetNodeValue(i,sp_varind[j]);// CB HS update
      //if(timeflag &! cp_vec[j]->mobil)
      //  conc = cp_vec[j]->getProcess()->GetNodeValue(i,sp_varind[j]-1);
      //else 
        conc = cp_vec[j]->getProcess()->GetNodeValue(i,sp_varind[j]); // all concentrations are available at new TL
      if(conc<0) 
        conc = 0;
      // distribute to density species
      for(k=0;k<8;k++)
        cv[k] += conc * ElementInSpecies[j][k];
    }
    // here, get amount of CO2 phase
    if(CO2){
      //cv_CO2 = pcs_vector[sp_pcsind[idxCO2]]->GetNodeValue(i,sp_varind[idxCO2]);// CB HS update
      //if(timeflag &! cp_vec[idxCO2]->mobil)
      //  cv_CO2 = cp_vec[idxCO2]->getProcess()->GetNodeValue(i,sp_varind[idxCO2]-1);
      //else 
        cv_CO2 = cp_vec[idxCO2]->getProcess()->GetNodeValue(i,sp_varind[idxCO2]); // all concentrations are available at new TL
    }
    else
      cv_CO2 = 0;
    // calculate and store water concentration for this node
	water_conc[i] = density::concentration_water(TT,PP,cv,cv_CO2); // CB removed from ifdef OGS_FEM_CAP
    // set the corrsponding node value of H2O transport process
    //pcs_vector[sp_pcsind[widx]]->SetNodeValue(i,1,water_conc[i]);// CB HS update
    cp_vec[widx]->getProcess()->SetNodeValue(i,1,water_conc[i]); 
  }  

  return;
}

/**************************************************************************
 ReacInt-Method:
 Task: Calculate the unit conversion factors: total or molal C
 Programing:
 //CB 01.2012 CB First implementation
 **************************************************************************/
void REACTINT::CalcUnitConversionFactors(long index, double *fl, double *fs, bool molal){
   
  if(molal){ // mol/m³-->mol/kgw
    *fl = MOLH2OPERKG / water_conc[index]; // mol/kgH2O * m³liquid/mol = m³l/kg
    *fs = (1 - node_porosity[index]) * MOLH2OPERKG / (water_conc[index] * node_porosity[index] * GetWaterSaturation(index) );
    if(*fs == 0) 
      *fs = (1 - node_porosity[index]) * MOLH2OPERKG / (water_conc[index] * node_porosity[index] * 1 );
  }
  else{      // Total moles mol/m³liquid -> mol/m³REV
    *fl = node_porosity[index] * GetWaterSaturation(index);
	if(*fl == 0.0) 
      *fl = node_porosity[index] * 1;
	*fs = 1 - node_porosity[index];
  }

}




/**************************************************************************
 ReacInt-Method:
 Task: Calculate water concentration at all nodes
 Programing:
 //CB 01.2012 CB First implementation
 **************************************************************************/
double REACTINT::GetCO2SolubilityDuan(long node){

  long j;
  double TT=298.15, PP = 1, Sal = 0, Sol = 0;
  int Naidx=-1;  
  int nc = (int)cp_vec.size();
  double unitfactor_l = 1;
  CompProperties *m_cp = NULL;

  TT = this->GetTemperature(node); // is in K
  PP = this->GetPressure(node); // is in ba
  unitfactor_l = MOLH2OPERKG / water_conc[node];
  //unitfactor_l = MOLH2OPERKG / 55335.831251734664;
  // Get water concentration index
  for(j=0;j<nc;j++){
    m_cp = cp_vec[j];    
    if(m_cp->compname.compare(SodiumSpeciesName)==0 ){
      Naidx = j;
      break;
    }
  }
  if(Naidx >= 0)
    Sal = cp_vec[Naidx]->getProcess()->GetNodeValue(node,sp_varind[Naidx]); // mol/m³l
  Sal *= unitfactor_l ;                          // mol/kgw 

  Sol = VLE::solubility_CO2(TT,PP,Sal) ; // mol/kgw
  Sol /= unitfactor_l;  
  return Sol ;
}




/**************************************************************************
FEMLib-Method: SettingInitialPorosityToElementValue()
Task: Initially, material parameters are not set for ele_value; and all the
      functions which needs porosity and permeability are adjusted to get 
      the value from ele_value, especially when PoroAndPerm is required. Thus
      preprocessing - SettingInitialPorosityAndPermToElementValue is necessary

Programing:
11.2010 AB Implementation
last modification:
**************************************************************************/
void REACTINT::SetInitialPorosityAndPermToElementValue()
{

    double n_initial, k_initial_xx;
    double k_initial_yy, k_initial_zz;
    int idx_n, idx_k, idx_k_yy, idx_k_zz, i, k, kk;

    //CElem* m_ele = NULL;
    CFEMesh* m_msh = fem_msh_vector[0];
    long n_ele = long(m_msh->ele_vector.size());

    CRFProcess *m_pcs_flow;
    m_pcs_flow = PCSGetFlow();

    // get medium property
    CMediumProperties *m_mat_mp = NULL;

    // return if porosity update is not required
    for(long n = 0; n < n_ele; n++){
        // get material group
        long n_group = m_msh->ele_vector[n]->GetPatchIndex();
        m_mat_mp = mmp_vector[n_group];
        if (m_mat_mp->porosity_model==13){
          poroupdate_flag = true; 
          break;
        }
    }
    if(poroupdate_flag==false)
      return;

    cout <<"SettingInitialPorosityAndPermToElementValue." << "\n";

    // first check if necessary component properties is defined
    k=kk=0;
    for(i=0;i<int(cp_vec.size());i++){
		  if(cp_vec[i]->transport_phase==1){ 
        if(cp_vec[i]->molar_weight>0) k++;
        if(cp_vec[i]->mineral_density>0) kk++;
		  }
	  }
    if(k==0){ 
      cout << " Warning in REACINT: molar_weight has not been defined for ANY solid phase species." << "\n";
      cout << "  No porosity update possible." << "\n";
    }
    if(kk==0){ 
      cout << " Warning in REACINT: mineral_density has not been defined for ANY solid phase species." << "\n";
      cout << "  No porosity update possible. " << "\n";
    }

    idx_n = m_pcs_flow->GetElementValueIndex("POROSITY");
    if(idx_n<0) cout << " Warning in REACINT: No POROSITY ele index found with PCS " << m_pcs_flow->getProcessType() << "." << "\n";
    idx_k = m_pcs_flow->GetElementValueIndex("PERMEABILITY");
    if(idx_k<0) cout << " Warning in REACINT: No PERMEABILITY ele index found with PCS " << m_pcs_flow->getProcessType()  << "." << "\n";
    idx_k_yy = m_pcs_flow->GetElementValueIndex("PERMEABILITY_YY");
    if(idx_k_yy<0) cout << " Warning in REACINT: No PERMEABILITY_YY ele index found with PCS " << m_pcs_flow->getProcessType()  << "." << "\n";
    idx_k_zz = m_pcs_flow->GetElementValueIndex("PERMEABILITY_ZZ");
    if(idx_k_zz<0) cout << " Warning in REACINT: No PERMEABILITY_ZZ ele index found with PCS " << m_pcs_flow->getProcessType()  << "." << "\n";

    for(long n = 0; n < n_ele; n++){

        // get material group
        long n_group = m_msh->ele_vector[n]->GetPatchIndex();
        m_mat_mp = mmp_vector[n_group];

        if((m_mat_mp->porosity_model==13)&& (aktueller_zeitschritt == 0)){

            n_initial = m_mat_mp->porosity_model_values[0];
            m_pcs_flow->SetElementValue( n, idx_n,   n_initial ); // old time level
            m_pcs_flow->SetElementValue( n, idx_n+1, n_initial ); // new time level

            //Isotropic permeability
            if (m_mat_mp->permeability_tensor_type == 0){
                
                k_initial_xx = m_mat_mp->permeability_tensor[0];
                m_pcs_flow->SetElementValue( n, idx_k,   k_initial_xx ); // old time level
                m_pcs_flow->SetElementValue( n, idx_k+1, k_initial_xx ); // new time level
            }
            //Orthogonal permeability
            else if (m_mat_mp->permeability_tensor_type == 1){

                if (m_mat_mp->GetGeoDimension() == 2){

                    k_initial_xx = m_mat_mp->permeability_tensor[0];
                    k_initial_yy = m_mat_mp->permeability_tensor[1];
                
                    m_pcs_flow->SetElementValue( n, idx_k,   k_initial_xx ); // old time level
                    m_pcs_flow->SetElementValue( n, idx_k+1, k_initial_xx ); // new time level

                    m_pcs_flow->SetElementValue( n, idx_k_yy,   k_initial_yy ); // old time level
                    m_pcs_flow->SetElementValue( n, idx_k_yy+1, k_initial_yy ); // new time level
                }
                else if (m_mat_mp->GetGeoDimension() ==  3){

                    k_initial_xx = m_mat_mp->permeability_tensor[0];
                    k_initial_yy = m_mat_mp->permeability_tensor[1];
                    k_initial_zz = m_mat_mp->permeability_tensor[2];

                    m_pcs_flow->SetElementValue( n, idx_k,   k_initial_xx ); // old time level
                    m_pcs_flow->SetElementValue( n, idx_k+1, k_initial_xx ); // new time level

                    m_pcs_flow->SetElementValue( n, idx_k_yy,   k_initial_yy ); // old time level
                    m_pcs_flow->SetElementValue( n, idx_k_yy+1, k_initial_yy ); // new time level

                    m_pcs_flow->SetElementValue( n, idx_k_zz,   k_initial_zz ); // old time level
                    m_pcs_flow->SetElementValue( n, idx_k_zz+1, k_initial_zz ); // new time level
                }
            }
        }
    }
    return;
}


/**************************************************************************
FEMLib-Method: PorosityVolumetricChemicalReactionupdate_kinetics
Task: Porosity variation owing to chemical reactions
      This fuction works when species/mineral concentration is given in [mol/m3 solid].

Programing:
11.2010 AB & CB Implementation
last modification:
**************************************************************************/
void REACTINT::PorosityVolumetricReactionUpdate()
{

  double porosity = 0.0; // n_init=0;
	double conc_old = 0.0; 
	double conc_new = 0.0; 
	double change_in_conc;
	double mineral_volume_fraction, n_previous;
  long node_idx, nnodes;
  double unitfactor, phasevolume;

	int i, j, k, l, idxC,idx_n, nImmobileComp;
	int  ns;//,no_processes,
	i = j = k = l = idxC = idx_n = nImmobileComp = 0;

  // the factors depend on units of concentration, molweight and mineraldensity 
  // I presume SI units, mol, kg, m³
  //this converts mol/m³ to m³ porosity change
  if(unitconversion) 
    unitfactor = 1;  // V=C*MW/rho*unitfactor : [-] = [mol/m³] * [kg/mol] * [m³/kg]  
  //this converts mol/L (or mol/kg) to L porosity change
  else 
    unitfactor = 1000; // V=C*MW/rho*unitfactor : [-] = [mol/kg] * [kg/mol] * [m³/kg] * [kg/cm³] 

  MeshLib::CElem* m_ele = NULL;
	CFEMesh* m_msh = fem_msh_vector[0];
  long n_ele = long(m_msh->ele_vector.size());

  cout <<"Updating porosity from Chemical reaction." << "\n";

  //vector<int>pcs_transport_comps_vector; // CB clear
  vector<int>immob_comps_vector; // CB clear
	vector<double>mineral_densities_vector;
  vector<double>mineral_molecular_weights_vector;

  CRFProcess* m_pcs_flow = NULL;
	m_pcs_flow = PCSGetFlow();
  idx_n = m_pcs_flow->GetElementValueIndex("POROSITY");

  // get the component data
  for(i=0;i<int(cp_vec.size());i++)
	{
		if(cp_vec[i]->transport_phase==1){ 
      //pcs_transport_comps_vector.push_back(sp_pcsind[i]);
      immob_comps_vector.push_back(i);
      mineral_molecular_weights_vector.push_back(cp_vec[i]->molar_weight); 
      mineral_densities_vector.push_back(cp_vec[i]->mineral_density);
			k ++;
		}
	}
	nImmobileComp = k;

  // CB: Instead of using the initial porosity always, we have to properly rescale concentrations 
  //     of species by the new porosity or solid phase volume after each time step. 
  //     see below ###
  //if(unitconversion) // mineral species are defined in mol/volume solid
  //  phasevolume = 1-n_init;
  //else               // mineral species are defined in mol/kg (or L) liquid
  //  phasevolume = n_init;

  // loop over all nodes
  nnodes = long(node_porosity.size());
  for(long n = 0; n < nnodes; n++){
      // get old porosity of node
      n_previous = node_porosity[n];
      mineral_volume_fraction = 0.0;

      if(unitconversion) // mineral species are defined in mol/volume solid
        phasevolume = 1-n_previous;
      else               // mineral species are defined in mol/kg (or L) liquid
        phasevolume = n_previous;

      // Here the porosity is updated for a component
      for(i = 0; i < nImmobileComp; i++) {
          change_in_conc = 0.0;
          l=immob_comps_vector[i]; 
          idxC = cp_vec[l]->getProcess()->GetNodeValueIndex(cp_vec[l]->getProcess()->pcs_primary_function_name[0]);
          conc_old = cp_vec[l]->getProcess()->GetNodeValue(n, idxC);// old timelevel
          conc_new = cp_vec[l]->getProcess()->GetNodeValue(n, (idxC+1));// new timelevel
          change_in_conc = conc_new - conc_old;
          //mineral_volume_fraction = change_in_conc[mol/m³solid]*molar_volume[m³min/mol]
          // the last factor depends on units of concentration, molweight and mineraldensity 
          if(mineral_densities_vector[i]>0)
            mineral_volume_fraction += change_in_conc* (mineral_molecular_weights_vector[i]/mineral_densities_vector[i]) *unitfactor;   
          else{
            cout << " Warning: mineral_density = 0 for a species in REACTINT::PorosityVolumetricReactionUpdate()" << "\n";
            cout << "  PorosityVolumetricReactionUpdate() not possible" << "\n";                  
          }
      }
      porosity = n_previous - (phasevolume * mineral_volume_fraction);
      // limit to reasonable values
      node_porosity[n] = MRange(1e-6,porosity,1.0) ;

      //some modification for Lagneau Benchmark probably required here
      //bool lagneau = false;
      if(KinReactData_vector.size()>0){
        if(KinReactData_vector[0]->lagneau ){
          if(n==1) cout << " Attention: Lagneau-Benchmark with modified poro update" << "\n";
          int sp=0;
          for(sp=0;sp<int(cp_vec.size());sp++){
            if (cp_vec[sp]->mobil==0){
              //l = sp_pcsind[sp]; 
              idxC = this->sp_varind[sp];
              break;
            }
          }
          phasevolume = n_previous;
          conc_old = cp_vec[sp]->getProcess()->GetNodeValue(n, idxC-1);// old timelevel
          conc_new = cp_vec[sp]->getProcess()->GetNodeValue(n, (idxC));// new timelevel
          change_in_conc = conc_new - conc_old;
          mineral_volume_fraction = change_in_conc * (cp_vec[sp]->molar_weight/cp_vec[sp]->mineral_density) * unitfactor;   
          node_porosity[n] = n_previous - (phasevolume * mineral_volume_fraction);
        }
      }
      // ### 
      // CB: Here, we should include a rescaling of concentrations in the phases,
      //     at least for solid phase concentrations
      ns = cp_vec.size();
      for(i = 0; i < ns; i++) {
          //l    = sp_pcsind[i];
          idxC = sp_varind[i];
          conc_new = cp_vec[i]->getProcess()->GetNodeValue(n, idxC);// new timelevel
          if(unitconversion){ // liquid and solid phase concentrations mol/m³phase
            if(cp_vec[i]->transport_phase==1)
              conc_new *= (1-n_previous)/(1-node_porosity[n]);
            // for dissolved species, is this handled by dn/dt term in transport equation??
            // then comment out next two lines 
            else if(cp_vec[i]->transport_phase==0)
              conc_new *= n_previous/node_porosity[n];
          }
          else // in this case, alls species are related to liquid phase mol/kgH2O
            conc_new *= n_previous/node_porosity[n];
          //pcs_vector[l]->SetNodeValue(n,idxC,conc_new);
          cp_vec[i]->getProcess()->SetNodeValue(n,idxC,conc_new);
      }

  } // loop over nodes

  // now set the new element porosities from node_porosity for mass transport and other processes
  // loop over all elements
  for(long n = 0; n < n_ele; n++){
    m_ele = m_msh->ele_vector[n];
    long nn = m_ele->GetNodesNumber(false);
    porosity  = 0;
    for(j = 0; j < nn; j++){
        node_idx = m_ele->GetNodeIndex(j);
        porosity  += node_porosity[node_idx];
    }
    //arithmetic average ... inverse distance?
    porosity  /= (double) nn;
    //set the element value
    m_pcs_flow->SetElementValue( n, idx_n+1, porosity );
  }
  
  //pcs_transport_comps_vector.clear();
  immob_comps_vector.clear();
  mineral_densities_vector.clear();
  mineral_molecular_weights_vector.clear();

  
   CTimeDiscretization *m_tim = NULL;
   m_tim = time_vector[0];
   int steps = int(m_tim->time_step_vector.size());
   bool plot = false;   
   if(aktueller_zeitschritt == 1 || aktueller_zeitschritt % 10 == 0 || steps == int(aktueller_zeitschritt))
      plot = true;

  // debug output
  if(KinReactData_vector.size()>0 && plot == true){
    if(KinReactData_vector[0]->debugoutflag){
      ofstream aus;
      string file = FileName + "_node_porosities.dump";
      aus.setf(ios::scientific,ios::floatfield);
      aus.precision(12);
      aus.open(file.c_str()); 
      for(i=0;i<nnodes;i++) aus << node_porosity[i] << "\n";
      aus.close();
    }
  }

  return ;
}

/**************************************************************************
FEMLib-Method: 
Task: not used at the moment, functionality is currently implemented in
   CKinReactData::CopyConcentrations
Programing:
2.2011 CB Implementation
last modification:
**************************************************************************/
void REACTINT::CopySymmetricConcentrationsInRadialModel()
{

 // long i,j, nnodes;
 // int nc = (int)cp_vec.size();

	//CFEMesh* m_msh = fem_msh_vector[0];
 // nnodes = long(m_msh->nod_vector.size());

 // CKinReactData *m_krd = NULL;
 // if(KinReactData_vector.size()>0)
	//  m_krd = KinReactData_vector[0];
 // if(m_krd == NULL){
 // 	// no KinReactData specified in *.krc file
	//  cout << "No CKinReactData available in CopySymmetricConcentrationsInRadialModel" << "\n";
 //   return;
 // }

 // double vecnodes[2];

 // vecnodes[0]=1;
 // vecnodes[1]=2;

 // for(i=0;j<nnodes;i++)
 //   if( i != vecnodes[i])
 //     for(j=0;j<nc;j++)
 //       pcs_vector[sp_pcsind[j]]->SetNodeValue( vecnodes[i], sp_varind[j], pcs_vector[sp_pcsind[j]]->GetNodeValue(i,sp_varind[j]) );

  return ;
}

/**************************************************************************
FEMLib-Method: 
Task: not used at the moment, functionality is currently implemented in
   CKinReactData::CopyConcentrations
Programing:
2.2011 CB Implementation
last modification:
**************************************************************************/
void REACTINT::DumpSolidSpeciesMoles(void){

   long i, j, nele, nnodes;
   const double * coord;
    MeshLib::CNode* m_nod = NULL;
   CFEMesh* m_msh = fem_msh_vector[0];            //SB: ToDo hart gesetzt
   ofstream aus;
   string filename = FileName + "_mineral_moles.tec";
   string eleType;

   if (dump_min_moles==false)
     return;

   bool plot = false;

   CTimeDiscretization *m_tim = NULL;
   m_tim = time_vector[0];
   int steps = int(m_tim->time_step_vector.size());

   nnodes = (long) m_msh->nod_vector.size();
   nele = (long) m_msh->ele_vector.size();

  CMediumProperties *m_mat_mp = NULL;
  
  for(long n = 0; n < nele; n++){
      long n_group = m_msh->ele_vector[n]->GetPatchIndex();
      m_mat_mp = mmp_vector[n_group];
      if (m_mat_mp->porosity_model==13){
        plot = true; 
        break;
      }
  }
  
  if(plot == true){
    if(aktueller_zeitschritt == 1 || aktueller_zeitschritt % ssp_outstep == 0 || steps == int(aktueller_zeitschritt))
      plot = true;
    else
      plot = false;
  }
  if(plot==false)
    return;

   if (m_msh->getNumberOfLines () > 0)
      eleType = "QUADRILATERAL";
   if (m_msh->getNumberOfQuads () > 0)
      eleType = "QUADRILATERAL";
   if (m_msh->getNumberOfHexs () > 0)
      eleType = "BRICK";
   if (m_msh->getNumberOfTris () > 0)
      eleType = "QUADRILATERAL";
   if (m_msh->getNumberOfTets () > 0)
      eleType = "TETRAHEDRON";
   if (m_msh->getNumberOfPrisms () > 0)
      eleType = "BRICK";

  if (aktueller_zeitschritt == 1)
     aus.open(filename.c_str());
  else
     aus.open(filename.c_str(), ios::app);

    //Header
    aus << "VARIABLES = " << "\"x\"" << " " << "\"y\"" << " " << "\"z\"" << " " ;
    for (i = 0; i < long(cp_vec.size()); i++)
      if(cp_vec[i]->transport_phase==1) 
        aus << "\"" << cp_vec[i]->compname << "\" ";
    aus << "\"" << "nodeporosity" << "\" " << "\n";
    aus << "ZONE T=" << "\"aktueller_zeitschritt=" << aktueller_zeitschritt
       << "\"";
    aus << ", N=" << nnodes << ", E=" << nele << " F=FEPOINT, ET="
       << eleType << "\n";
    //nodedata
    aus.setf(ios::scientific,ios::floatfield);
    aus.precision(12);
    for (i = 0; i < nnodes; i++)
    {
       m_nod = m_msh->nod_vector[i];
       coord = m_nod->getData();
       aus << coord[0] << " " << coord[1] << " " << coord[2] << " ";
       for (j = 0; j < long(cp_vec.size()); j++){
         if(cp_vec[j]->transport_phase==1) // plot total mol / m³ aquifer = c*(1-n)
           aus << cp_vec[j]->getProcess()->GetNodeValue(i, sp_varind[j]) * (1-node_porosity[i]) << " "; 
       }
       aus << node_porosity[i] << "\n";
    }
    //eledata
    for (i = 0l; i < nele; i++)
       m_msh->ele_vector[i]->WriteIndex_TEC(aus);

  aus.close();

  return;
}

/**************************************************************************
FEMLib-Method: 
Task: 

Programing:
2.2011 CB Implementation
last modification:
**************************************************************************/
void REACTINT::DumpAllVariables()
{

  int idx;
  long nnodes; 
  string file;
  ofstream _dump;
  bool plot = false;
  
  if(dump_all_pcs==false)
   return;

  CTimeDiscretization *m_tim = NULL;
  m_tim = time_vector[0];
  int step = int(m_tim->time_step_vector.size());

  if(aktueller_zeitschritt % this->pcs_outstep == 0 || step == int(aktueller_zeitschritt))
    plot = true;

  if(plot == false)
    return;

  nnodes = (long) fem_msh_vector[0]->nod_vector.size();

  for(long i=0;i<long(pcs_vector.size());i++){
    for(int j = 0; j <  pcs_vector[i]->pcs_number_of_primary_nvals; j++){
      // each PV is stored twice: old TL, new TL
      idx = pcs_vector[i]->GetNodeValueIndex(pcs_vector[i]->nod_val_name_vector[j*2])+1;
      // Restartfiles ; each PV name is stored twice: old TL, new TL
      file = FileName + "_" + pcs_vector[i]->nod_val_name_vector[j*2] + ".dump";
      _dump.setf(ios::scientific,ios::floatfield);
      _dump.precision(12);
      _dump.open(file.c_str());
      // header
      _dump << "#0#0#0#1#100000#0#4.2.13 #########################################" << "\n";
      _dump << "1 1 4" << "\n" << "1 1" << "\n";
      _dump << pcs_vector[i]->nod_val_name_vector[j*2] << ", [-]" << "\n";
      // data
      for(long k = 0; k < nnodes; k++)
        _dump << k << " " << pcs_vector[i]->GetNodeValue(k, idx) << " " << "\n" << flush;
      _dump.close();
    }
  }

}


/**************************************************************************
Task: Kill simulation when enough NAPL is infiltrated

Programing:
   3/2013   CB   Implementation                                          */
/**************************************************************************/
void REACTINT::DumpMassIntegrals()
{

  size_t  k;
  long j, i;
  long nod;
  double  volume ;//mass,
  long nelems, nnodes; 
  string file;
  ofstream _dump;
 
  CFEMesh* m_msh = NULL;

 std::vector<double> massintegrals;
 for(k=0;k<cp_vec.size();k++)
   massintegrals.push_back(0);


 // m_tim = time_vector[0];
  m_msh = fem_msh_vector[0]; //SB: ToDo hart gesetzt
 // m_mat_mp = mmp_vector[0];
  nelems = (long) m_msh->ele_vector.size();
  //nnodes = (long) m_msh->nod_vector.size();


// Tecfile
file = FileName + "_masses_divide_by_phasevol.tec";
 _dump.setf(ios::scientific,ios::floatfield);
_dump.precision(12);

if(aktueller_zeitschritt == 1){ // new file
  _dump.open(file.c_str());

  // header
  _dump << "VARIABLES  = \"TIME\",\"" ;
  for(k = 0; k < cp_vec.size(); k++) 
    _dump << cp_vec[k]->compname << "\",\"";
  _dump << "\n"<< flush;

  // integrate initial masses
  for(i = 0; i < nelems; i++){
    volume = m_msh->ele_vector[i]->GetVolume();
    nnodes = m_msh->ele_vector[i]->GetNodesNumber(false);
    //mass = 0;
    for(j = 0; j < nnodes; j++){
      nod = m_msh->ele_vector[i]->GetNodeIndex(j);
      for(k = 0; k < cp_vec.size(); k++){
        massintegrals[k] += cp_vec[k]->getProcess()->GetNodeValue(nod, sp_varind[k]-1) / nnodes * volume;
      }
    }
  }
  // output of initial masses
  _dump << 0.0  ;
  for(k = 0; k < cp_vec.size(); k++) 
    _dump << " " << massintegrals[k] ;
  _dump << "\n"<< flush;

}
else { // running model, append
  _dump.open(file.c_str(), ios::app);
}

// integrate and output for current time step
for(k = 0; k < cp_vec.size(); k++) 
  massintegrals[k] = 0;
// integrate
for(i = 0; i < nelems; i++){
  volume = m_msh->ele_vector[i]->GetVolume();
  nnodes = m_msh->ele_vector[i]->GetNodesNumber(false);
  //mass = 0;
  for(j = 0; j < nnodes; j++){
    nod = m_msh->ele_vector[i]->GetNodeIndex(j);
    for(k = 0; k < cp_vec.size(); k++){
      massintegrals[k] += cp_vec[k]->getProcess()->GetNodeValue(nod, sp_varind[k]) / nnodes * volume;
    }
  }
}

// output
_dump << aktuelle_zeit  ;
for(k = 0; k < cp_vec.size(); k++) 
  _dump << " " << massintegrals[k] ;
_dump << "\n"<< flush;

// close
_dump.close();

massintegrals.clear();

}


void REACTINT::PrintGCCdata(void){

  return;

  string file;
  ofstream _dump;

  file = FileName + "_GCC_data.txt";
  _dump.setf(ios::scientific, ios::floatfield);
  _dump.precision(12);
  _dump.open(file.c_str());

  // integrate initial masses
  
  double P0 = 0.1;
  double T0 = 273.15;
  double T, P;

  for (int i = 0; i < 100; i++){
    T = T0 + (double)i * 1 ;
    P = P0 ;
    _dump << T << " " << P0 << " " << IF97::cp(T, P) << endl; // kJ/kg/K 
  

  
  }
  _dump.close();


};





/**************************************************************************
FEMLib-Method:
Task: Permeability will be calculated using either Kozeny-Carman or 
      Verma-Pruess formulation. And obviously it writes the new permeability
      on ele_value. It will be called in post coupling process thereby
      it will be executed once.

Programing:
2010.11 AB Implementation
last modification:
**************************************************************************/
void REACTINT::PermeabilityPorosityUpdate()
{

    //static double tensor[1];
    int idx_k = -1, idx_n = -1, idx_k_yy = -1, idx_k_zz = -1;
    double k_init = 0.0, n_init = 0.0;
    double n_t = 0.0, k_t = 0.0, k_t_new = 0.0;
    double k_t_yy = 0.0, k_t_zz = 0.0, k_t_new_yy = 0.0;
    double k_init_yy =0.0, k_init_zz = 0.0, k_t_new_zz = 0.0;
    long n;

    //CElem* m_ele = NULL;
    CFEMesh* m_msh = fem_msh_vector[0];
    long n_ele = long(m_msh->ele_vector.size());
    CRFProcess *m_pcs_flow;
    m_pcs_flow = PCSGetFlow();

    // get variable index
    idx_n = m_pcs_flow->GetElementValueIndex("POROSITY");
    idx_k = m_pcs_flow->GetElementValueIndex("PERMEABILITY");
    idx_k_yy = m_pcs_flow->GetElementValueIndex("PERMEABILITY_YY");
    idx_k_zz = m_pcs_flow->GetElementValueIndex("PERMEABILITY_ZZ");

    // get medium property
    CMediumProperties *m_mat_mp = NULL;

    for(n = 0; n < n_ele; n++){

        // get material group
        long n_group = m_msh->ele_vector[n]->GetPatchIndex();
        m_mat_mp = mmp_vector[n_group];

        // Permeability updating is not required
        if (((m_mat_mp->permeability_tensor_type!= 0) || (m_mat_mp->permeability_tensor_type!= 1))&& m_mat_mp->permeability_model != 8 )
            return;
            
        if(n==0) cout <<"Updating Permeability:" << "\n";

            // get current & initial porosity: n_t, n_init.
            n_init = m_mat_mp->porosity_model_values[0];
            n_t = m_pcs_flow->GetElementValue( n, idx_n+1 ); // new time level, recently set in PorosityVolumetricReactionUpdate
            //double test1 = m_pcs_flow->GetElementValue( n, idx_n );
            //if(n==1)
            //  cout << test1 << " " << n_init << " " << n_t << "\n";

            
            if (m_mat_mp->permeability_tensor_type == 0){

                k_init = m_mat_mp->permeability_tensor[0];
                
                // get and save old permeability
                k_t = m_pcs_flow->GetElementValue( n, idx_k + 1 );
                m_pcs_flow->SetElementValue( n, idx_k, k_t);
                //double test = m_pcs_flow->GetElementValue( n, idx_k);
            

                /*/ ***************************************************** /*/
                // If Kozeny-Carman formulation is choosen
                if(m_mat_mp->permeability_porosity_updating_type == 0){

                    // Calculate for the new perm. k_t_new
                    k_t_new = m_mat_mp->KozenyCarmanNew( k_init,n_init,n_t );
                }
                
                // If Verma-Pruess formulation is choosen
                else if(m_mat_mp->permeability_porosity_updating_type == 1){

                    // Calculate for the new perm. k_t_new
                    k_t_new = m_mat_mp->VermaPruess( k_init,n_init,n_t );
                }
                /*/ ***************************************************** /*/

                // save new permeability: in index+1
                m_pcs_flow->SetElementValue( n, idx_k+1, k_t_new	);
                //double test_2 = m_pcs_flow->GetElementValue( n, idx_k+1);
            
                /*/ ***************************************************** /*/
            }

            if (m_mat_mp->permeability_tensor_type == 1){

              if (m_mat_mp->GetGeoDimension()==2){

                    k_init = m_mat_mp->permeability_tensor[0];
                    k_init_yy = m_mat_mp->permeability_tensor[1];

                    // get and save old permeability
                    k_t = m_pcs_flow->GetElementValue( n, idx_k + 1 );
                    k_t_yy = m_pcs_flow->GetElementValue( n, idx_k_yy + 1 );
                    m_pcs_flow->SetElementValue( n, idx_k, k_t);
                    m_pcs_flow->SetElementValue( n, idx_k_yy, k_t_yy);


                    /*/ ***************************************************** /*/
                    // If Kozeny-Carman formulation is choosen
                    if(m_mat_mp->permeability_porosity_updating_type == 0){

                        k_t_new = m_mat_mp->KozenyCarmanNew( k_init,n_init,n_t );
                        k_t_new_yy = m_mat_mp->KozenyCarmanNew( k_init_yy,n_init,n_t );
                    }
                    
                    // If Verma-Pruess formulation is choosen
                    else if(m_mat_mp->permeability_porosity_updating_type == 1){

                        k_t_new = m_mat_mp->VermaPruess( k_init,n_init,n_t );
                        k_t_new_yy = m_mat_mp->VermaPruess( k_init_yy,n_init,n_t );
                    }
                    /*/ ***************************************************** /*/

                    // save new permeability: in index+1
                    m_pcs_flow->SetElementValue( n, idx_k+1, k_t_new	);
                    m_pcs_flow->SetElementValue( n, idx_k_yy+1, k_t_new_yy	);
                
                    /*/ ***************************************************** /*/
                }

              if (m_mat_mp->GetGeoDimension() ==3){

                    k_init = m_mat_mp->permeability_tensor[0];
                    k_init_yy = m_mat_mp->permeability_tensor[1];
                    k_init_zz = m_mat_mp->permeability_tensor[2];

                    // get and save old permeability
                    k_t = m_pcs_flow->GetElementValue( n, idx_k + 1 );
                    k_t_yy = m_pcs_flow->GetElementValue( n, idx_k_yy + 1 );
                    k_t_zz = m_pcs_flow->GetElementValue( n, idx_k_zz + 1 );
                    m_pcs_flow->SetElementValue( n, idx_k, k_t);
                    m_pcs_flow->SetElementValue( n, idx_k_yy, k_t_yy);
                    m_pcs_flow->SetElementValue( n, idx_k_zz, k_t_zz);


                    /*/ ***************************************************** /*/
                    // If Kozeny-Carman formulation is choosen
                    if(m_mat_mp->permeability_porosity_updating_type == 0){

                        k_t_new = m_mat_mp->KozenyCarmanNew( k_init,n_init,n_t );
                        k_t_new_yy = m_mat_mp->KozenyCarmanNew( k_init_yy,n_init,n_t );
                        k_t_new_zz = m_mat_mp->KozenyCarmanNew( k_init_zz,n_init,n_t );
                    }
                    
                    // If Verma-Pruess formulation is choosen
                    else if(m_mat_mp->permeability_porosity_updating_type == 1){

                        k_t_new = m_mat_mp->VermaPruess( k_init,n_init,n_t );
                        k_t_new_yy = m_mat_mp->VermaPruess( k_init_yy,n_init,n_t );
                        k_t_new_zz = m_mat_mp->VermaPruess( k_init_zz,n_init,n_t );
                    }
                    /*/ ***************************************************** /*/

                    // save new permeability: in index+1
                    m_pcs_flow->SetElementValue( n, idx_k+1, k_t_new	);
                    m_pcs_flow->SetElementValue( n, idx_k_yy+1, k_t_new_yy	);
                    m_pcs_flow->SetElementValue( n, idx_k_zz+1, k_t_new_zz	);
                
                    /*/ ***************************************************** /*/
                }
            }
    }
    return ;
}





/**************************************************************************
 ReacInt-Method:
 Task: Convert name (chemical formula) of species to indices
    for water concentration calculation
    0-Li, 1-Na, 2-K, 3-Mg, 4-Ca, 5-Cl, 6-SO4, 7-CO3, 8-H, 9-O
 Programing:
 //DL 01.2011 DL First implementation
 **************************************************************************/

std::vector<int> REACTINT::formula2index(std::string formula){
	const string Chemical_Element[10]={"Li", "Na", "K", "Mg", "Ca", "Cl", "S", "C", "H", "O"};
	int i,j,n,symb_asc,is_bracket;
	vector<int> id_iz,type_x;
	vector<int> bracket_ia, bracket_ib, bracket_iz, bracket_ia0, bracket_ib0; // CB clear
	vector<int> element_ia, element_ib, element_iz;
	vector<int>  number_ia,  number_ib,  number_iz;
	vector<string> element_name;

	id_iz.clear();
	type_x.clear();
	n=int(formula.size());
	
	element_ia.clear();
	element_ib.clear();
	element_iz.clear();
	element_name.clear();
	
	number_ia.clear();
	number_ib.clear();
	number_iz.clear();

	bracket_ia.clear(); //to store begin position of bracket
	bracket_ib.clear(); //to store end position of bracket
	bracket_iz.clear();
	bracket_ia0.clear();
	bracket_ib0.clear();

	// 0-A-Z, 1-a-z, 2-0-9, 3-[(, 4-]), 5-others
	for(i=0;i<n;i++){
		symb_asc=int(formula[i]);

		if(symb_asc>=65 && symb_asc<=90)  //A-Z
			type_x.push_back(0);
		else if(symb_asc>=97 && symb_asc<=122) //a-z
			type_x.push_back(1);
		else if(symb_asc>=48 && symb_asc<=57)  //0-9
			type_x.push_back(2);
		else if(symb_asc==40 || symb_asc==91)  // [(
			type_x.push_back(3);
		else if(symb_asc==41 || symb_asc==93)  // ])
			type_x.push_back(4);
		else 
			type_x.push_back(5);
	}
	type_x.push_back(-1); //set string end mark


	// 0-A-Z, 1-a-z, 2-0-9, 3-[(, 4-]), 5-others
	//search for chemical elements
	for(i=0;i<n;i++){
		if(type_x[i]==0)
			if(type_x[i+1]==1){
				element_name.push_back(formula.substr(i,2));
				element_ia.push_back(i);
				element_ib.push_back(i+1);
			}
			else {
				element_name.push_back(formula.substr(i,1));
				element_ia.push_back(i);
				element_ib.push_back(i);
			}
		else if(type_x[i]==2){
			if(i>0 && type_x[i-1]!=2)
				number_ia.push_back(i);
			if(i>0 && type_x[i]==2 && type_x[i+1]!=2)
				number_ib.push_back(i);
			}
		else if(type_x[i]==3)
			bracket_ia0.push_back(i);
		else if(type_x[i]==4)
			bracket_ib0.push_back(i);
	}

	//to check the bracket symb is correct or not, e.g.  (( ... ))
	for(i=0;i<(int)bracket_ia0.size();i++){
		is_bracket=0;
		for(j=0;j<(int)bracket_ib0.size();j++){
			if(i==(int)bracket_ia0.size()-1){
				if(bracket_ib0[j]>bracket_ia0[i])
					is_bracket=1;
			}
			else{
				if(bracket_ib0[j]>bracket_ia0[i] && bracket_ib0[j]<bracket_ia0[i+1])
					is_bracket=1;	
			}
		}
		if(is_bracket==1)
			bracket_ia.push_back(bracket_ia0[i]);
	}
	for(i=0;i<(int)bracket_ib0.size();i++){
		is_bracket=0;
		for(j=0;j<(int)bracket_ia0.size();j++){
			if(i==0){
				if(bracket_ia0[j]<bracket_ib0[i])
					is_bracket=1;
			}
			else{
				if(bracket_ia0[j]<bracket_ib0[i] && bracket_ia0[j]>bracket_ib0[i-1])
					is_bracket=1;	
			}
		}
		if(is_bracket==1)
			bracket_ib.push_back(bracket_ib0[i]);
	}



	for(i=0;i<(int)number_ia.size();i++)
		number_iz.push_back(atoi(formula.substr(number_ia[i],number_ib[i]-number_ia[i]+1).c_str()));

	for(i=0;i<(int)element_name.size();i++)
		if(type_x[element_ib[i]+1]==2)
			for(j=0;j<(int)number_iz.size();j++){
				if(element_ib[i]+1==number_ia[j])
					element_iz.push_back(number_iz[j]);
			}
		else
			element_iz.push_back(1);

	for(i=0;i<(int)bracket_ia.size();i++)
		if(type_x[bracket_ib[i]+1]==2)
			for(j=0;j<(int)number_iz.size();j++){
				if(bracket_ib[i]+1==number_ia[j])
					bracket_iz.push_back(number_iz[j]);
			}
		else
			bracket_iz.push_back(1);

	for(i=0;i<(int)bracket_ia.size();i++)
		for(j=0;j<(int)element_ia.size();j++)
			if(element_ia[j]>bracket_ia[i] && element_ib[j]<bracket_ib[i])
				element_iz[j] *= bracket_iz[i];

	for(i=0;i<10;i++){
		id_iz.push_back(0);
		for(j=0;j<int(element_ia.size());j++)
			if(Chemical_Element[i]==element_name[j])
				id_iz[i] += element_iz[j];
	}

	//for(i=0;i<(int)id_iz.size();i++)
	//	cout << " " << i << " " << id_iz[i] << "\n";
	return id_iz;
}


//split string line to pieces, and store in a vector
vector<std::string> REACTINT::string2vector(std::string line){
	stringstream in;
	std::string sp;
	vector<std::string> pies;
	pies.clear();
	in.str(line);
	while(1){
		if(in.eof()) break;
		sp="";
		in >> sp;	
		if(sp!="") pies.push_back(sp);
	}
	return pies;
}


/**************************************************************************/
/* Return the volume fraction of a particular phase at a node             */
/* 0 pore space, 1 solid phase, 2 bio phase                               */
/* DS-TBC                                                                 */
/* 09/2009     CB         Introduced new C++ concept, Data structures     */
/**************************************************************************/
double GetNodePhaseVolume(long node, double theta, int phase){

  CMediumProperties *m_mat_mp = NULL;
  MeshLib::CNode* m_nod = NULL;
  MeshLib::CElem* m_ele = NULL;
  //OK411 CRFProcess *m_pcs = NULL;
  CFEMesh* m_msh = fem_msh_vector[0]; //SB: ToDo hart gesetzt

  long idx=0, i, el, elem ,group; //OK411
  const double * coord;
  double distance, weight, sum_w ;
  const double*  grav_c;
  double vol = 0, poro = 0;

  // get Indices for phase 1 or 2, only if heterogeneous porosity model = 11, i.e. vol_mat_model = vol_bio_model = 2
  group = 0; //SB todo group = m_ele->GetPatchIndex(); Todo CB
  m_mat_mp = mmp_vector[group];
  if(m_mat_mp->vol_bio_model == 2 && m_mat_mp->vol_mat_model == 2 ){
    CFEMesh* _mmesh = m_mat_mp->getMesh();
    switch(phase) {
      case 1: //solid phase
        // Get VOL_MAT index
        for(idx=0;idx<(int)_mmesh->mat_names_vector.size();idx++){
            if(_mmesh->mat_names_vector[idx].compare("VOL_MAT")==0) break;        
        }
        break;
      case 2: //bio phase
        // Get VOL_BIO index
        for(idx=0;idx<(int)_mmesh->mat_names_vector.size();idx++){
          if(_mmesh->mat_names_vector[idx].compare("VOL_BIO")==0) break;
        }
        break;
      default: 
        break;
    }
  }

  // initialize data structures
  sum_w = 0;

  // Get node coordinates
  m_nod = m_msh->nod_vector[node];
  coord = m_nod->getData();

  for(el=0;el<(int)m_nod->getConnectedElementIDs().size();el++){
    // initialize for each connected element
    distance = weight = poro = 0; 
    // Get the connected element
    elem = m_nod->getConnectedElementIDs()[el]; // element index
    m_ele = m_msh->ele_vector[elem];  
    //get the phase volume of current element elem
    group = m_ele->GetPatchIndex(); 
    m_mat_mp = mmp_vector[group];
    switch(phase) {
      case 0: //pore space
        poro = m_mat_mp->Porosity(elem, theta); // CB Now provides also heterogeneous porosity, model 11
        break;
      case 1: //solid phase
        if(m_mat_mp->vol_mat_model == 1) // homogeneous
   	      poro = m_mat_mp->vol_mat;
        else if(m_mat_mp->vol_mat_model == 2)  // CB heterogeneous
          poro = m_ele->mat_vector(idx) ;
        else
          cout << "Warning! No valid VOL_MAT model in CKinReact::GetPhaseVolumeAtNode, vol_mat_model =" << m_mat_mp->vol_mat_model << "\n";
        break;
      case 2: //bio phase
        if(m_mat_mp->vol_bio_model == 1) // homogeneous
          poro = m_mat_mp->vol_bio;
        else if(m_mat_mp->vol_bio_model == 2)  // CB heterogeneous
          poro = m_ele->mat_vector(idx) ;
        else
          cout << "Warning! No valid VOL_BIO model in CKinReact::GetPhaseVolumeAtNode, vol_bio_model =" << m_mat_mp->vol_bio_model << "\n";
        break;
      case 3: // NAPL phase (refers to REV)
	       poro = 1.0;
        break;
      default:
        cout << "Error in CKinReact::GetPhaseVolumeAtNode: no valid phase" << "\n";
        break;
    }
    // calculate distance node <-> element center of gravity  
    grav_c = m_ele->GetGravityCenter();
    for(i=0;i<3;i++)
      distance += pow((coord[i]-grav_c[i]),2);
    // linear inverse distance weight = 1/(distance)
    distance = sqrt(distance);// for quadratic interpolation uncomment this line
    weight=(1/distance);
    sum_w += weight;
    // add the weighted phase volume
    vol += poro*weight;
  } // loop over connected elements

  // normalize weighted sum by sum_of_weights sum_w
  vol *= 1/sum_w;
  if(sum_w==0) 
	  vol = 0.99;
  return vol;
}


/*****************************************************************************************/
/* Calculate the reference volume of a phase at a node                                   */
/* DS-TBC                                                                                */
/* 02/2006     SB         Introduced new C++ concept, Data structures                    */
/* 08/2008     DS         Consider saturation of water phase in case of multiphase flow  */
/* 09/2009     CB         Heterogeneous Porosities update                                */
/*****************************************************************************************/
double REACTINT::GetWaterSaturation( long index ){

  double saturation=1;
  int timelevel = 1;
  CRFProcess *m_pcs=NULL;
  int idx;
  // special case for initial computation
  // saturation is not yet set for new time level
  if(aktueller_zeitschritt==0) timelevel = 0; 
  
  switch(flowtype){
    case 0:
      break;
    case 1: //Groundwater Flow //Liquid_Flow
      m_pcs = PCSGetFlow(); 
      if (m_pcs->getProcessType() == FiniteElement::LIQUID_FLOW && m_pcs->napl_dissolution)
      {
        // Sat of water phase
        idx = m_pcs->GetNodeValueIndex("SATURATION1");
        saturation = m_pcs->GetNodeValue(index, idx);
      }
      break;
    case 66: //Overland Flow
      break;
    case 5: //Air Flow
      break;
    case 11: //Componental Flow
      break;
    case 1212: //Multiphase Flow
      m_pcs = PCSGetFlow();
      idx = m_pcs->GetNodeValueIndex("SATURATION1"); // Sat of water phase
      saturation = m_pcs->GetNodeValue(index, idx + timelevel); 
      break;
    case 12: // Two_phase_Flow
      m_pcs = PCSGetFlow();
      idx = m_pcs->GetNodeValueIndex("SATURATION1"); // Sat of water phase
      saturation = m_pcs->GetNodeValue(index, idx + timelevel); 
      break;
    case 1313: // PS_GLOBAL
      m_pcs = PCSGetFlow();
      idx = m_pcs->GetNodeValueIndex("SATURATION1"); // Sat of water phase
      saturation = m_pcs->GetNodeValue(index, idx + timelevel); 
      break;    
    case 22: // Richards flow
      m_pcs = PCSGetFlow();
      idx = m_pcs->GetNodeValueIndex("SATURATION1"); // Sat of water phase
      saturation = m_pcs->GetNodeValue(index, idx + timelevel); 
      break;
    default:
      break;
  }
  return saturation;
}
/**************************************************************************
FEMLib-Method:
Task:
   Dynamische Fluessigkeits-Viskositaet nach Yaws et al. (1976)
   als Funktion von Temperatur
   in Reid et al. (1988), S. 441/455
   Eqn.(3): ln(my) = A + B/T + CT + DT^2
Programing:
08/2004 OK MFP implementation
           based on CalcFluidViscosityMethod8 by OK (06/2001)
last modification:
**************************************************************************/
double REACTINT::LiquidViscosity_Yaws_1976(double T)
{
   double ln_my,my;
   double A,B,C,D;

   A = -2.471e+01;
   B =  4.209e+03;
   C =  4.527e-02;
   D = -3.376e-5;

   //A = -11.6225;     //values for water viscosity after YAWS et al. (1976)
   //B =  1.9490e+03;
   //C =  2.1641e-02;
   //D = -1.5990e-05;

   ln_my = A + B/T + C*T + D*T*T;
   //my = pow(10.0, ln_my);
   my = exp(ln_my);             // CB exp --> pow(10)     /* in cP */
   my *= 1.e-3;                               /* in Pa s */
   return my;
}


double REACTINT::LiquidDensity_Busch(double T)
{
   double density;
   double rho_0, T_0, drho_dT;

   rho_0 = 998.203;
   T_0 = 293;
   drho_dT = -2.5e-4;

   // rho(T) = rho_0*(1+beta_T*(T-T_0))
   density = rho_0 * (1.0 + drho_dT * (max(T,0.0) - T_0) );

   return density;
}



// External Functions
//----for pressure update at fixed phases conditions----
// 2012.06.13 DL
void VLE_CalcNewPressure(double T, double &P, double &V_gas, double &V_liquid, double CO2_gas_old, double CO2_gas_new, double CO2_liquid_old, double CO2_liquid_new, double H2O_liquid_old, double H2O_liquid_new, double &rho){
	Phase_Properties2 gas,liquid,solid;
	int f=1;
	
	gas.CO2=CO2_gas_old;
	gas.H2O=0.0;
	gas.NaCl=0.0;

	liquid.CO2=CO2_liquid_old;
	liquid.H2O=H2O_liquid_old;
	liquid.NaCl=0.0;

	solid.CO2=0.0;
	solid.H2O=0.0;
	solid.NaCl=0.0;

	VLE_isobaric_fixphase(T,P,gas,liquid,solid,f);

	gas.CO2=CO2_gas_new;
	liquid.CO2=CO2_liquid_new;
	liquid.H2O=H2O_liquid_new;

	VLE_isochoric_fixphase(T,P,gas,liquid,solid,f);

	V_gas=gas.volume;
	V_liquid=liquid.volume;
    rho = VLE::density_CO2(T,P); //to use mixture fluid EoS
}

void VLE_isobaric_fixphase(double T, double P, Phase_Properties2 &vapor, Phase_Properties2 &liquid, Phase_Properties2 & /*solid*/, int /*f*/){
	double wH2O=0, mCO2=0, mNaCl=0;
	wH2O=liquid.H2O/55.51;
	mCO2=liquid.CO2/wH2O;
	mNaCl=liquid.NaCl/wH2O;

	liquid.temperature=T;
	liquid.pressure=P;
    liquid.density=density::CO2brine(T,P,mNaCl,mCO2);
	liquid.viscosity=-1.0;
	liquid.mass=liquid.H2O*18.01528+liquid.CO2*44.009+liquid.NaCl*58.443;
	liquid.volume=liquid.mass/liquid.density;

	vapor.temperature=T;
	vapor.pressure=P;
	vapor.density=VLE::density_CO2(T,P); //to use mixture fluid EoS
	vapor.viscosity=-1.0;
	vapor.mass=vapor.CO2*44.009+vapor.H2O*18.01528;
	vapor.volume=vapor.mass/vapor.density;

}

void VLE_isochoric_fixphase(double T, double &P, Phase_Properties2 &vapor, Phase_Properties2 &liquid, Phase_Properties2 &solid, int f){
	int i=0,i_max=31;
	double V=0, V0=0, P1=0, P2=0;
	double err=1.0e-8;

	V0 =vapor.volume+liquid.volume;
	VLE_isobaric_fixphase(T, P, vapor, liquid, solid, f);
	V=vapor.volume+liquid.volume;

	if(abs(V-V0)<err)
		return;
	else if(V>V0){
		P1=P;
		P2=5.0*P;
	}
	else if(V<V0){
		P1=0.2*P;
		P2=P;
	}
 for(i=0;i<i_max;i++){
		P=0.5*(P1+P2);
		VLE_isobaric_fixphase(T, P, vapor, liquid, solid, f);
		V=vapor.volume+liquid.volume;
		if(abs(V0-V)<err)
			break;
		else if(V<V0) 
			P2=P;
		else if(V>V0)
			P1=P;
	}
 if(i==i_max){
   cout << "VLE_isochoric_fixphase() max.it.=" << i << " reached, dV=" << abs(V-V0) << " rho=";
   cout << vapor.density*1000 << "\n";
 }
}

/**************************************************************************/
/* ROCKFLOW - Funktion: ReactionDeactivation()                            */
/*                                                                        */
/*                                                                        */
/* Task:                                                                  */
/* Deactivates kinetic reaction calculation at individual nodes based on  */
/* evaluation of reaction rates of the previous time step in local        */
/* neighborhoods around nodes, or based on comparison to previous         */
/*   concentrations at the nod                                            */
/*                                                                        */
/* Programmaenderungen:                                                   */
/* 08/2009     CB         First implementation                            */
/*                                                                        */
/**************************************************************************/
void REACTINT::ReactionDeactivation(long nonodes)
{
  long node;
  int react_t = 10;                     // reactions are calculated at all nodes every react_t timesteps
  size_t sp;
  const size_t Number_of_Components(cp_vec.size());
  if (ReactDeactRelative)
    cout << " Reaction deactivation mode " << ReactDeactMode << " with relative concentration change! " << "\n";
  else
    cout << " Reaction deactivation  mode " << ReactDeactMode << " with absolute concentration change! " << "\n";

  // reactivate all nodes every n time steps
  if (((aktueller_zeitschritt) % react_t) == 0) // when evaluated before loop over nodes, i.e. prepare for this time step
    //if(((aktueller_zeitschritt+1) % react_t ) == 0) {  // when evaluated after  loop over nodes, i.e. prepare for next time step
  for (node = 0; node < nonodes; node++)
  {
    ReactDeact[node] = false;
    if (ReactDeactMode != 3) // for all timesteps, anyway prepare the concentrationmatrix for the next time step, i.e. save the current concentrations after Transport
    for (sp = 0; sp < Number_of_Components; sp++)
      concentrationmatrix[node][sp] = cp_vec[sp]->getProcess()->GetNodeValue(node, sp_varind[sp]);
  }
  // CB Now check if node can be deactivated for the next time step
  else
  {
    // only for mode 1: first calculate for each node the sum of rates and store in vector,
    // which then is accessed in next loop over nodes, to evaluate the sum of sum of rates
    // over the neighbours
    if (ReactDeactMode == 1)
    for (node = 0; node < nonodes; node++) {
      React_dCdT[node] = 0;
      for (sp = 0; sp < Number_of_Components; sp++) {
        //This is C of last time step after reactions, i.e. the old time level for this time step
        double Concentration = cp_vec[sp]->getProcess()->GetNodeValue(node, (sp_varind[sp] - 1));
        //This is C of previous time step after transport only, which was stored in matrix
        double Concentration_old = concentrationmatrix[node][sp];
        // relative concentration change
        double maxi = 1;
        if (ReactDeactRelative)
          maxi = DMAX(Concentration_old, Concentration);
        if (maxi < ReactDeactCThresh)
          React_dCdT[node] += 0.0;
        else {
          React_dCdT[node] += fabs((Concentration_old - Concentration) / maxi) / dt; // normalized by current local concentration
          // and now prepare concentrationmatrix for next time step, i.e. save current concentrations after Transport
          concentrationmatrix[node][sp] = cp_vec[sp]->getProcess()->GetNodeValue(node, sp_varind[sp]);
        }
      }
    }

    // this is the check, if a node may be deactivated for this time step
    for (node = 0; node < nonodes; node++) {
      double sumReact_dCdT = 0;
      switch (ReactDeactMode) {
      case 1:       // loop over no of connected nodes and their respective neighbours
        for (size_t nn = 0; nn < ReactNeighborhood[node].size(); nn++) {
          int node_idx = ReactNeighborhood[node][nn];
          sumReact_dCdT += React_dCdT[node_idx];
        }
        if (ReactDeactRelative)
          sumReact_dCdT /= double(ReactNeighborhood[node].size());
        break;
      case 2:       // compare with C after transport of last timestep, loop over all components
        for (sp = 0; sp < Number_of_Components; sp++) {
          //This is C of current time step after transport
          double Concentration = cp_vec[sp]->getProcess()->GetNodeValue(node, sp_varind[sp]);
          //This is C of previous time step after transport only
          double Concentration_old = concentrationmatrix[node][sp];
          double maxi = 1;
          if (ReactDeactRelative)
            maxi = DMAX(Concentration_old, Concentration);
          if (maxi > ReactDeactCThresh)
            sumReact_dCdT += fabs((Concentration - Concentration_old) / maxi);
          // and now prepare the concentrationmatrix for the next time step, i.e. save the current concentrations after Transport
          concentrationmatrix[node][sp] = Concentration;
        }
        break;
      case 3:       // compare with C after reaction of last timestep, loop over all components
        for (sp = 0; sp < Number_of_Components; sp++){
          //This is C of current time step after transport
          double Concentration = cp_vec[sp]->getProcess()->GetNodeValue(node, sp_varind[sp]);
          //This is C of previous time step after transport & reaction
          double Concentration_old = cp_vec[sp]->getProcess()->GetNodeValue(node, (sp_varind[sp] - 1));
          double maxi = 1;
          if (ReactDeactRelative)
            maxi = DMAX(Concentration_old, Concentration);
          if (maxi > ReactDeactCThresh)
            sumReact_dCdT += fabs((Concentration - Concentration_old) / maxi);
        }
        break;
      default:
        break;
      }

      // check if deactivation criterion is met
      if (ReactDeactRelative) sumReact_dCdT /= double(Number_of_Components);
      if (sumReact_dCdT < ReactDeactEpsilon) ReactDeact[node] = true;  // negligible change, deactivate the node for the next time step
      else ReactDeact[node] = false;  // sufficient change, reactivate the node for the next time step
    }  // nnodes                                
  }  //else

  // Resets the reaction rates vector, only in case of model 1
  if (ReactDeactMode == 1)
  for (node = 0; node < nonodes; node++)
    React_dCdT[node] = 0;

}

/**************************************************************************/
/* ROCKFLOW - Funktion: ReactDeactReset_dCdT()                            */
/*                                                                        */
/*                                                                        */
/* Task:                                                                  */
/* Sets the C after reaction of last time step as new C after reaction    */
/* for deactivated nodes                                                          */
/*                                                                        */
/* Programmaenderungen:                                                   */
/* 11/2009     CB         First implementation                            */
/*                                                                        */
/**************************************************************************/
void REACTINT::ReactDeactSetOldReactionTerms(long nonodes)
{
  long node;
  int sp;
  int Number_of_Components = (int)cp_vec.size();
  double Concentration;

  for (node = 0; node < nonodes; node++)
  {
    if (ReactDeact[node] == true)
    for (sp = 0; sp < Number_of_Components; sp++)
    {
      // Get the C after reactions of last time step (old time level, index = 0)
      Concentration = cp_vec[sp]->getProcess()->GetNodeValue(node, sp_varind[sp] - 1);
      // Set this C as the new concentration after reaction
      cp_vec[sp]->getProcess()->SetNodeValue(node, sp_varind[sp],
        Concentration);
    }
  }
}

/**************************************************************************/
/* ROCKFLOW - Funktion: ReactDeactPlotFlagsToTec()                        */
/*                                                                        */
/*                                                                        */
/* Task:                                                                  */
/* Prints flags for reaction deactivation in tecplot format               */
/* for all nodes                                                          */
/*                                                                        */
/* Programmaenderungen:                                                   */
/* 11/2009     CB         First implementation                            */
/*                                                                        */
/**************************************************************************/
void REACTINT::ReactDeactPlotFlagsToTec()
{
  CFEMesh* m_msh = fem_msh_vector[0]; //SB: ToDo hart gesetzt
  std::string eleType;

  if (m_msh->getNumberOfLines() > 0)
    eleType = "QUADRILATERAL";
  if (m_msh->getNumberOfQuads() > 0)
    eleType = "QUADRILATERAL";
  if (m_msh->getNumberOfHexs() > 0)
    eleType = "BRICK";
  if (m_msh->getNumberOfTris() > 0)
    eleType = "QUADRILATERAL";
  if (m_msh->getNumberOfTets() > 0)
    eleType = "TETRAHEDRON";
  if (m_msh->getNumberOfPrisms() > 0)
    eleType = "BRICK";

  if (ReactDeactFlag)
  {
    std::string filename(FileName + "_Deactivated_nodes.tec");
    std::ofstream aus;
    if (aktueller_zeitschritt == 1)
      aus.open(filename.c_str());
    else
      aus.open(filename.c_str(), ios::app);

    const size_t nnodes(m_msh->nod_vector.size());
    const size_t nele(m_msh->ele_vector.size());

    aus << "VARIABLES = " << "\"x\"" << " " << "\"y\"" << " " << "\"z\""
      << "\"active\"" << "\n";
    aus << "ZONE T=" << "\"aktueller_zeitschritt=" << aktueller_zeitschritt
      << "\"";
    aus << ", N=" << nnodes << ", E=" << nele << " F=FEPOINT, ET="
      << eleType << "\n";

    for (size_t i = 0; i < nnodes; i++)
    {
      double const* coord = (m_msh->nod_vector[i])->getData(); // Coordinates(coord);
      aus << coord[0] << " " << coord[1] << " " << coord[2] << " ";
      if (is_a_CCBC[i] == true)
        aus << 0 << "\n";
      else if (ReactDeact[i] == true)
        aus << 0 << "\n";
      else
        aus << 1 << "\n";
    }
    for (size_t i = 0; i < nele; i++)
      m_msh->ele_vector[i]->WriteIndex_TEC(aus);

    aus.close();
  }
}


/**************************************************************************
Reaction-Method:
Task: Copy concentrations after reaction step to symmetric nodes
Implemented only for radial models with one hex element in theta direction
Then it is sufficient to simulate flow and conservative transport at all nodes,
but reactions only at one side of the element. Then copy to other side of element.
Programing:
//SB 01.2011 SB First implementation
**************************************************************************/
void REACTINT::CopyConcentrations(void){

  long i, l, k, nnodes, nele;
  size_t j;
  CFEMesh* m_msh = fem_msh_vector[0]; //SB: ToDo hart gesetzt
  CRFProcess *m_pcs = NULL;
  //CNode* m_nod0 = NULL, *m_nod1 = NULL, *m_nod2 = NULL, *m_nod3 = NULL;
  MeshLib::CElem* m_ele = NULL;
  vec<MeshLib::CNode*> ele_nodes(8);
  Math_Group::vec <long>vec_nod_index(8);
  const double * coord0, *coord1, *coord2, *coord3;
  bool err_flag = true; // if true, all is OK, if false, do not use this option
  int copy_direction = -1;
  long i_copy_from, i_copy_to;
  int nidx0;
  long nidx = 0;

  if (copy_concentrations == false)
    return;

  cout << " CopyConcentrations " << "\n";
  // get number of nodes
  nnodes = (long)m_msh->nod_vector.size();
  nele = (long)m_msh->ele_vector.size();

  // first time called, Check if mesh is appropriate
  if (copy_nodes.size() < 1){
    // BATCH
    if (batch){
      if (nnodes == 2 && nele == 1 && m_msh->ele_vector[0]->GetElementType() == 1){
        // for a single line element
        if (is_a_CCBC[0] == true && is_a_CCBC[1] == false)
          batch = true;
        else if (is_a_CCBC[0] == false && is_a_CCBC[1] == true)
          batch = true;
        else{
          cout << "\n" << " Copy Concentrations only for a single (batch) line element with one deactivated node!"
            << " - Option switched off " << "\n";
          batch = copy_concentrations = false;
        }
      }
    }
    // RADIAL
    else if ((radial) || (TwoDinThreeD)){
      // Only for hex-elements
      for (i = 0; i<nele; i++) {
        if (m_msh->ele_vector[i]->GetElementType() != 3){
          cout << "\n" << " Copy Cons only for hex elements ! - Option switched off " << "\n";
          radial = TwoDinThreeD = copy_concentrations = false;
          break;
        }
      }

      // Now test for radial flow model along x-Axis
      // Test is if model is symmetric to x-Axis; not tested is, if it is one column
      m_ele = m_msh->ele_vector[0];      // get first element
      m_ele->GetNodes(ele_nodes);        // get element nodes
      coord0 = ele_nodes[0]->getData(); // Get coordinates of first four element nodes ("upper" part of element)
      coord1 = ele_nodes[1]->getData(); // CB 03/11
      coord2 = ele_nodes[2]->getData();
      coord3 = ele_nodes[3]->getData();

      /*       *3-----*2
      |\     |\
      | *0---+-*1
      *-|----* |
      \|     \|
      *------*        */

      // check orientation of element: all should have same z-value
      if ((coord0[2] != coord1[2]) || (coord0[2] != coord2[2]) || (coord0[2] != coord3[2])){
        cout << " Error: Not identical z coordinates of element nodes of element 0 " << "\n";
        err_flag = false;
      }
      // node0 and node3 should have same x value, as well as node1 and node2
      if ((coord0[0] != coord3[0]) || (coord1[0] != coord2[0])){
        cout << " Error: Not identical x coordinates of element nodes (0/3) or (1/2) of element 0 " << "\n";
        err_flag = false;
      }
      // this is different for radial and TwoDinThreeD
      if (radial){
        // node0 and node3 should have same y=-y value, as well as node1 and node2
        if ((coord0[1] != -1.0*coord3[1]) || (coord1[1] != -1.0*coord2[1])){ // this is just for radial flow models
          cout << " Error: Not y_i = -y_j for element nodes (0/3) or (1/2) of element 0 " << "\n";
          err_flag = false;
        }
      }
      else if (TwoDinThreeD){
        // node0 and node1 should have same y value, as well as node2 and node3
        if ((coord0[1] != coord1[1]) || (coord2[1] != coord3[1])){ // this is for regular hex grid flow models
          cout << " Error: Not y_i = -y_j for element nodes (0/3) or (1/2) of element 0 " << "\n";
          err_flag = false;
        }
      }
      // everything fine so far?
      if (err_flag == false) {
        cout << "\n" << " Copy Cons - Grid not appliccable ! - Option switched off " << "\n";
        radial = TwoDinThreeD = copy_concentrations = false;
      }

      // Now check, if the correct nodes are switched off and determine copy direction, only for radial
      err_flag = true;       // reset for reuse
      if (radial){
        if (is_a_CCBC[ele_nodes[0]->GetIndex()] == true) { // element node 0 is not calculated
          if (is_a_CCBC[ele_nodes[1]->GetIndex()] != true) err_flag = false;
          if (is_a_CCBC[ele_nodes[2]->GetIndex()] != false) err_flag = false;
          if (is_a_CCBC[ele_nodes[3]->GetIndex()] != false) err_flag = false;
          copy_direction = 1;
        }
        if (is_a_CCBC[ele_nodes[3]->GetIndex()] == true) { // element node 0 is not calculated
          if (is_a_CCBC[ele_nodes[2]->GetIndex()] != true) err_flag = false;
          if (is_a_CCBC[ele_nodes[1]->GetIndex()] != false) err_flag = false;
          if (is_a_CCBC[ele_nodes[0]->GetIndex()] != false) err_flag = false;
          copy_direction = 0;
        }
        // everything fine so far?
        if (err_flag == false) {
          cout << "\n" << " Copy Cons - Please switch off one model side for reaction calculation ! - Option switched off " << "\n";
          radial = TwoDinThreeD = copy_concentrations = false;
        }
      }
    }// ((radial) || (TwoDinThreeD)){

  }// end of mesh and geometry checks

  // now quit, if somethings wrong
  if (copy_concentrations == false)
    return;



  // first time called, build data structures and set up index vector
  if (copy_nodes.size() < 1){

    // Make vector for all nodes for mapping
    for (i = 0; i<nnodes; i++)
      copy_nodes.push_back(i); // map each node to itself by default

    // Now, prepare copy node relations
    if (batch){
      for (i = 0; i<nnodes; i++)
      if (is_a_CCBC[i] == true && i>0) copy_nodes[i] = i - 1;
      else if (is_a_CCBC[i] == true) copy_nodes[i] = i + 1;
    }// batch
    else if (radial){
      // Go through all elements and set copy values respectively 
      for (i = 0; i<nele; i++) {
        m_ele = m_msh->ele_vector[i]; // get element
        m_ele->GetNodeIndeces(vec_nod_index);
        if (copy_direction == 1){ // element nodes 0, 1, 4, 5 not calculated
          // copy concentrations 3 -> 0
          i_copy_from = vec_nod_index[3];
          i_copy_to = vec_nod_index[0];
          copy_nodes[i_copy_to] = i_copy_from;
          // copy concentrations 2 -> 1
          i_copy_from = vec_nod_index[2];
          i_copy_to = vec_nod_index[1];
          copy_nodes[i_copy_to] = i_copy_from;
          // copy concentrations 7 -> 4
          i_copy_from = vec_nod_index[7];
          i_copy_to = vec_nod_index[4];
          copy_nodes[i_copy_to] = i_copy_from;
          // copy concentrations 6 -> 5
          i_copy_from = vec_nod_index[6];
          i_copy_to = vec_nod_index[5];
          copy_nodes[i_copy_to] = i_copy_from;
        }
        if (copy_direction == 0){ // element nodes 2, 3, 6, 7 not calculated
          // copy concentrations 0 -> 3
          i_copy_from = vec_nod_index[0];
          i_copy_to = vec_nod_index[3];
          copy_nodes[i_copy_to] = i_copy_from;
          // copy concentrations 1 -> 2
          i_copy_from = vec_nod_index[1];
          i_copy_to = vec_nod_index[2];
          copy_nodes[i_copy_to] = i_copy_from;
          // copy concentrations 4 -> 7
          i_copy_from = vec_nod_index[4];
          i_copy_to = vec_nod_index[7];
          copy_nodes[i_copy_to] = i_copy_from;
          // copy concentrations 5 -> 6
          i_copy_from = vec_nod_index[5];
          i_copy_to = vec_nod_index[6];
          copy_nodes[i_copy_to] = i_copy_from;
        }
        if (i == 0) {
          cout << "\n" << "copy_nodes[i] = " << "\n";
          for (int ii = 0; ii<8; ii++) cout << "          " << ii << " = " << copy_nodes[ii] << "\n";
        }
      }
    } // radial
    else if (TwoDinThreeD){
      for (i = 0; i<nnodes; i++){
        if (is_a_CCBC[i] == false){ // reactions are quantified at this node
          // get node coordinates  
          coord0 = m_msh->nod_vector[i]->getData();
          // check all neighbour nodes
          for (j = 0; j<m_msh->nod_vector[i]->getNumConnectedNodes(); j++){
            nidx = m_msh->nod_vector[i]->getConnectedNodes()[j];
            coord1 = m_msh->nod_vector[nidx]->getData();
            // check orientation of nodes: should have same x + z and y of neighbour should be larger
            if ((coord0[0] == coord1[0]) && (coord0[2] == coord1[2]) && (coord0[1] < coord1[1])){
              // Found! 
              // Check if reaction is really switched off at this neighbour
              if (is_a_CCBC[nidx] == true) //set the neighbour index
                //copy_nodes[i]=nidx;
                copy_nodes[nidx] = i;
              break;
            }
          } // done
          // now check for success:
          if (copy_nodes[nidx] == nidx)
            cout << "Warning in CKinReactData::CopyConcentrations:!"
            << " No Neighbour found for copy nodes " << i << "\n";
        }// is_a_CCBC[i] == false -> node with reactions
      }// node loop
    }// TwoDinThreeD

  }// if first time


  // Now, Copy concentrations in all cases
  // Go through all transport processes and copy values for all nodes
  // Where to put this? -> ToDo CB in new concept
  // this is the same for batch or radial or 2Din3D
  for (i = 0; i<(int)pcs_vector.size(); i++){
    m_pcs = pcs_vector[i];
    if (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT) { // if it is a MASS_TRANSPORT process
      l = 0;
      // get index of concentration
      for (j = 0; j<(size_t)m_pcs->pcs_number_of_primary_nvals; j++){
        nidx0 = m_pcs->GetNodeValueIndex(m_pcs->pcs_primary_function_name[j]) + 1; // new time level
        for (k = 0; k<nnodes; k++){
          i_copy_to = k;
          i_copy_from = copy_nodes[i_copy_to];
          if (i_copy_to != i_copy_from){ // only do for differing node indicess 
            m_pcs->SetNodeValue(i_copy_to, nidx0, m_pcs->GetNodeValue(i_copy_from, nidx0));
            l++;
          }
        }
      }
    }
  }
  cout << "     at " << l << " of " << nnodes << " nodes." << "\n";
  // end function
}

/**************************************************************************
Reaction-Method:
Task: Calculate impact of CO2 concentration on fluid density in density model 23  - CO2 disolved in salt water
Programing:
// ABM / JOD  11.2015:  
**************************************************************************/

double REACTINT::CalcDeltaDensityFromCO2(double CO2Dissolved, vector<double> &CO2_vector, long Index)
{
	double density_new, density_old;
	double unitfactor_l = 1.0, unitfactor_s = 1.0;
	double mNaCl;

	if (unitconversion)
		CalcUnitConversionFactors(0, &unitfactor_l, &unitfactor_s, true);// mol/m³l --> mol/kg
	// for salt concentration in molality
	mNaCl = c_salt * unitfactor_l * 1000;
	// ----------------------------------------
	//old time step
	if (CO2_vector[Index] >= 1.0e-12)
		density_old = density::CO2brine(GetTemperature(0), GetPressure(0), mNaCl, CO2_vector[Index] * unitfactor_l); // pressure and temperature from node 0
	else
		density_old = 0;
	//new time step
	CO2_vector[Index] = max(1.0e-12, CO2Dissolved);   //  store CO2 concentration value for next time step
	density_new = density::CO2brine(GetTemperature(0), GetPressure(0), mNaCl, CO2_vector[Index] * unitfactor_l); // pressure and temperature from node 0
	
	return 1000 * (density_new - density_old);  // Delta_density
}