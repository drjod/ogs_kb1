/**************************************************************************
   FEMLib - Object: Source Terms ST
   Task: class implementation
   Programing:
   01/2003 OK Implementation
   last modified
**************************************************************************/
#ifndef rf_st_new_INC
#define rf_st_new_INC

// FEM
#include "DistributionInfo.h" // TF
#include "GeoInfo.h"                              // TF
#include "LinearFunctionData.h" // TF
#include "ProcessInfo.h"                          // TF
#include "fem_ele.h"
#include "Constrained.h"

#include "wellDoubletControl.h"
#include "OGS_contraflow.h"
#include <list>

class CNodeValue;
class CGLPolyline;
class CGLLine;
class Surface;
class Polyline;

namespace process                                 //WW
{
class CRFProcessDeformation;
};

namespace MathLib
{
template<class T1, class T2> class InverseDistanceInterpolation;
}

namespace MeshLib
{
class CFEMesh;
class MeshNodesAlongPolyline;
}


class SourceTerm;


struct Threshold
{
	enum Type { no, lower, upper};
	enum Scheme {_explicit, _implicit};

	Type type;
	Scheme scheme;
	double value;
	double delta;
	std::string process;
	int verbosity;
};

struct StorageRate
{
	std::string process;
	bool apply; // true: is flux from transport, false: is not (default)
	double absMaximum;
	double inputValue;
	int verbosity;
	std::string inlet_geometry_name;
	std::string outlet_geometry_name;
	std::vector<long> inlet_msh_node_numbers;
	std::vector<long> outlet_msh_node_numbers;
	std::vector<double> inlet_msh_node_areas;
	std::vector<double> outlet_msh_node_areas;
	double inlet_totalArea, outlet_totalArea;
};



typedef struct
{
	std::vector<double> value_reference;
	std::vector<double> last_source_value;
	//double value_store[10][5000];
	double** value_store;                 //[PCS_NUMBER_MAX*2]First ref no processes(two fields per process..time, value), second ref no values
} NODE_HISTORY;

class CSourceTerm : public ProcessInfo, public GeoInfo, public DistributionInfo
{
	GeoInfo* geoInfo_connected;
	GeoInfo* geoInfo_threshold;  // JOD 2018-02-20
	GeoInfo* geoInfo_storageRateInlet;  // JOD 2018-02-22
	GeoInfo* geoInfo_storageRateOutlet;
	//GeoInfo* geoInfo_wellDoublet_HE;  // heat exchanger JOD 2018-2-13
	GeoInfo* geoInfo_wellDoublet_well1_aquifer;  // warm well -  JOD 2018-6-13
	GeoInfo* geoInfo_wellDoublet_well2_aquifer;  // cold well
	GeoInfo* geoInfo_wellDoublet_well1_liquidBC;  // warm well
	GeoInfo* geoInfo_wellDoublet_well2_liquidBC;  // cold well

	//std::string well1_geometry_name_HE; // heat exchanger JOD 2018-2-13
	std::string well1_geometry_name_aquifer;  // JOD 2018-08-09
	std::string well2_geometry_name_aquifer;
	std::string well1_geometry_name_liquidBC;
	std::string well2_geometry_name_liquidBC;


	OGS_WDC* ogs_WDC;  // pointer to vector entry in pcs
	OGS_contraflow* ogs_contraflow;  // pointer to vector entry in pcs
	int wdc_connector_materialGroup;
	double wdc_connector_normaldirectionVector[3];
public:
	CSourceTerm();
	CSourceTerm(const SourceTerm* st);
	~CSourceTerm();

	std::ios::pos_type Read(std::ifstream* in,
	                        const GEOLIB::GEOObjects & geo_obj,
	                        const std::string & unique_name);
	void Write(std::fstream*);



	void SetNOD2MSHNOD(std::vector<long> & nodes, std::vector<long> & conditional_nodes);

	/**
	 * @param nodes
	 * @param conditional_nodes
	 */
	// 09/2010 TF
	void SetNOD2MSHNOD(const std::vector<size_t>& nodes,
	                   std::vector<size_t>& conditional_nodes) const;

	double GetNodeLastValue(long n, int idx); // used only once in a global in rf_st_new
	                                          // used only once in a global in rf_st_new
	void SetNodePastValue(long n, int idx, int pos, double value);
	// used only once in a global in rf_st_new
	void SetNodeLastValue(long n, int idx, double value);
	double GetNodePastValue(long, int, int); // used only once in a global in rf_st_new
	                                         // used only once in a global in rf_st_new
	double GetNodePastValueReference(long n, int idx);
	// used only in sourcetermgroup
	void SetSurfaceNodeVectorConditional(std::vector<long> & sfc_nod_vector,
	                                     std::vector<long> & sfc_nod_vector_cond);
	void SetSurfaceNodeVectorConnected(std::vector<long> & sfc_nod_vector,
		std::vector<long> & sfc_nod_vector_cond);                       // JOD 2/2015
	void SetPolylineNodeVectorConnected(std::vector<long> & sfc_nod_vector,
		std::vector<long> & sfc_nod_vector_cond);                       // JOD 10/2018
	// used only once in sourcetermgroup
	void InterpolatePolylineNodeValueVector(CGLPolyline* m_ply,
	                                        std::vector<double> & Distribed,
	                                        std::vector<double> & ply_nod_vector);
	void InterpolatePolylineNodeValueVector(
	        std::vector<double> const & nodes_as_interpol_points,
	        std::vector<double>& node_values) const;

	void SetNodeValues(const std::vector<long> & nodes, const std::vector<long> & nodes_cond,
	                   const std::vector<double> & node_values, int ShiftInNodeVector); // used only in sourcetermgroup

	void SetNOD();

	//23.02.2009. WW
	void DirectAssign(const long ShiftInNodeVector);

	// KR / NB
	
	// including climate data into source terms
	const MathLib::InverseDistanceInterpolation<GEOLIB::PointWithID*, GEOLIB::Station*> *getClimateStationDistanceInterpolation() const { return this->_distances; };
	const std::vector<GEOLIB::Station*>& getClimateStations() const { return this->_weather_stations; }; //NB


	//03.2010. WW
	std::string DirectAssign_Precipitation(double current_time);

	double getCoupLeakance () const;

	const std::vector<double>& getDistribution () const { return DistribedBC; }

	/**
	 * REMOVE CANDIDATE only for compatibility with old GEOLIB version
	 * @return
	 */
	const std::string & getGeoName() const;

	double getGeoNodeValue() const { return geo_node_value; } //KR

	int CurveIndex;
	std::vector<int> element_st_vector;

	double st_rill_height, coup_given_value, coup_residualPerm; // JOD 
	double sorptivity, constant, rainfall, rainfall_duration, moistureDeficit /*1x*/;
	bool node_averaging, distribute_volume_flux; // JOD 
	bool no_surface_water_pressure, explicit_surface_water_pressure; // JOD 
	
	bool isCoupled () const { return _coupled; }
	bool isConnected() const { return connected_geometry; }  // JOD 2/2015
	bool hasThreshold() const { return threshold_geometry; }
	bool calculatedFromStorageRate() const { return storageRate_geometry; }

	double getNormalDepthSlope () const { return normaldepth_slope; }
	bool everyoneWithEveryone; // take all nodes from surface and connect with all nodes of other surface  JOD 2015-11-18

	// constrain a ST by other process
	bool isConstrainedST() const { return _isConstrainedST; }
	Constrained const & getConstrainedST(std::size_t constrainedID) const { return _constrainedST[constrainedID]; }
	std::size_t getNumberOfConstrainedSTs() const { return _constrainedST.size(); }
	bool isCompleteConstrainST(std::size_t constrainedID) const { return _constrainedST[constrainedID]._isCompleteConstrained; }
	bool getCompleteConstrainedSTStateOff(std::size_t constrainedID) const { return _constrainedST[constrainedID]._completeConstrainedStateOff; }
	void setCompleteConstrainedSTStateOff(bool status, std::size_t i) { _constrainedST[i]._completeConstrainedStateOff = status; }

	void pushBackConstrainedSTNode(std::size_t constrainedID, bool constrainedStatus) { _constrainedST[constrainedID]._constrainedNodes.push_back(constrainedStatus); }
	bool getConstrainedSTNode(std::size_t constrainedID, int position) const { return _constrainedST[constrainedID]._constrainedNodes[position]; }
	void setConstrainedSTNode(std::size_t constrainedID, bool node_status, int position) { _constrainedST[constrainedID]._constrainedNodes[position] = node_status; }
	std::size_t getNumberOfConstrainedSTNodes(std::size_t constrainedID) const { return _constrainedST[constrainedID]._constrainedNodes.size(); }
	int getConstrainedSTNodesIndex(int position) const { return _constrainedSTNodesIndices[position]; }
	void setConstrainedSTNodesIndex(int st_node_index, int position) { _constrainedSTNodesIndices[position]=st_node_index; }


	std::size_t getSTVectorGroup() const { return _st_vector_group; }
	void setSTVectorGroup(int group) { _st_vector_group = group; }

	const std::string& getFunctionName () const { return fct_name; }
	int getFunctionMethod () const { return fct_method; }

	const std::vector<int>& getPointsWithDistribedST () const { return PointsHaveDistribedBC; }
	const std::vector<double>& getDistribedST() const { return DistribedBC; }
	bool isAnalytical () const { return analytical; }
	bool isPressureBoundaryCondition () const { return pressureBoundaryCondition; }
	void SetPressureBoundaryConditionFlag (bool pBC) { pressureBoundaryCondition = pBC;}
	double GetAnalyticalSolution(long location);
	size_t getNumberOfTerms () const { return number_of_terms; }
	void setMaxNumberOfTerms (size_t max_no_terms) { _max_no_terms = max_no_terms; }
	void setNumberOfAnalyticalSolutions (size_t no_an_sol) { _no_an_sol = no_an_sol; }
	double GetRelativeInterfacePermeability(CRFProcess* m_pcs,  // JOD
                                              double head,
                                              long msh_node);

	double CheckThreshold(const double &value, const CNodeValue* cnodev) const;  // JOD 2018-1-31
	double CalculateFromStorageRate(const double &value, const CNodeValue* cnodev) const;
	double apply_wellDoubletControl(const double &value, const CNodeValue* cnodev,
			const double& aktuelle_zeit, const CRFProcess* m_pcs);  // JOD 2018-06-14
	double apply_contraflow(const double &value, const double& aktuelle_zeit, CRFProcess* m_pcs, double* eqs_rhs);  // JOD 2019-7-30
	bool channel, channel_width, air_breaking;
	double air_breaking_factor, air_breaking_capillaryPressure, air_closing_capillaryPressure;
	int geo_node_number;

	double* nodes;
	std::vector<int> node_number_vector;
	std::vector<double> node_value_vector;
	std::vector<int> node_renumber_vector;
	std::string tim_type_name;
    std::string interpolation_method;			//BG
    int TimeInterpolation;						//BG

	std::string pcs_type_name_cond;
	std::string pcs_pv_name_cond;
	
	int getTimeContrCurve() {return time_contr_curve; } //SB:02.2014 get bc ativity controlled curve
	std::string getTimeContrFunction() {return time_contr_function; } //SB:02.2014 get bc ativity controlled curve
  std::string getFunction() { return fct_name; } //SB:02.2014 get bc ativity controlled curve

	int getSubDomainIndex () const { return _sub_dom_idx; }

	std::string fname;

	  double getTransferCoefficient() { return transfer_coefficient; } //TN
	  double getValueSurrounding() { return value_surrounding;} //TN
	  std::vector<double> get_node_value_vectorArea(){ return node_value_vectorArea;} //TN

	  std::vector<long> st_node_ids;
	  // NNNC
	  bool connected_geometry;			// SB 02/2015
	  std::string connected_geometry_name;
	  std::vector<long> connected_nodes_idx;
	  int connected_geometry_verbose_level;
	  double connected_geometry_exchange_term;  // leakance
	  double connected_geometry_offset;   // not used
	  int connected_geometry_mode;
	  long connected_geometry_ref_element_number;        //  JOD 2015-11-18 - mode 2
	  double connected_geometry_minimum_velocity_abs;    //                      
	  double connected_geometry_reference_direction[3];  //

	  bool threshold_geometry;  // JOD 2018-02-20
	  std::string threshold_geometry_name;
	  long msh_node_number_threshold;

	  bool storageRate_geometry;  // JOD 2018-02-22

private:                                          // TF, KR
	void ReadDistributionType(std::ifstream* st_file);
	void ReadGeoType(std::ifstream* st_file,
	                 const GEOLIB::GEOObjects& geo_obj,
	                 const std::string& unique_name);

	void CreateHistoryNodeMemory(NODE_HISTORY* nh);
	void DeleteHistoryNodeMemory();

	double geo_node_value;
	// active state is controlled by time curve SB
	int time_contr_curve;
	std::string time_contr_function;

	/**
	 * is the source term coupled with another source term
	 */
	bool _coupled;

	double normaldepth_slope;             // used only once in a global in rf_st_new

	/// Subdomain index for excvation simulation
	// 14.12.2010. WW
	int _sub_dom_idx;

	int fct_method;
	std::string fct_name;

	LinearFunctionData* dis_linear_f;     //24.8.2011. WW

	bool analytical;                      //2x?
	bool pressureBoundaryCondition; // pressure load boundary condition

	size_t number_of_terms;
	size_t _max_no_terms;                 // used only once in a global in rf_st_new
	size_t _no_an_sol;
	int analytical_material_group;        // used only once in a global in rf_st_new
	int resolution;                       // used only once in a global in rf_st_new
	double analytical_diffusion;          // used only once in a global in rf_st_new
	double analytical_porosity;           // used only once in a global in rf_st_new
	double analytical_tortousity;         // used only once in a global in rf_st_new
	double analytical_linear_sorption_Kd; // used only once in a global in rf_st_new
	double analytical_matrix_density;     // used only once in a global in rf_st_new
	double factor;

	  double transfer_coefficient; //TN - for DIS_TYPE TRANSFER_SURROUNDING
	  double value_surrounding; //TN - for DIS_TYPE TRANSFER_SURROUNDING

	std::string nodes_file;
	int msh_node_number;
	std::string msh_type_name;
	std::vector<int> PointsHaveDistribedBC;
	std::vector<double> DistribedBC;
	std::vector<double> node_value_vectorArea;
	std::vector<double*> normal2surface;
	std::vector<double*> pnt_parameter_vector;
	// 03.2010. WW
	long start_pos_in_st;
	double* GIS_shape_head;               // 07.06.2010. WW
	std::vector<double> precip_times;
	std::vector<std::string> precip_files;

	friend class CSourceTermGroup;
	friend class process::CRFProcessDeformation; //WW
	friend class ::CRFProcess;            //WW

	std::string geo_name;
	double _coup_leakance;

	// including climate data into source terms
	MathLib::InverseDistanceInterpolation<GEOLIB::PointWithID*, GEOLIB::Station*> *_distances; // NB
	std::vector<GEOLIB::Station*> _weather_stations; //NB

	bool _isConstrainedST;
	std::vector<Constrained> _constrainedST;
	//std::vector<bool> _constrainedSTNodes;
	std::vector<int> _constrainedSTNodesIndices;

	std::size_t _st_vector_group;

	Threshold threshold; // JOD 2018-1-31
	StorageRate storageRate;
	bool ignore_axisymmetry; // JOD 2018-08-15
	bool assign_to_element_edge;  // JOD 2018-9-28
};

class CSourceTermGroup
{
public:
	CSourceTermGroup()                    //WW
	{
	}
	void Set(CRFProcess* m_pcs, const int ShiftInNodeVector, std::string this_pv_name = "");
	//WW    std::vector<CNodeValue*>group_vector;
	/**
	 * \brief process type for the physical process
	 * possible values are
	 * <table>
	 * <tr><td>LIQUID_FLOW</td> <td>H process (incompressible flow)</td></tr>
	 * <tr><td>GROUNDWATER_FLOW</td> <td>H process (incompressible flow)</td></tr>
	 * <tr><td>RIVER_FLOW</td> <td>H process (incompressible flow)</td></tr>
	 * <tr><td>RICHARDS_FLOW</td> <td>H process (incompressible flow)</td></tr>
	 * <tr><td>OVERLAND_FLOW</td> <td>process (incompressible flow)</td></tr>
	 * <tr><td>GAS_FLOW</td> <td>H process (compressible flow)</td></tr>
	 * <tr><td>TWO_PHASE_FLOW</td> <td>H2 process (incompressible/compressible flow)</td></tr>
	 * <tr><td>COMPONENTAL_FLOW</td> <td>H2 process (incompressible/compressible flow)</td></tr>
	 * <tr><td>HEAT_TRANSPORT</td> <td>T process (single/multi-phase flow)</td></tr>
	 * <tr><td>DEFORMATION</td> <td>M process (single/multi-phase flow)</td></tr>
	 * <tr><td>MASS_TRANSPORT</td> <td>C process (single/multi-phase flow)</td></tr>
	 * </table>
	 */
	std::string pcs_name;
	std::string pcs_type_name;            //OK
	std::string pcs_pv_name;              //OK
	MeshLib::CFEMesh* m_msh;
	MeshLib::CFEMesh* m_msh_cond;
	//WW    std::vector<CSourceTerm*>st_group_vector; //OK
	//WW double GetConditionalNODValue(int,CSourceTerm*); //OK
	//WW double GetRiverNODValue(int,CSourceTerm*, long msh_node); //MB
	//WW double GetCriticalDepthNODValue(int,CSourceTerm*, long msh_node); //MB
	//WW double GetNormalDepthNODValue(int,CSourceTerm*, long msh_node); //MB JOD
	//WW Changed from the above
	//    double GetAnalyticalSolution(CSourceTerm *m_st,long node_number, std::string process);//CMCD
	// TRANSFER OF DUAL RICHARDS
	std::string fct_name;                 //YD
private:
	// JOD
	void SetPNT(CRFProcess* m_pcs, CSourceTerm* m_st, const int ShiftInNodeVector);
	// JOD
	void SetLIN(CRFProcess* m_pcs, CSourceTerm* m_st, const int ShiftInNodeVector);
	//OK
	void SetPLY(CSourceTerm* st, int ShiftInNodeVector);
	// JOD
	void SetDMN(CSourceTerm* m_st, const int ShiftInNodeVector);
	// JOD
	void SetSFC(CSourceTerm* m_st, const int ShiftInNodeVector);
	// JOD
	void SetCOL(CSourceTerm* m_st, const int ShiftInNodeVector);

	// JOD
//      void SetPolylineNodeVector(CGLPolyline* m_ply, std::vector<long>&ply_nod_vector);

	void SetPolylineNodeVectorConditional(CSourceTerm* st,
	                                      std::vector<long>&ply_nod_vector,
	                                      std::vector<long>&ply_nod_vector_cond);

	/**
	 * 09/2010 TF
	 * @param st
	 * @param ply_nod_vector
	 * @param ply_nod_vector_cond
	 */
	void SetPolylineNodeVectorConditional(CSourceTerm* st,
	                                      std::vector<size_t>& ply_nod_vector,
	                                      std::vector<size_t>& ply_nod_vector_cond);

      void SetPolylineNodeValueVector(CSourceTerm* st, CGLPolyline * old_ply,
    		  const std::vector<long>& ply_nod_vector,
    		  std::vector<long>& ply_nod_vector_cond, std::vector<double>& ply_nod_val_vector);

	/**
	 * 09/2010 / 03/2011 TF
	 * @param st
	 * @param ply_nod_vector
	 * @param ply_nod_vector_cond
	 * @param ply_nod_val_vector
	 */
	void SetPolylineNodeValueVector(CSourceTerm* st,
	                                std::vector<long> const & ply_nod_vector,
	                                std::vector<long> const & ply_nod_vector_cond,
	                                std::vector<double>& ply_nod_val_vector) const;

	// JOD
	void SetSurfaceNodeVector(Surface* m_sfc, std::vector<long>&sfc_nod_vector);
	void SetSurfaceNodeVector(GEOLIB::Surface const* sfc,
		std::vector<std::size_t> & sfc_nod_vector);
	void SetSurfaceNodeValueVector( CSourceTerm* m_st,
	                                Surface* m_sfc,
	                                std::vector<long> const &sfc_nod_vector,
	                                std::vector<double>&sfc_nod_val_vector);
	void AreaAssembly(const CSourceTerm* const st, const std::vector<long>& ply_nod_vector_cond,
	                  std::vector<double>&  ply_nod_val_vector) const;
	void DistributeVolumeFlux(CSourceTerm* st, std::vector<long> const & ply_nod_vector, // 5.3.07 
		                      std::vector<double>& ply_nod_val_vector);
  
	// CB JOD MERGE //  
	//void CSourceTermGroup::WriteNodeConnections(); 
	void WriteNodeConnections(); // WTP

};

extern CSourceTermGroup* STGetGroup(std::string pcs_type_name,std::string pcs_pv_name);
extern std::list<CSourceTermGroup*> st_group_list;

/**
 * read source term file
 * @param file_base_name base file name (without extension) containing the source terms
 * @param geo_obj object of class GEOObjects managing the geometric entities
 * @param unique_name unique name to access the geometric entities in geo_obj
 * @return true if source terms found in file, else false
 */
bool STRead(const std::string& file_base_name,
            const GEOLIB::GEOObjects& geo_obj,
            const std::string& unique_name);

extern void STWrite(std::string);
#define ST_FILE_EXTENSION ".st"
//extern void STCreateFromPNT();
extern std::vector<CSourceTerm*> st_vector;
extern void STDelete();
void STCreateFromLIN(std::vector<CGLLine*>);
CSourceTerm* STGet(std::string);
extern void STGroupDelete(std::string pcs_type_name,std::string pcs_pv_name);
extern void STGroupsDelete(void);                 //Haibing;
extern size_t aktueller_zeitschritt;
extern double aktuelle_zeit;
extern std::vector<std::string>analytical_processes;
//OK
extern CSourceTerm* STGet(const std::string&, const std::string&, const std::string&);

// WW moved here
//CMCD, WW
extern double GetAnalyticalSolution(long node_number, CSourceTerm* m_st);
//extern  void GetRiverNODValue(double& value, CNodeValue* cnodev, const CSourceTerm* m_st);
extern double GetConditionalNODValue(CSourceTerm* m_st, CNodeValue* cnodev);
//MB
extern void GetCriticalDepthNODValue(double& value, CSourceTerm*, long msh_node);
// JOD
extern void GetCouplingNODValue(double& value, CSourceTerm* m_st, CNodeValue* cnodev);
// JOD
extern void GetCouplingNODValueNewton(double& value, CSourceTerm* m_st, CNodeValue* cnodev);
#if !defined(USE_PETSC)
// && !defined(NEW_EQS) // && defined(other parallel libs)//03~04.3012. WW
//#ifndef NEW_EQS                                   //WW. 06.11.2008
                                                  //MB JOD
extern void GetNormalDepthNODValue(double& value, CSourceTerm*, long msh_node);
// JOD
extern void GetCouplingNODValuePicard(double& value, CSourceTerm* m_st, CNodeValue* cnodev);
extern void GetCouplingNODValueConvectiveForm(double& value, CSourceTerm* m_st, CNodeValue* cnodev);  // JOD 2020-03-25
#endif
// JOD
extern double CalcCouplingValue(double factor,
                                double h_this,
                                double h_cond,
                                double z_cond,
                                CSourceTerm* m_st);
// JOD
extern void GetCouplingNODValueMixed(double& value, CSourceTerm* m_st, CNodeValue* cnodev);
// JOD
extern void GetCouplingFieldVariables(CRFProcess* m_pcs_this,  
                                      CRFProcess* m_pcs_cond,
									  double* h_this,
                                      double* h_cond,
                                      double* h_shifted,
                                      double* h_averaged,
                                      double* z_this,
                                      double* z_cond,
                                      CSourceTerm* m_st,
                                      CNodeValue* cnodev, 
									  long msh_node_number, 
									  long msh_node_number_cond, 
									  double gamma);

// JOD
extern void GetPhilipNODValue(double& value, const CSourceTerm* m_st);
// JOD
extern void GetGreenAmptNODValue(double& value, CSourceTerm* m_st, long msh_node);
// JOD
extern void GetNODValue(double& value, CNodeValue* cnodev,CSourceTerm* m_st);
void IncorporateConnectedGeometries(double& value, CNodeValue* cnodev, CSourceTerm* m_st);// JOD 2/2015
extern void GetNODHeatTransfer(double& value, CSourceTerm* st, long geo_node); //TN

extern void DomainIntegration(CRFProcess* m_pcs,
	                       const std::vector<long> & nodes_in_dom,
	                       std::vector<double> & node_value_vector);

extern void FaceIntegration(MeshLib::CFEMesh* m_msh,
                     std::vector<long> const & nodes_on_sfc,
                     std::vector<double> & node_value_vector,
					 Surface* m_surface, FiniteElement::DistributionType disType, int ele_gauss_points);
	
extern void EdgeIntegration(MeshLib::CFEMesh* m_msh,
	                     const std::vector<long> & nodes_on_ply,
	                     std::vector<double> & node_value_vector, 
			     FiniteElement::DistributionType dis_type, 
			     FiniteElement::PrimaryVariable prim_val, 
			     bool flag_ignore_axisymmetry, bool flag_is_bc);


double get_average(CRFProcess* m_pcs, std::vector<long> vec, long ndx);
#endif
