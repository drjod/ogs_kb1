/*
  Class to define a problem to be solved.
  Programming WW 08.07.2008
*/
#ifndef problem_INC
#define problem_INC


#include <vector>
class CRFProcess;

// GEOLIB
#include "GEOObjects.h"
#include "msh_mesh.h"

extern int logging;

class Problem;
//class CMediumProperties;
class Properties;
//---------------------------------------------------------------------
//Pointers to member functions
typedef  double (Problem::*ProblemMemFn)(void);
#define Call_Member_FN(object,ptrToMember)  ((object)->*(ptrToMember))
//---------------------------------------------------------------------
class Problem
{
   public:
      Problem(char* filename = NULL);
      ~Problem();
      void Euler_TimeDiscretize();
      void RosenBrock_TimeDiscretize() {};
      //
	void SetActiveProcesses();
	void SetTimeActiveProcesses(); // JT2012
	void PCSRestart();
	//
	bool CouplingLoop();
	void PostCouplingLoop();
	void PreCouplingLoop(CRFProcess *m_pcs = NULL);
    bool PrintTimes(){return print_times;};
	// Copy u_n for auto time stepping
    double* GetBufferArray(const bool is_x_k = false)
    {
      if(is_x_k)
        return buffer_array1;
      else			
        return buffer_array;
    }
	int GetCPLMaxIterations() { return cpl_overall_max_iterations; }

	/**
       * get the geometric objects stored in GEOLIB::GEOObjects
       * @return a pointer to an instance of class GEOLIB::GEOObjects
       */
      const GEOLIB::GEOObjects* getGeoObj () const;
      /**
       * Get the name of the project. The name is used by GEOLIB::GEOObjects
       * to access the geometric data.
       * @return the name to acces geometric data
       */
      const std::string& getGeoObjName () const;

#ifdef BRNS
      // BRNS-Coupling: For writing spatially resolved reaction rates at the final iteration,
      // we need to get the timing information.

      double getCurrentTime();
      double getEndTime();
#endif                                         //BRNS

	/// Set rank and size for MPI
    void setRankandSize(const int rank, const int size)
	{
        mrank = rank;  
        msize = size;  
	} 
   private:
      // Time:
      double start_time;
      double end_time;
	double current_time;
	double* buffer_array;
	double* buffer_array1;
	int step_control_type;
	bool last_dt_accepted;		//JT
	bool force_post_node_copy; //JT
	// Mixed time step WW
	double dt0;                           // Save the original time step size

	// Controls
	bool external_coupling_exists;
	int cpl_overall_max_iterations;
	int cpl_overall_min_iterations;
	int loop_process_number;
	size_t max_time_steps;
	//
      int lop_coupling_iterations;
      bool CalcVelocities;
      bool conducted;

      // Print flag
      bool print_result;
      bool print_times;
      // Processes

      std::vector<CRFProcess*> total_processes;
      std::vector<CRFProcess*> transport_processes;
	std::vector<CRFProcess*> multiphase_processes;
	std::vector<CRFProcess*> singlephaseflow_process;
	ProblemMemFn* active_processes;
	std::vector<int> active_process_index;
	std::vector<int> coupled_process_index;
	//
	bool* exe_flag;
	inline int AssignProcessIndex(CRFProcess* m_pcs, bool activefunc = true);
	//
      void PCSCreate();
      // Perform processes:
      inline double LiquidFlow();
      inline double RichardsFlow();
      inline double TwoPhaseFlow();
      inline double MultiPhaseFlow();
      inline double PS_Global();                  // 03 2009 PCH
	  inline double MULTI_COMPONENTIAL_FLOW();            // 02 2011 AKS/NB
	  inline double TNEQ();						//07/2013 HS/TN
	  inline double TES();						//07/2013 HS/TN
      inline double GroundWaterFlow();
      inline double ComponentalFlow();
      inline double OverlandFlow();
      inline double AirFlow();
      inline double HeatTransport();
      inline double FluidMomentum();
      inline double RandomWalker();
      inline double MassTrasport();
      inline double PostMassTrasport();
      inline double Deformation();
      // Accessory
	void LOPExecuteRegionalRichardsFlow(CRFProcess* m_pcs_global, int loop_process_number);
	void LOPCalcELEResultants();
	inline void ASMCalcNodeWDepth(CRFProcess* m_pcs);
	void PCSCalcSecondaryVariables();
	bool Check();                     //OK
	void TestOutputEclipse(CRFProcess* m_pcs);	//BG
	void TestOutputDuMux(CRFProcess* m_pcs);	//BG
	void OutputMassOfComponentInModel(std::vector<CRFProcess*> flow_pcs, CRFProcess *transport_pcs); //BG
	void OutputMassOfGasInModel(CRFProcess* m_pcs); //BG

	/**
	 * pointer to an instance of class GEOObjects,
	 * that manages geometric entities
	 */
	GEOLIB::GEOObjects* _geo_obj;         // TF
	/**
	 * project/file name for geometry file,
	 * used to access data in data manager GEOObjects
	 */
	std::string _geo_name;                // TF

	/// rank for MPI
	int mrank;
    /// processor number for MPI  
	int msize;

	static const size_t max_processes = 16;
};

extern bool MODCreate();                          //OK

void GetHeterogeneousFields();
void SetDistributedELEProperties(Properties*, const std::string&, const std::string&, const int&);
void WriteTecplotDistributedProperties(const Properties* const prop);

long GetNearestHetVal2(long EleIndex,
							MeshLib::CFEMesh* m_msh,
                              std::vector <double> xvals,
                              std::vector <double> yvals,
                              std::vector <double> zvals,
                              std::vector <double> mmpvals);
double GetAverageHetVal2(long EleIndex,
						MeshLib::CFEMesh* m_msh,
                         std::vector <double> xvals,
                         std::vector <double> yvals,
                         std::vector <double> zvals,
                         std::vector <double> mmpvals);

#endif
