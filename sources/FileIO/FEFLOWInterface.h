
#ifndef FEFLOW_INTERFACE_H
#define FEFLOW_INTERFACE_H

#include "GEOObjects.h"
#include "msh_mesh.h"

class QDomElement;
class QString;

class FEFLOWInterface
{
public:
	// CLASS
	typedef struct
	{
		int problem_class;
		int time_mode;
		int orientation;
		int dimension;
		int n_layers3d;
		int saturation_flag;
		int save_fsize_rreal;
		int save_fsize_creal;
	} FEFLOW_FEM_CLASS;

	// DIMENSION
	typedef struct
	{
		long n_nodes;
		long n_elements;
		int n_nodes_of_element;
		long n_steps;
		long icrank;
		int upwind;
		long obs;
		int optim;
		int aquifer_type;
		int nwca;
		int np_cor;
		int adaptive_mesh;
		int sp_fem_pcs_id;
		int sorption_type;
		int reaction_type;
		int dispersion_type;
	} FEFLOW_FEM_DIM;

	/// Constructor
	FEFLOWInterface(GEOLIB::GEOObjects* geoObjects) : _geoObjects(geoObjects) {}

	MeshLib::CFEMesh* readFEFLOWModelFile(const std::string &filename);

private:
	//void readSuperMesh(std::ifstream &feflow_file, FEFLOW_FEM_CLASS &fem_class, FEFLOW_FEM_DIM &fem_dim, std::vector<GEOLIB::Point*> *points, std::vector<GEOLIB::Polyline*> *lines);
	void readSuperMesh(std::ifstream &feflow_file,
	                   FEFLOW_FEM_CLASS &fem_class,
	                   std::vector<GEOLIB::Point*>** points,
	                   std::vector<GEOLIB::Polyline*>** lines);
	void readPoints(QDomElement &nodesEle, const std::string &tag, int dim, std::vector<GEOLIB::Point*> &points);
	void setMaterialID(MeshLib::CFEMesh* m_msh, std::vector<GEOLIB::Polyline*>* lines);


	GEOLIB::GEOObjects* _geoObjects;
};

#endif
