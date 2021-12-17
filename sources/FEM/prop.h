#ifndef FEM_PROP_H_
#define FEM_PROP_H_

#include <string>
#include "msh_mesh.h"

// interface class for CMediumProperties, CSolidProperties
// it is used with distributed properties (via files)

class Properties  // JOD 2021-12-16
{
protected:
	MeshLib::CFEMesh* _mesh; //OK

public:
	std::string permeability_file;        //SB //OK/MB string permeability_dis_type_file,BW: X Direction; taken by JOD 2020-3-20
	std::string permeability_Y_file;        //SB //OK/MB string permeability_dis_type_file,BW: Y Direction
	std::string permeability_Z_file;        //SB //OK/MB string permeability_dis_type_file,BW: Z Direction;
	std::string porosity_file;            //OK/MB
	std::string geo_area_file;            //OK
	std::string file_name_conductivity;	 // JOD 2021-12-16
	std::string file_name_capacity;	 // JOD 2021-12-16

	virtual int getNumber() const {return -1;}
public:
	MeshLib::CFEMesh* getMesh(void) const { return _mesh; }
    void setMesh( MeshLib::CFEMesh* m_msh) { _mesh = m_msh; }
};

#endif /* FEM_PROP_H_ */
