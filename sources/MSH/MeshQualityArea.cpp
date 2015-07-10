/*
 * MeshQualityArea.cpp
 *
 * 2011/03/17 KR Initial Implementation
 */

#include "MeshQualityArea.h"
#include "mathlib.cpp"

namespace MeshLib
{
MeshQualityArea::MeshQualityArea(CFEMesh const* const mesh)
	: MeshQualityChecker(mesh)
{}

void MeshQualityArea::check()
{
	// get all elements of mesh
	const std::vector<MeshLib::CElem*>& msh_elem(_mesh->getElementVector());

	size_t nElems(msh_elem.size());
	for (size_t k(0); k < nElems; k++) {
		MshElemType::type elem_type(msh_elem[k]->GetElementType());
		if (elem_type == MshElemType::LINE) {
			_mesh_quality_measure[k] = -1.0;
			continue;
		}

		double area(std::numeric_limits<double>::max());
		if (elem_type == MshElemType::TRIANGLE || elem_type != MshElemType::QUAD) {
			area = msh_elem[k]->calcVolume();
			if (area < sqrt(fabs(std::numeric_limits<double>::min()))) errorMsg(msh_elem[k], k);
		} else {
			size_t nFaces(msh_elem[k]->GetFacesNumber());

			int face_node_index[4];
			for (size_t i = 0; i < nFaces; i++) {
				size_t nNodes = msh_elem[k]->GetElementFaceNodes(i, face_node_index);

				double subarea(0);
				if (nNodes == 3) // face is a triangle
					subarea = ComputeDetTri(msh_elem[k]->GetNode(0)->getData(),
									msh_elem[k]->GetNode(1)->getData(),
									msh_elem[k]->GetNode(2)->getData());
				else if (nNodes == 4) // face is a quad
					subarea = ComputeDetTri(msh_elem[k]->GetNode(0)->getData(),
									msh_elem[k]->GetNode(1)->getData(),
									msh_elem[k]->GetNode(2)->getData()) +
									ComputeDetTri(
									msh_elem[k]->GetNode(2)->getData(),
									msh_elem[k]->GetNode(3)->getData(),
									msh_elem[k]->GetNode(0)->getData());
				else std::cout << "Error in MeshQualityArea::check()" << std::endl;

				if (subarea < sqrt(fabs(std::numeric_limits<double>::min())))
					errorMsg(msh_elem[k], k);
				if (subarea < area) area = subarea;
			}
		}
		// update _min and _max values
		if (_min > area) _min = area;
		if (_max < area) _max = area;
		_mesh_quality_measure[k] = area;
	}
}

} // end namespace MeshLib
