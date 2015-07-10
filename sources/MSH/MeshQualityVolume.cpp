/*
 * MeshQualityVolume.cpp
 *
 *  Created on: Mar 3, 2011
 *      Author: TF
 */

#include "MeshQualityVolume.h"

namespace MeshLib
{
MeshQualityVolume::MeshQualityVolume(
        CFEMesh const* const mesh) :
	MeshQualityChecker(mesh)
{ }

void MeshQualityVolume::check()
{
	// get all elements of mesh
	const std::vector<MeshLib::CElem*>& msh_elem(_mesh->getElementVector());

	size_t error_count(0);

	for (size_t k(0); k < msh_elem.size(); k++)
	{
		MshElemType::type elem_type (msh_elem[k]->GetElementType());
		if (elem_type == MshElemType::LINE
		    || elem_type == MshElemType::TRIANGLE
		    || elem_type == MshElemType::QUAD)
		{
            _mesh_quality_measure[k] = -1.0;
            continue;
        }

        double volume (msh_elem[k]->calcVolume());
        if (volume > _max)
            _max = volume;
        if (volume < sqrt(fabs(std::numeric_limits<double>::min()))) {
			errorMsg(msh_elem[k], k);
			error_count++;
		} else if (volume < _min)
            _min = volume;
        _mesh_quality_measure[k] = volume;
	}

	std::cout << "MeshQualityVolume::check() minimum: " << _min
	          << ", max_volume: " << _max << std::endl;
	if (error_count > 0)
		std::cout << "Warning: " << error_count << " elements with zero volume found." <<
		std::endl;
}

} // end namespace MeshLib
