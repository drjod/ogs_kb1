/*
 * MeshQualityChecker.cpp
 *
 *  Created on: Dec 8, 2010
 *      Author: TF
 */

#include "MeshQualityChecker.h"
#include "msh_elem.h"
#include <cmath>

namespace MeshLib
{
MeshQualityChecker::MeshQualityChecker(CFEMesh const* const mesh) :
	_min (std::numeric_limits<double>::max()), _max (std::numeric_limits<double>::min()), _mesh (mesh)
{
	if (_mesh)
		_mesh_quality_measure.resize ((_mesh->getElementVector()).size(), -1.0);
}

BASELIB::Histogram<double> MeshQualityChecker::getHistogram (size_t nclasses) const
{
	if (nclasses == 0) {
		// simple suggestion: number of classes with Sturges criterion
		nclasses = static_cast<size_t>(1 + 3.3 * log (static_cast<float>((_mesh->getElementVector()).size())));
	}

	return BASELIB::Histogram<double>(getMeshQuality(), nclasses, true);
}

void MeshQualityChecker::errorMsg (CElem* elem, size_t idx) const
{
	std::cout <<
	"Error in MeshQualityChecker::check() - Calculated value of element is below double precision minimum."
			  << std::endl;
	std::cout << "Points of " << MshElemType2String(elem->GetElementType()) << "-Element " <<
	idx << ": " << std::endl;
	for (int i(0); i < elem->GetVertexNumber(); i++)
		std::cout << "\t Node " << i << " " <<
		GEOLIB::Point((elem->GetNode(i))->getData()) << std::endl;
}

std::vector<double> const&
MeshQualityChecker::getMeshQuality () const
{
	return _mesh_quality_measure;
}

double MeshQualityChecker::getMinValue() const
{
	return _min;
}

double MeshQualityChecker::getMaxValue() const
{
	return _max;
}

} // end namespace MeshLib
