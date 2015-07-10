/*
 * MeshQualityEquiAngleSkew.cpp
 *
 *  Created on: Mar 17, 2011
 *      Author: TF
 */

#include "MeshQualityEquiAngleSkew.h"

// MathLib
#include "MathTools.h"

#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

namespace MeshLib
{
MeshQualityEquiAngleSkew::MeshQualityEquiAngleSkew(CFEMesh const* const mesh) :
	MeshQualityChecker(mesh), M_PI_THIRD (M_PI / 3.0), TWICE_M_PI (2 * M_PI)
{}

MeshQualityEquiAngleSkew::~MeshQualityEquiAngleSkew()
{}

void MeshQualityEquiAngleSkew::check ()
{
	// get all elements of mesh
	const std::vector<MeshLib::CElem*>& msh_elem(_mesh->getElementVector());

	for (size_t k(0); k < msh_elem.size(); k++)
	{
		switch (msh_elem[k]->GetElementType())
		{
		case MshElemType::LINE:
			_mesh_quality_measure[k] = -1.0;
			break;
		case MshElemType::TRIANGLE:
			_mesh_quality_measure[k] = checkTriangle (msh_elem[k]);
			break;
		case MshElemType::QUAD:
			_mesh_quality_measure[k] = checkQuad (msh_elem[k]);
			break;
		case MshElemType::TETRAHEDRON:
			_mesh_quality_measure[k] = checkTetrahedron (msh_elem[k]);
			break;
		case MshElemType::HEXAHEDRON:
			_mesh_quality_measure[k] = checkHexahedron (msh_elem[k]);
			break;
		case MshElemType::PRISM:
			_mesh_quality_measure[k] = checkPrism (msh_elem[k]);
			break;
		default:
			break;
		}
	}
}

double MeshQualityEquiAngleSkew::checkTriangle (CElem const* const elem) const
{
	double const* const node0 (elem->GetNode(0)->getData());
	double const* const node1 (elem->GetNode(1)->getData());
	double const* const node2 (elem->GetNode(2)->getData());

	double min_angle (M_PI_2), max_angle (0.0);
	getMinMaxAngleFromTriangle (node0, node1, node2, min_angle, max_angle);

	return 1.0 -
	       std::max((max_angle - M_PI_THIRD) / (M_PI - M_PI_THIRD),
	                (M_PI_THIRD - min_angle) / (M_PI_THIRD));
}

double MeshQualityEquiAngleSkew::checkQuad (CElem const* const elem) const
{
	double const* const node0 (elem->GetNode(0)->getData());
	double const* const node1 (elem->GetNode(1)->getData());
	double const* const node2 (elem->GetNode(2)->getData());
	double const* const node3 (elem->GetNode(3)->getData());

	double min_angle (TWICE_M_PI);
	double max_angle (0.0);

	getMinMaxAngleFromQuad (node0, node1, node2, node3, min_angle, max_angle);

	return 1.0 -
	       std::max((max_angle - M_PI_2) / (M_PI - M_PI_2), (M_PI_2 - min_angle) / (M_PI_2));
}

double MeshQualityEquiAngleSkew::checkTetrahedron (CElem const* const elem) const
{
	double const* const node0 (elem->GetNode(0)->getData());
	double const* const node1 (elem->GetNode(1)->getData());
	double const* const node2 (elem->GetNode(2)->getData());
	double const* const node3 (elem->GetNode(3)->getData());

	double min_angle (M_PI_2);
	double max_angle (0.0);

	// first triangle (0,1,2)
	getMinMaxAngleFromTriangle(node0, node1, node2, min_angle, max_angle);
	// second triangle (0,1,3)
	getMinMaxAngleFromTriangle(node0, node1, node3, min_angle, max_angle);
	// third triangle (0,2,3)
	getMinMaxAngleFromTriangle(node0, node2, node3, min_angle, max_angle);
	// fourth triangle (1,2,3)
	getMinMaxAngleFromTriangle(node1, node2, node3, min_angle, max_angle);

	return 1.0 - std::max((max_angle - M_PI_2) / (M_PI - M_PI_THIRD),
	                      (M_PI_THIRD - min_angle) / (M_PI_THIRD));
}

double MeshQualityEquiAngleSkew::checkHexahedron (CElem const* const elem) const
{
	double const* const node0 (elem->GetNode(0)->getData());
	double const* const node1 (elem->GetNode(1)->getData());
	double const* const node2 (elem->GetNode(2)->getData());
	double const* const node3 (elem->GetNode(3)->getData());
	double const* const node4 (elem->GetNode(4)->getData());
	double const* const node5 (elem->GetNode(5)->getData());
	double const* const node6 (elem->GetNode(6)->getData());
	double const* const node7 (elem->GetNode(7)->getData());

	double min_angle (2 * M_PI);
	double max_angle (0.0);

	// first surface (0,1,2,3)
	getMinMaxAngleFromQuad (node0, node1, node2, node3, min_angle, max_angle);
	// second surface (0,3,7,4)
	getMinMaxAngleFromQuad (node0, node3, node7, node4, min_angle, max_angle);
	// third surface (4,5,6,7)
	getMinMaxAngleFromQuad (node4, node5, node6, node7, min_angle, max_angle);
	// fourth surface (5,1,2,6)
	getMinMaxAngleFromQuad (node5, node1, node2, node6, min_angle, max_angle);
	// fifth surface (5,1,0,4)
	getMinMaxAngleFromQuad (node5, node1, node0, node4, min_angle, max_angle);
	// sixth surface (6,2,3,7)
	getMinMaxAngleFromQuad (node6, node2, node3, node7, min_angle, max_angle);

	return 1.0 -
	       std::max((max_angle - M_PI_2) / (M_PI - M_PI_2), (M_PI_2 - min_angle) / (M_PI_2));
}

double MeshQualityEquiAngleSkew::checkPrism (CElem const* const elem) const
{
	double const* const node0 (elem->GetNode(0)->getData());
	double const* const node1 (elem->GetNode(1)->getData());
	double const* const node2 (elem->GetNode(2)->getData());
	double const* const node3 (elem->GetNode(3)->getData());
	double const* const node4 (elem->GetNode(4)->getData());
	double const* const node5 (elem->GetNode(5)->getData());

	double min_angle_tri (2 * M_PI);
	double max_angle_tri (0.0);

	// first triangle (0,1,2)
	getMinMaxAngleFromTriangle (node0, node1, node2, min_angle_tri, max_angle_tri);
	// second surface (3,4,5)
	getMinMaxAngleFromTriangle (node3, node4, node5, min_angle_tri, max_angle_tri);

	double tri_criterion (1.0 - std::max((max_angle_tri - M_PI_2) / (M_PI - M_PI_THIRD),
	                                     (M_PI_THIRD - min_angle_tri) / (M_PI_THIRD)));

	double min_angle_quad (2 * M_PI);
	double max_angle_quad (0.0);
	// surface (0,3,4,1)
	getMinMaxAngleFromQuad (node0, node3, node4, node1, min_angle_quad, max_angle_quad);
	// surface (2,5,3,0)
	getMinMaxAngleFromQuad (node2, node5, node3, node0, min_angle_quad, max_angle_quad);
	// surface (1,2,5,4)
	getMinMaxAngleFromQuad (node1, node2, node5, node4, min_angle_quad, max_angle_quad);

	double quad_criterion (1.0 - std::max((max_angle_quad - M_PI_2) / (M_PI - M_PI_2),
	                                      (M_PI_2 - min_angle_quad) / (M_PI_2)));

	return std::min (tri_criterion, quad_criterion);
}

void MeshQualityEquiAngleSkew::getMinMaxAngleFromQuad (
        double const* const n0, double const* const n1,
        double const* const n2, double const* const n3,
        double &min_angle, double &max_angle) const
{
	double angle (MathLib::getAngle (n3, n0, n1));
	if (angle < min_angle)
		min_angle = angle;
	if (angle > max_angle)
		max_angle = angle;

	angle = MathLib::getAngle (n0, n1, n2);
	if (angle < min_angle)
		min_angle = angle;
	if (angle > max_angle)
		max_angle = angle;

	angle = MathLib::getAngle (n1, n2, n3);
	if (angle < min_angle)
		min_angle = angle;
	if (angle > max_angle)
		max_angle = angle;

	angle = MathLib::getAngle (n2, n3, n0);
	if (angle < min_angle)
		min_angle = angle;
	if (angle > max_angle)
		max_angle = angle;
}

void MeshQualityEquiAngleSkew::getMinMaxAngleFromTriangle(double const* const n0,
                                                          double const* const n1,
                                                          double const* const n2,
                                                          double &min_angle,
                                                          double &max_angle) const
{
	double angle (MathLib::getAngle (n2, n0, n1));
	if (angle < min_angle)
		min_angle = angle;
	if (angle > max_angle)
		max_angle = angle;

	angle = MathLib::getAngle (n0, n1, n2);
	if (angle < min_angle)
		min_angle = angle;
	if (angle > max_angle)
		max_angle = angle;

	angle = MathLib::getAngle (n1, n2, n0);
	if (angle < min_angle)
		min_angle = angle;
	if (angle > max_angle)
		max_angle = angle;
}
} // end namespace MeshLib
