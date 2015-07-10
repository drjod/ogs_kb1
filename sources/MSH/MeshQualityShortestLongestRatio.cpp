/*
 * MeshQualityShortestLongestRatio.cpp
 *
 *  Created on: Mar 3, 2011
 *      Author: TF
 */

#include "MeshQualityShortestLongestRatio.h"

namespace MeshLib
{
MeshQualityShortestLongestRatio::MeshQualityShortestLongestRatio(
        CFEMesh const* const mesh) :
	MeshQualityChecker(mesh)
{
}

void MeshQualityShortestLongestRatio::check()
{
	// get all elements of mesh
	const std::vector<MeshLib::CElem*>& msh_elem(_mesh->getElementVector());

	for (size_t k(0); k < msh_elem.size(); k++)
	{
		switch (msh_elem[k]->GetElementType())
		{
		case MshElemType::LINE:
			_mesh_quality_measure[k] = 1.0;
			break;
		case MshElemType::TRIANGLE: {
			GEOLIB::Point* a(new GEOLIB::Point(
			                         (msh_elem[k]->GetNode(0))->getData()));
			GEOLIB::Point* b(new GEOLIB::Point(
			                         (msh_elem[k]->GetNode(1))->getData()));
			GEOLIB::Point* c(new GEOLIB::Point(
			                         (msh_elem[k]->GetNode(2))->getData()));
			_mesh_quality_measure[k] = checkTriangle(a, b, c);
			delete a;
			delete b;
			delete c;
			break;
		}
		case MshElemType::QUAD: {
			GEOLIB::Point* a(new GEOLIB::Point(
			                         (msh_elem[k]->GetNode(0))->getData()));
			GEOLIB::Point* b(new GEOLIB::Point(
			                         (msh_elem[k]->GetNode(1))->getData()));
			GEOLIB::Point* c(new GEOLIB::Point(
			                         (msh_elem[k]->GetNode(2))->getData()));
			GEOLIB::Point* d(new GEOLIB::Point(
			                         (msh_elem[k]->GetNode(3))->getData()));
			_mesh_quality_measure[k] = checkQuad(a, b, c, d);
			delete a;
			delete b;
			delete c;
			delete d;
			break;
		}
		case MshElemType::TETRAHEDRON: {
			GEOLIB::Point* a(new GEOLIB::Point(
			                         (msh_elem[k]->GetNode(0))->getData()));
			GEOLIB::Point* b(new GEOLIB::Point(
			                         (msh_elem[k]->GetNode(1))->getData()));
			GEOLIB::Point* c(new GEOLIB::Point(
			                         (msh_elem[k]->GetNode(2))->getData()));
			GEOLIB::Point* d(new GEOLIB::Point(
			                         (msh_elem[k]->GetNode(3))->getData()));
			_mesh_quality_measure[k] = checkTetrahedron(a, b, c, d);
			delete a;
			delete b;
			delete c;
			delete d;
			break;
		}
		case MshElemType::PRISM: {
			std::vector<GEOLIB::Point*> pnts;
			for (size_t j(0); j < 6; j++)
				pnts.push_back(new GEOLIB::Point(
				                       (msh_elem[k]->GetNode(j))->getData()));
			_mesh_quality_measure[k] = checkPrism(pnts);
			for (size_t j(0); j < 6; j++)
				delete pnts[j];
			break;
		}
		case MshElemType::HEXAHEDRON: {
			std::vector<GEOLIB::Point*> pnts;
			for (size_t j(0); j < 8; j++)
				pnts.push_back(new GEOLIB::Point(
				                       (msh_elem[k]->GetNode(j))->getData()));
			_mesh_quality_measure[k] = checkHexahedron(pnts);
			for (size_t j(0); j < 8; j++)
				delete pnts[j];
			break;
		}
		default:
			std::cout <<
			"MeshQualityShortestLongestRatio::check () check for element type "
			          << MshElemType2String(msh_elem[k]->GetElementType())
			          << " not implemented" << std::endl;
		}
	}
}

double MeshQualityShortestLongestRatio::checkTriangle (GEOLIB::Point const* const a,
                                                       GEOLIB::Point const* const b,
                                                       GEOLIB::Point const* const c) const
{
	double len0 (sqrt(MathLib::sqrDist (b,a)));
	double len1 (sqrt(MathLib::sqrDist (b,c)));
	double len2 (sqrt(MathLib::sqrDist (a,c)));

	if (len0 < len1 && len0 < len2)
	{
		if (len1 < len2)
			return len0 / len2;
		else
			return len0 / len1;
	}
	else
	{
		if (len1 < len2)
		{
			if (len0 < len2)
				return len1 / len2;
			else
				return len1 / len0;
		}
		else
		{
			if (len0 < len1)
				return len2 / len1;
			else
				return len2 / len0;
		}
	}
}

double MeshQualityShortestLongestRatio::checkQuad (GEOLIB::Point const* const a,
                                                   GEOLIB::Point const* const b,
                                                   GEOLIB::Point const* const c,
                                                   GEOLIB::Point const* const d) const
{
	double sqr_lengths[4] = {MathLib::sqrDist (b,a),
		                 MathLib::sqrDist (c,b),
		                 MathLib::sqrDist (d,c),
		                 MathLib::sqrDist (a,d)};

	// sort lengths - since this is a very small array we use bubble sort
	for (size_t i(0); i < 4; i++)
		for (size_t j(i + 1); j < 4; j++)
			if (sqr_lengths[i] >= sqr_lengths[j])
				std::swap (sqr_lengths[i], sqr_lengths[j]);

	return sqrt(sqr_lengths[0]) / sqrt(sqr_lengths[3]);
}

double MeshQualityShortestLongestRatio::checkTetrahedron (GEOLIB::Point const* const a,
                                                          GEOLIB::Point const* const b,
                                                          GEOLIB::Point const* const c,
                                                          GEOLIB::Point const* const d) const
{
	double sqr_lengths[6] = {MathLib::sqrDist (b,a), MathLib::sqrDist (c,b),
		                 MathLib::sqrDist (c,a), MathLib::sqrDist (a,d),
		                 MathLib::sqrDist (b,d), MathLib::sqrDist (c,d)};

	// sort lengths - since this is a very small array we use bubble sort
	for (size_t i(0); i < 6; i++)
		for (size_t j(i + 1); j < 6; j++)
			if (sqr_lengths[i] >= sqr_lengths[j])
				std::swap (sqr_lengths[i], sqr_lengths[j]);

	return sqrt(sqr_lengths[0]) / sqrt(sqr_lengths[5]);
}

double MeshQualityShortestLongestRatio::checkPrism (std::vector<GEOLIB::Point*> const & pnts) const
{
	double sqr_lengths[9] = {MathLib::sqrDist (pnts[0],pnts[1]),
		                 MathLib::sqrDist (pnts[1],pnts[2]),
		                 MathLib::sqrDist (pnts[2],pnts[0]),
		                 MathLib::sqrDist (pnts[3],pnts[4]),
		                 MathLib::sqrDist (pnts[4],pnts[5]),
		                 MathLib::sqrDist (pnts[5],pnts[3]),
		                 MathLib::sqrDist (pnts[0],pnts[3]),
		                 MathLib::sqrDist (pnts[1],pnts[4]),
		                 MathLib::sqrDist (pnts[2],pnts[5])};

	// sort lengths - since this is a very small array we use bubble sort
	for (size_t i(0); i < 9; i++)
		for (size_t j(i + 1); j < 9; j++)
			if (sqr_lengths[i] >= sqr_lengths[j])
				std::swap (sqr_lengths[i], sqr_lengths[j]);

	return sqrt(sqr_lengths[0]) / sqrt(sqr_lengths[8]);
}

double MeshQualityShortestLongestRatio::checkHexahedron (std::vector<GEOLIB::Point*> const & pnts)
const
{
	double sqr_lengths[12] = {MathLib::sqrDist (pnts[0],pnts[1]),
		                  MathLib::sqrDist (pnts[1],pnts[2]),
		                  MathLib::sqrDist (pnts[2],pnts[3]),
		                  MathLib::sqrDist (pnts[3],pnts[0]),
		                  MathLib::sqrDist (pnts[4],pnts[5]),
		                  MathLib::sqrDist (pnts[5],pnts[6]),
		                  MathLib::sqrDist (pnts[6],pnts[7]),
		                  MathLib::sqrDist (pnts[7],pnts[4]),
		                  MathLib::sqrDist (pnts[0],pnts[4]),
		                  MathLib::sqrDist (pnts[1],pnts[5]),
		                  MathLib::sqrDist (pnts[2],pnts[6]),
		                  MathLib::sqrDist (pnts[3],pnts[7])};

	// sort lengths - since this is a very small array we use bubble sort
	for (size_t i(0); i < 12; i++)
		for (size_t j(i + 1); j < 12; j++)
			if (sqr_lengths[i] >= sqr_lengths[j])
				std::swap (sqr_lengths[i], sqr_lengths[j]);

	return sqrt(sqr_lengths[0]) / sqrt(sqr_lengths[11]);
}
} // end namespace MeshLib
