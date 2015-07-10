/*
 * MeshQualityShortestLongestRatio.h
 *
 *  Created on: Mar 3, 2011
 *      Author: TF
 */

#ifndef MESHQUALITYSHORTESTLONGESTRATIO_H_
#define MESHQUALITYSHORTESTLONGESTRATIO_H_

#include "MeshQualityChecker.h"

namespace MeshLib
{
class MeshQualityShortestLongestRatio : public MeshQualityChecker
{
public:
	MeshQualityShortestLongestRatio(CFEMesh const* const mesh);
	virtual ~MeshQualityShortestLongestRatio () {}

	virtual void check ();

private:
	double checkTriangle (GEOLIB::Point const* const a,
	                      GEOLIB::Point const* const b,
	                      GEOLIB::Point const* const c) const;
	double checkQuad (GEOLIB::Point const* const a,
	                  GEOLIB::Point const* const b,
	                  GEOLIB::Point const* const c,
	                  GEOLIB::Point const* const d) const;
	double checkTetrahedron (GEOLIB::Point const* const a,
	                         GEOLIB::Point const* const b,
	                         GEOLIB::Point const* const c,
	                         GEOLIB::Point const* const d) const;
	double checkPrism (std::vector<GEOLIB::Point*> const & pnts) const;
	double checkHexahedron (std::vector<GEOLIB::Point*> const & pnts) const;
};
}

#endif /* MESHQUALITYSHORTESTLONGESTRATIO_H_ */
