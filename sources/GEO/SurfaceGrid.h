/*
 * SurfaceGrid.h
 *
 *  Created on: Feb 9, 2012
 *      Author: TF
 */

#ifndef SURFACEGRID_H_
#define SURFACEGRID_H_

#include <vector>

// GEOLIB
#include "AxisAlignedBoundingBox.h"

namespace GEOLIB {

// forward declaration
class Surface;
class Triangle;

class SurfaceGrid : public AABB {
public:
	SurfaceGrid(Surface const*const sfc);
	virtual ~SurfaceGrid();

	bool isPntInSurface(const double* pnt, double eps = 0) const;

private:
	double _step_sizes[3];
	double _inverse_step_sizes[3];
	size_t _n_steps[3];
	std::vector<GEOLIB::Triangle const*>* _triangles_in_grid_box;
};

} // end namespace GEOLIB

#endif /* SURFACEGRID_H_ */
