/*
 * Surface.h
 *
 *  Created on: Jan 22, 2010
 *      Author: TF
 */

#ifndef SURFACE_H_
#define SURFACE_H_

#include <vector>

#include "AxisAlignedBoundingBox.h"
#include "GeoObject.h"
#include "Point.h"
#include "Polyline.h"
#include "Triangle.h"
#include "SurfaceGrid.h"

namespace GEOLIB
{
/**
 * \ingroup GEOLIB
 *
 * \brief A Surface is represented by Triangles. It consists of a reference
 * to a vector of (pointers to) points (m_sfc_pnts) and a vector that stores
 * the Triangles consisting of points from m_sfc_pnts.
 * */
class Surface : public GeoObject
{
public:
	Surface (const std::vector<Point*> &pnt_vec);
	virtual ~Surface ();

	/**
	 * adds three indices describing a triangle and updates the bounding box
	 * */
	void addTriangle (size_t pnt_a, size_t pnt_b, size_t pnt_c);

	/// Triangulates a new surface based on closed polyline.
	static Surface* createSurface(const Polyline &ply);

	void deleteSurfaceGrid();

	/**
	 * returns the number of triangles describing the Surface
	 * */
	size_t getNTriangles () const;

	/** \brief const access operator for the access to the i-th Triangle of the surface.
	 */
	const Triangle* operator[] (size_t i) const;

	/**
	 * is the given point in the bounding volume of the surface
	 */
	bool isPntInBV (const double* pnt, double eps = std::numeric_limits<double>::epsilon()) const;

	/**
	 * is the given point pnt located in the surface
	 * @param pnt the point
	 * @return true if the point is contained in the surface
	 */
	bool isPntInSfc (const double* pnt,
		double eps = std::numeric_limits<double>::epsilon()) const;

	const std::vector<Point*>* getPointVec() const { return &_sfc_pnts; }

	/**
	 * method allows access to the internal axis aligned bounding box
	 * @return axis aligned bounding box
	 */
	AABB const & getAABB () const { return _bv; }

	void initSurfaceGrid();
	void calculateTriangleNormals() const;
	double const* getTriangleNormal(const std::size_t triangle_id) const;
	int getTriangleIDOfPoint(const double* pnt) const;

protected:
	/** a vector of pointers to Points */
	const std::vector<Point*> &_sfc_pnts;
	/** position of pointers to the geometric points */
	std::vector<Triangle*> _sfc_triangles;
	/** bounding volume is an axis aligned bounding box */
	AABB _bv;
	/** a helper structure to improve the search */
	SurfaceGrid* _sfc_grid;
};
}

#endif /* SURFACE_H_ */
