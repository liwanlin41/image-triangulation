#ifndef segment_cuh
#define segment_cuh

#include "point.cuh"

/**
 * represent a point (point1, point2) in the plane
 */

class Segment {
	private:
		Point *endpoint1, *endpoint2;
	public:
		__device__ Segment(Point *a, Point *b);
		__device__ double length();

		// return true if this segment intersects other in exactly one point
		__device__ bool intersects(Segment other);
		// return intersection point if this intersects other
		// undefined behavior if intersection does not exist or if segments overlap with positive length
		__device__ Point getIntersection(Segment other);

};

#endif 
