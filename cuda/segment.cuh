#ifndef segment_cuh
#define segment_cuh

#include <math.h>
#include "point.cuh"
#include "matrix.cuh"

/**
 * represent a point (point1, point2) in the plane
 */

class Segment {
	private:
		Point *endpoint1, *endpoint2;
	public:
		__device__ Segment(); // this just exists to create arrays
		__device__ Segment(Point *a, Point *b);
		__host__ __device__ double length();

		// return true if this segment intersects other in exactly one point
		__device__ bool intersects(Segment &other);
		// return intersection point if this intersects other
		// undefined behavior if intersection does not exist or if segments overlap with positive length
		__device__ Point getIntersection(Segment &other);
		
		// return the 2x1 unit normal to this segment which lies on the right
		// going from endpoint 1 to endpoint 2
		__device__ Matrix unitNormal();
		// return the 2x1 normal to this segment that has length |segment|/2
		__device__ Matrix scaledNormal();

		// helper function that determines parameters at intersection point of e and f,
		// storing as t1, t2; intersection is x0 * (t1 / det) + x1 * (1 - t1/det)
		friend __device__ void parametrize(Segment &e, Segment &f, double *t1, double *t2, double *det);
		friend class Pixel;
};

#endif 
