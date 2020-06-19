#ifndef pixel_cuh
#define pixel_cuh

#include <nvfunctional>
#include <assert.h>
#include "point.cuh"
#include "segment.cuh"
#include "triangle.cuh"

/**
 * class representing an image pixel as a square, with pixel
 * (x,y) being a unit square centered at (x, y)
 */

class Pixel {
	private:
		int color; // represent luminance of pixel from 0 to 255; may change to support color images later
		int x, y; // center of pixel
		Point corners[4]; // corners of pixel in ccw order
	public:
		// create a pixel centered at (x, y) with color c
		Pixel(int x, int y, int c);
		__host__ __device__ int getColor();
		// compute the length of the intersection of Segment e
		// with this pixel (may be 0); if x, y are given, store the coordinates
		// of the midpoint of the segment in x, y; DOES NOT SET x, y if length is 0
		__device__ double intersectionLength(Segment &e, double *x = NULL, double *y = NULL);
		// return true if this pixel contains point p (inside or on the boundary)
		__device__ bool containsPoint(Point &p);

		// return approximate line integral contribution of this pixel when multiplied
		// with integrand func (exact for func linear); in other words, integrate
		// color times func over the subsegment of e lying in this pixel
		__device__ double lineIntegral(nvstd::function<double(double, double)> &func, Segment &e);

		// approximate line integral contribution over a triangle
		__device__ double lineIntegral(nvstd::function<double(double, double)> &func, Triangle &t);
		
		// compute area of the intersection of this pixel with Triangle t (may be 0);
		// store the boundary of the intersection in polygon if available along
		// with the size of the boundary
		__device__ double intersectionArea(Triangle &t, Point* polygon = NULL, int* size = NULL);

		// return approximate double integral contribution of this pixel when multiplied
		// with integrand func; in other words, integrate color times func over the part
		// of this pixel that lies inside triangle t;
		// approximation by average value of x, y in intersection of pixel and triangle
		__device__ double doubleIntegral(nvstd::function<double(double, double)> &func, Triangle &t);
};

#endif
