#ifndef pixel_cuh
#define pixel_cuh

#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#define CUDA_DEV __device__
#else
#define CUDA_HOSTDEV
#define CUDA_DEV
#endif

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
		CUDA_HOSTDEV int getColor();
		// compute the length of the intersection of Segment e
		// with this pixel (may be 0); if x, y are given, store the coordinates
		// of the midpoint of the segment in x, y; DOES NOT SET x, y if length is 0
		CUDA_DEV double intersectionLength(Segment &e, double *x = NULL, double *y = NULL);
		// return true if this pixel contains point p (inside or on the boundary)
		CUDA_DEV bool containsPoint(Point &p);

		// compute area of the intersection of this pixel with Triangle t (may be 0);
		// store the boundary of the intersection in polygon if available along
		// with the size of the boundary
		CUDA_DEV double intersectionArea(Triangle t, Point* polygon = NULL, int* size = NULL);
};

#endif
