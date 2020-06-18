#ifndef point_cuh
#define point_cuh

#include <math.h>
#include <iostream>
#include <iomanip>

using namespace std;

/**
 * represent a point (x,y) on the plane
 * mutable; allows perturbing vertices
 */

class Point {
	private:
		double x, y;
		// determine if point is on edge of image and thus cannot move
		bool borderX, borderY;
	public:
		__host__ __device__ Point(double x, double y, bool borderX = false, bool borderY = false);
		__host__ __device__ double getX() const;
		__host__ __device__ double getY() const;
		// return true if point was constructed on a vertical image edge
		__device__ bool isBorderX() const;
		__device__ bool isBorderY() const;
		__host__ __device__ double distance(Point &other);
		__device__ void move(double deltaX, double deltaY);
		__device__ bool operator==(const Point &other) const;
		__device__ bool operator!=(const Point &other) const;

		friend ostream& operator<<(ostream& os, const Point &p);
};

#endif
