#ifndef triangle_cuh
#define triangle_cuh

#include <math.h>
#include "point.cuh"
#include "segment.cuh"
#include "matrix.cuh"

/**
 * represent a triangle whose vertices are movable
 */

class Triangle {
	private:
		Point* vertices[3];
	public:
		// construct from three point pointers and orient in ccw direction
		Triangle(Point *a, Point *b, Point *c);
		__host__ __device__ double getArea();
		// get signed area based on the order of vertices
		// with ccw direction positive
		__host__ __device__ double getSignedArea();
		// get the change in area when the pth vertex is moving at velocity (vx, vy)
		double dA(int &p, double vx, double vy);
		// get the gradient in the x direction for pth vertex
		double gradX(int &p);
		// get the gradient in the y direction for pth vertex
		double gradY(int &p);

		// determine if triangle contains point p
		__device__ bool contains(Point &p);

		// get vertices of triangle, store in a, b, c
		__device__ void copyVertices(Point *a, Point *b, Point *c);

		// static signed area function
		__host__ __device__ static double getSignedArea(Point *a, Point *b, Point *c);

		friend ostream& operator<<(ostream& os, const Triangle &t);
		
		friend class ConstantApprox;
		friend class Pixel;
};

#endif
