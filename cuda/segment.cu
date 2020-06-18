#include "segment.cuh"

// helper function for determining if t in [a, b] where order of a, b is unknown

__device__ bool isBetween(const double &t, const double &a, const double &b) {
	return (a <= t && t <= b) || (b <= t && t <= a);
}

__device__ Segment::Segment() {}

__device__ Segment::Segment(Point *a, Point *b) : endpoint1(a), endpoint2(b) {}

double Segment::length() {
	double x1 = endpoint1->getX();
	double y1 = endpoint1->getY();
	double x2 = endpoint2->getX();
	double y2 = endpoint2->getY();
	return pow(pow(x1-x2,2) + pow(y1-y2,2), 0.5);
}

__device__ void parametrize(Segment &e, Segment &f, double *t1, double *t2, double *det) {
	// parametrize and represent as matrix equation to be solved: 
	// t1 * x0 + (1-t1) * x1 = t2 * x2 + (1-t2) * x3
    // (x0-x1) * t1 + (x3-x2) * t2 = x3 - x1
	double arr[4];
	// first column (matches t1)
	arr[0] = e.endpoint1->getX() - e.endpoint2->getX();
	arr[2] = e.endpoint1->getY() - e.endpoint2->getY();
	// second column (matches t2)
	arr[1] = f.endpoint2->getX() - f.endpoint1->getX();
	arr[3] = f.endpoint2->getY() - f.endpoint1->getY();
	// create 2x2 matrix
	Matrix mat(arr, 2, 2);
	
	double determinant = mat.determinant();
	Matrix adj = mat.adjugate();
	// target vector
	Matrix target(f.endpoint2->getX() - e.endpoint2->getX(), f.endpoint2->getY() - e.endpoint2->getY());
	// compute scaled solution 
	Matrix result = adj.multiply(target);
	*det = determinant;
	*t1 = result.get(0,0);
	*t2 = result.get(1,0);
}

__device__ bool Segment::intersects(Segment &other) {
	double t1 = 0;
	double t2 = 0;
	double det = 0;
	parametrize(*this, other, &t1, &t2, &det);
	if(det == 0) return false;
	// otherwise need to check t1, t2 are both between 0 and determinant
	return isBetween(t1, 0, det) && isBetween(t2, 0, det);
}

__device__ Point Segment::getIntersection(Segment &other) {
	double t1 = 0;
	double t2 = 0;
	double det = 0;
	parametrize(*this, other, &t1, &t2, &det);
	double x = (endpoint1->getX() * t1 + endpoint2->getX() * (det - t1)) / det;
	double y = (endpoint1->getY() * t1 + endpoint2->getY() * (det - t1)) / det;

	return Point(x, y);
}

__device__ Matrix Segment::unitNormal() {
	double deltaX = endpoint2->getX() - endpoint1->getX();
	double deltaY = endpoint2->getY() - endpoint1->getY();
	double unitX = deltaX / length();
	double unitY = deltaY / length();
	// rotate pi/2 clockwise
	return Matrix(unitY, -unitX);
}

__device__ Matrix Segment::scaledNormal() {
	double deltaX = endpoint2->getX() - endpoint1->getX();
	double deltaY = endpoint2->getY() - endpoint1->getY();
	return Matrix(deltaY/2, -deltaX/2);
}