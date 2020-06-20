#ifndef constant_cuh
#define constant_cuh

#include <assert.h>
#include <map>
#include <vector>
#include <array>
#include "triangle.cuh"
#include "CImg.h"
#include "parallelInt.cuh"

using namespace std;
using namespace cimg_library;

/**
 * compute a piecewise constant coarse triangular approximation of an image
 */

class ConstantApprox {
	private:
		static const ApproxType APPROXTYPE = constant;
		int maxX, maxY; // dimensions of image
		double stepSize; // size of gradient descent step
		Pixel *pixArr; // pixel (x, y) is pixArr[x * maxY + y]
		vector<Point> *points; // store vertices of triangulation
		Triangle *triArr; // store triangles of triangulation
		int numTri; // number of triangles
		map<Point*, double> gradX; // map points to gradient x values
		map<Point*, double> gradY; // map points to gradient y values
		double *colors; // colors[i] is the color of triangle triArr[i]
		double *results; // used to hold temporary results on device (minimize costly allocations)

	public:
		// initialize a constant approximation triangulation on img
		// with input triangulation points, faces
		ConstantApprox(CImg<unsigned char> *img, vector<Point> *points, vector<array<int, 3>> &triangleInd, double step);
		// deallocate all the shared space
		~ConstantApprox();
		// compute energy of triangulation at this point in time
		double computeEnergy();
		// compute gradient at this instant, updating gradX and gradY
		void computeGrad();
		// store gradient values in gradX, gradY of energy over triangle triArr[t]
		// of the point in t with index movingPt
		void gradient(Triangle *triArr, int t, int movingPt, double *gradX, double *gradY);
		// move points according to gradient values and return true
		// if movement was successful, i.e. no triangle inverts under the process
		bool gradUpdate();
		// in the case that a triangle inverts, undo all the changes at this step
		// and halve the stepSize
		void undo();
		// update approximation value for each triangle
		void updateApprox();
		// run the entire procedure for either maxIter iterations or
		// until change in energy is at most eps
		void run(int maxIter = 1000, double eps = 0.001);

		// return data for image to be displayed
		double getStep();
		// get points of triangulation
		vector<Point> getVertices();
		// get colors for each triangle
		vector<array<double, 3>> getColors();
};

#endif
