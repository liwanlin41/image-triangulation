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
		Point *points; // store vertices of triangulation
		int numPoints;
		vector<array<int, 3>> faces; // hold triangle connections
		Triangle *triArr; // store triangles of triangulation
		int numTri; // number of triangles
		map<Point*, double> gradX; // map points to gradient x values
		map<Point*, double> gradY; // map points to gradient y values
		double *imageInt; // hold integrals of image dA over each triangle
		double *grays; // grays[i] is the luminance of triangle triArr[i]
		/*
		double *reds; // red value of triangle triArr[i]
		double *greens;
		double *blues;
		*/
		double ds; // step size for integration approximation
		ParallelIntegrator integrator; // do all the integrations

		map<array<int, 2>, vector<int>> edgeBelonging; // map edge to indices of triangles containing edge
		// edge represented by sorted indices of points

		// compute energy change associated with subdividing each edge at its midpoint
		void computeEdgeEnergies(double *edgeEnergies);

	public:
		// create an approximation instance on img with given stepsize and sampling rate
		ConstantApprox(CImg<unsigned char> *img, double step, double ds = 0.1);
		// deallocate all the shared space
		~ConstantApprox();

		// initialize the triangulation on this approximation using a coarse grid,
		// sampling once every pixelRate pixels
		void initialize(int pixelRate);
		// initialize a constant approximation triangulation on img
		// with input triangulation points, faces
		void initialize(vector<Point> *points, vector<array<int, 3>> &triangleInd);

		// compute energy of triangulation at this point in time
		double computeEnergy();
		// compute gradient at this instant, updating gradX and gradY
		void computeGrad();
		// store gradient values in gradX, gradY of energy over triangle triArr[t]
		// of the point in t with index movingPt; imInt = integral of image * dA
		// over triangle to avoid computing multiple times
		void gradient(int t, int movingPt, double imInt, double *gradX, double *gradY);
		// move points according to gradient values and return true
		// if movement was successful, i.e. no triangle inverts under the process
		bool gradUpdate();
		// in the case that a triangle inverts, undo all the changes at this step
		// and halve the stepSize
		void undo();
		// update approximation value for each triangle
		void updateApprox();
		// run one full step of the procedure while tracking energy values
		// return stepsize used in this step
		double step(double &prevEnergy, double &newEnergy);
		// run the entire procedure for either maxIter iterations or
		// until change in energy is at most eps
		void run(int maxIter = 1000, double eps = 0.001);

		// handle adaptive retriangulation to finer mesh

		// greedily subdivide the top n edges
		void subdivide(int n);

		// return data for image to be displayed
		double getStep();
		// get points of triangulation
		vector<Point> getVertices();
		// get edges of triangulation
		vector<array<int, 3>> getFaces();
		// get colors for each triangle
		vector<array<double, 3>> getColors();
};

#endif
