#ifndef approx_h
#define approx_h

#include <map>
#include <vector>
#include <array>
#include <set>
#include "triangle.cuh"
#include "parallelInt.cuh"
#include "CImg.h"
#include "cuda.h"
#include "cuda_runtime.h"

using namespace std;
using namespace cimg_library;

/**
 * base class for triangular image approximation
 */

class Approx {
	protected:
		int maxX, maxY; // dimensions of image
		double stepSize; // size of gradient descent step
		double originalStep; // starting step size
		Pixel *pixArr; // pixel (x, y) is pixArr[x * maxY + y]
		Point *points; // store vertices of triangulation
		int numPoints;
		vector<array<int, 3>> faces; // hold triangle connections
		Triangle *triArr; // store triangles of triangulation
		int numTri; // number of triangles
		map<Point*, double> gradX; // map points to gradient x values
		map<Point*, double> gradY; // map points to gradient y values
		double ds; // step size for integration approximation
		ParallelIntegrator integrator; // do all the integrations

		map<array<int, 2>, vector<int>> edgeBelonging; // map edge to indices of triangles containing edge
		// edge represented by sorted indices of points

        // initialize triangulation and integrator
		// initialize the triangulation on this approximation using a coarse grid,
		// sampling once every pixelRate pixels
		void initialize(ApproxType approxtype, int pixelRate);
		// initialize a constant approximation triangulation on img
		// with input triangulation points, faces
		void initialize(ApproxType approxtype, vector<Point> &points, vector<array<int, 3>> &triangleInd);

	public:
		// create an approximation instance on img with given stepsize and sampling rate
		Approx(CImg<unsigned char> *img, double step, double ds = 0.1);
		// deallocate all the shared space
		~Approx();

        // free and reassign space that is not handled in approx 
        virtual void reallocateSpace() = 0;

        // create starting approximation
        virtual void initialize(int pixelRate) = 0;
        virtual void initialize(vector<Point> &points, vector<array<int, 3>> &triangleInd) = 0;

		// compute energy of triangulation at this point in time
		virtual double computeEnergy() = 0;
		// compute gradient at this instant, updating gradX and gradY
		virtual void computeGrad() = 0;
		// move points according to gradient values and return true
		// if movement was successful, i.e. no triangle inverts under the process
		bool gradUpdate();
		// in the case that a triangle inverts, undo all the changes at this step
		// and halve the stepSize
		void undo();
		// update approximation value for each triangle
		virtual void updateApprox() = 0;
		// run one full step of the procedure while tracking energy values
		// return stepsize used in this step
		// stringent determines whether any energy increase is tolerated,
		// defaults true (no tolerance)
		double step(double &prevEnergy, double &newEnergy, bool stringent = true);
		// run the entire procedure for either maxIter iterations or
		// until change in energy is at most eps
		void run(int maxIter = 100, double eps = 0.001);

		// handle adaptive retriangulation to finer mesh

		// greedily subdivide the top n edges
		void subdivide(int n);
		// compute energy change associated with subdividing each edge at its midpoint
		// {i, j, e} in vector means edge (i, j) has energy e
		virtual void computeEdgeEnergies(vector<array<double, 3>> *edgeEnergies) = 0;
		// update subdivided mesh
		void updateMesh(vector<Point> *newPoints, vector<array<int, 3>> *newFaces, set<int> *discardedFaces);

		// return data for image to be displayed
		double getStep();
		// get points of triangulation
		vector<Point> getVertices();
		// get edges of triangulation
		vector<array<int, 3>> getFaces();
        virtual vector<array<double, 3>> getColors() = 0;
};

#endif