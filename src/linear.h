#ifndef linear_h
#define linear_h

#include <assert.h>
#include <map>
#include <vector>
#include <array>
#include <set>
#include "approx.h"

using namespace std;
using namespace cimg_library;

/**
 * compute a piecewise constant coarse triangular approximation of an image
 */

class Approx;

class LinearApprox : public Approx {
	private:
		static const ApproxType APPROXTYPE = linear;

	public:
		// create an approximation instance on img with given stepsize and sampling rate
		LinearApprox(CImg<unsigned char> *img, double step, double ds = 0.1);
		// deallocate all the shared space
		~LinearApprox();
		ApproxType getApproxType();

		void reallocateSpace();

		void initialize(int pixelRate);
		void initialize(vector<Point> &points, vector<array<int,3>> &faces);

		void computeGrad();
		void updateApprox();
		double computeEnergy();

		// store gradient values in gradX, gradY of energy over triangle triArr[t]
		// of the point in t with index movingPt; imInt = integral of image * dA
		// over triangle to avoid computing multiple times
		void gradient(int t, int movingPt, double imInt, double *gradX, double *gradY);

		void computeEdgeEnergies(vector<array<double, 3>> *edgeEnergies);

		vector<array<double, 3>> getColors();
        /*
        vector<Point> testPoints();
        vector<array<int, 3>> testFaces();
        */
};

#endif