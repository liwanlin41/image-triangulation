#ifndef constant_h
#define constant_h

#include <vector>
#include <array>
#include <map>
#include <assert.h>
#include "triangle.hpp"
#include "lineIntegral.hpp"
#include "doubleIntegral.hpp"
#include "matrix.hpp"
#include "CImg.h"

using namespace std;
using namespace cimg_library;

/**
 * compute a piecewise constant coarse triangular approximation of an image
 */

class ConstantApprox {
    private:
        int maxX, maxY; // dimensions of image
        double stepSize; // size of gradient descent step
        vector<vector<Pixel>> *image; // image.at(x).at(y) is the pixel located at (x, y)
        vector<Point> points; // store vertices of triangulation
        vector<Triangle> triangles; // store triangles of triangulation
        map<Point*, double> gradX; // map points to gradient x values
        map<Point*, double> gradY; // map points to gradient y values
        map<Triangle*, double> approx; // map triangles to approximation values
        vector<array<int, 3>> triangleInd; // triangleInd[i] are the indices in points of vertices of triangles[i]
    public:
        // create a constant approximation triangulation on img with n+4 vertices
        // and a given gradient stepsize
        // TODO: figure out how to initialize this
        ConstantApprox(vector<vector<Pixel>> *img, int n, double step);
        // compute energy of triangulation at this point in time
        double computeEnergy();
        // compute gradient at this instant, updating gradX and gradY
        void computeGrad();
        // store gradient values in gradX, gradY of
        // energy over Triangle t of movement 
        // of the point with index movingPt
        void gradient(Triangle &t, int movingPt, double *gradX, double *gradY);
        // move points according to gradient values and
        // return true if movement was successful, i.e. no triangle inverts
        // under this process
        bool gradUpdate();
        // in the case that a triangle inverts, undo all the
        // changes at this step and halve the stepSize
        void undo();
        // update approximation value for each triangle
        void updateApprox();
        // run the entire procedure for either maxIter iterations or until
        // change in energy is at most eps
        void run(int maxIter = 1000, double eps = 0.0001);

        // return data for image to be displayed
        
        // get points of triangulation
        vector<Point> getVertices();
        // get edges of triangulation as vertex indices
        vector<array<int, 3>> getFaces();
        // get colors for each triangle
        vector<array<double,3>> getColors();
};

#endif