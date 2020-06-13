#ifndef constant_h
#define constant_h

#include <vector>
#include <map>
#include <assert.h>
#include "triangle.hpp"
#include "lineIntegral.hpp"
#include "doubleIntegral.hpp"
#include "matrix.hpp"

using namespace std;

/**
 * compute a piecewise constant coarse triangular approximation of an image
 */

class ConstantApprox {
    private:
        int maxX, maxY; // dimensions of image
        double stepSize;
        vector<Point> corners; // corners of image
        vector<vector<Pixel>> *image; // image.at(x).at(y) is the pixel located at (x, y)
        vector<Point> points; // store movable vertices of triangulation
        vector<Triangle> triangles; // store triangles of triangulation
        map<Triangle*, vector<Point*>> trianglesToPoints; // map triangles to points of triangle
        map<Point*, double> gradX; // map points to gradient x values
        map<Point*, double> gradY; // map points to gradient y values
        map<Triangle*, double> approx; // map triangles to approximation values
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
        // energy over Triangle t (with provided vertices for space) of movement 
        // of the point with index movingPt
        void gradient(Triangle &t, vector<Point*> &vertices, int movingPt, double *gradX, double *gradY);
        // move points according to gradient values and
        // return the index of the first point whose movement
        // causes a triangle to invert (-1 if none)
        int gradUpdate();
        // in the case that a triangle inverts, undo all the
        // changes at this step and halve the stepSize
        void undo(int ind);
        // update approximation value for each triangle
        void updateApprox();
        // run the entire procedure for either maxIter iterations or until
        // largest gradient norm is at most eps
        void run(int maxIter = 10000, double eps = 0.001);

        // check to ensure dictionary trianglesToPoints is referenced correctly
        void crossCheck();
};

#endif