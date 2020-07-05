#ifndef parallel_int_cuh
#define parallel_int_cuh
#include <stdio.h>
#include <assert.h>
#include "pixel.cuh"
#include "triangle.cuh"

using namespace std;

// define different supported approximation types
enum ApproxType{constant, linear, quadratic};

// parallelize integral computation
/*
class ParallelIntegrator {
    private:
        // two large arrays to do computations
        double *arr;
        double *helper;
        ApproxType approx; // type of approximation to do
        Pixel *pixArr; // reference to the image being approximated
        Triangle *triArr; // reference to triangles of mesh
        int maxX, maxY; // size of image
        // true if computations are exact rather than approximate
        bool computeExact;
    public:
        // create a parallel integrator
        // space indicates the amount of computation space needed
        // default to using approximate integrals
        ParallelIntegrator(Pixel *pix, Triangle *tri, int xMax, int yMax, ApproxType a, long long space, bool exact = false);
        // free allocated space
        ~ParallelIntegrator();
        // copy assignment operator
        ParallelIntegrator& ParallelIntegrator::operator=(const ParallelIntegrator &other);

        // actual integrals

        // compute exact line integral (v dot n) * f phi ds over triangle triArr[t]
        // where FEM basis phi is dependent on approx
        // consider when point with index pt in triArr[t] is moving at velocity (1,0) if isX and (0,1) if !isX
        double lineIntEval(int t, int pt, bool isX);
        // approximate line integral using one sample every ds length
        double lineIntApprox(int t, int pt, bool isX, double ds);

        // compute double integral f phi dA
        double
};
*/

// compute the line integral (v dot n) * f phi ds over triangle triArr[t] where FEM basis phi is dependent on ApproxType
// consider when point with index pt in triArr[t] is moving at velocity (1,0) if isX and (0,1) if !isX
// workingTriangle is shared memory to save space when getting vertices of triArr[t]
double lineIntEval(ApproxType approx, Pixel *pixArr, int &maxX, int &maxY, Triangle *triArr, int &t, int &pt, bool isX, double *results0, double *results1, Point *workingTriangle);
double lineIntApprox(ApproxType approx, Pixel *pixArr, int &maxY, Triangle *triArr, int &t, int &pt, bool isX, double *results0, double *results1, double ds, Point *workingTri);

// compute the double integral f phi dA over triangle triArr[t] where FEM basis phi is dependent on approx
// image pixel color channel determined by channel
double doubleIntEval(ApproxType approx, Pixel *pixArr, int &maxX, int &maxY, Triangle *triArr, int &t, double *results0, double *results1, ColorChannel channel = GRAY);
double doubleIntApprox(ApproxType approx, Pixel *pixArr, int &maxY, Triangle *tri, double *results0, double *results1, double &ds, Point *workingTri, ColorChannel channel = GRAY);

// given the array of pixels and the array of triangles with color colors[i] on triangle i,
// return the energy of the approximation
double constantEnergyEval(Pixel *pixArr, int &maxX, int &maxY, Triangle *triArr, double *colors, int &numTri, double *results0, double *results1);
// compute a gridded approximate energy at intervals of ds
double constantEnergyApprox(Pixel *pixArr, int &maxY, Triangle *triArr, double *colors, int &numTri, double *results0, double *results1, double ds, Point *workingTri);

#endif
