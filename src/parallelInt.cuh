#ifndef parallel_int_cuh
#define parallel_int_cuh
#include <stdio.h>
#include <assert.h>
#include "pixel.cuh"
#include "triangle.cuh"
#include "cuda_runtime.h"

using namespace std;

// define different supported approximation types
enum ApproxType{constant, linear, quadratic};

// parallelize integral computation

class ParallelIntegrator {
    private:
        // thread setup
        static const int threadsX = 32;
        static const int threadsY = 16;
        static const int threads1D = 1024; // NOTE: changing this will require changes in sumBlock
        dim3 threads2D;
        // two large arrays to do computations
        double *arr;
        double *helper;
        ApproxType approx; // type of approximation to do
        Pixel *pixArr; // reference to the image being approximated
        Point *curTri; // hold vertices of current working triangle
        int maxX, maxY; // size of image
        // true if computations are exact rather than approximate
        bool computeExact;

        // sum the first size values of arr
        double sumArray(int size);

    public:
        ParallelIntegrator();
        // initialize parallel integrator
        // space indicates the amount of computation space needed
        // default to using approximate integrals
        // return true if successful
        bool initialize(Pixel *pix, int xMax, int yMax, ApproxType a, long long space, bool exact = false);
        // free allocated space
        ~ParallelIntegrator();

        // actual integrals

        // for energy, it makes sense to separate approximation types into different functions
        // because the inputs will be greatly different (one array for each coefficient)

        // compute energy integral over tri
        double constantEnergyEval(Triangle *tri, double color, double ds);
        // compute energy by exact integral
        double constantEnergyExact(Triangle *tri, double color);
        // approximate energy using barycentric sampling
        double constantEnergyApprox(Triangle *tri, double color, double ds);

        // compute integral of (f-g)^2 dA over triangle abc where
        // approximation is coeff[0] * a + coeff[1] * b + coeff[2] * c
        double linearEnergyApprox(Point *a, Point *b, Point *c, double *coeffs, double ds);

        // compute exact line integral (v dot n) * f phi ds over triangle tri
        // where FEM basis phi is dependent on approx
        // consider when point with index pt in triArr[t] is moving at velocity (1,0) if isX and (0,1) if !isX
        double lineIntEval(Triangle *tri, int pt, bool isX, double ds);
        // compute by exact integral
        double lineIntExact(Triangle *tri, int pt, bool isX);
        // approximate line integral using one sample every ds length
        double lineIntApprox(Triangle *tri, int pt, bool isX, double ds);

        // compute double integral f phi dA over triangle triArr[t] where FEM basis phi depends on approx
        double doubleIntEval(Triangle *tri, double ds, ColorChannel channel = GRAY);
        // compute by exact integral
        double doubleIntExact(Triangle *tri, ColorChannel channel);
        // approximate by grid with side length ds
        double doubleIntApprox(Triangle *tri, double ds, ColorChannel channel);

        double doubleIntApprox(Point *a, Point *b, Point *c, double ds, ColorChannel channel, int basisIndex);
};

#endif
