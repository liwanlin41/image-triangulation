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
        Triangle *triArr; // reference to triangles of mesh
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
        bool initialize(Pixel *pix, Triangle *tri, int xMax, int yMax, ApproxType a, long long space, bool exact = false);
        // free allocated space
        ~ParallelIntegrator();

        // actual integrals

        // for energy, it makes sense to separate approximation types into different functions
        // because the inputs will be greatly different (one array for each coefficient)

        // compute energy integral over triArr[t]
        double constantEnergyEval(double color, int t, double ds);
        // compute energy by exact integral
        double constantEnergyExact(double color, int t);
        // approximate energy using barycentric sampling
        double constantEnergyApprox(double color, int t, double ds);

        // compute exact line integral (v dot n) * f phi ds over triangle triArr[t]
        // where FEM basis phi is dependent on approx
        // consider when point with index pt in triArr[t] is moving at velocity (1,0) if isX and (0,1) if !isX
        double lineIntEval(int t, int pt, bool isX, double ds);
        // compute by exact integral
        double lineIntExact(int t, int pt, bool isX);
        // approximate line integral using one sample every ds length
        double lineIntApprox(int t, int pt, bool isX, double ds);

        // compute double integral f phi dA over triangle triArr[t] where FEM basis phi depends on approx
        double doubleIntEval(int t, double ds, ColorChannel channel = GRAY);
        // compute by exact integral
        double doubleIntExact(int t, ColorChannel channel);
        // approximate by grid with side length ds
        double doubleIntApprox(int t, double ds, ColorChannel channel);
};

#endif
