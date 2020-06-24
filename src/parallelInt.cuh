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

// compute the line integral (v dot n) * f phi ds over triangle triArr[t] where FEM basis phi is dependent on ApproxType
// consider when point with index pt in triArr[t] is moving at velocity (1,0) if isX and (0,1) if !isX
// workingTriangle is shared memory to save space when getting vertices of triArr[t]
double lineIntEval(ApproxType approx, Pixel *pixArr, int &maxX, int &maxY, Triangle *triArr, int &t, int &pt, bool isX, double *results, Point *workingTriangle);

double lineIntTest(ApproxType approx, Pixel *pixArr, int &maxX, int &maxY, Triangle *triArr, int &t, int &pt, bool isX, double *results, Point *working);

// compute the double integral f phi dA over triangle triArr[t] where FEM basis phi is dependent on approx
double doubleIntEval(ApproxType approx, Pixel *pixArr, int &maxX, int &maxY, Triangle *triArr, int &t, double *results);

// given the array of pixels and the array of triangles with color colors[i] on triangle i,
// return the energy of the approximation
double constantEnergyEval(Pixel *pixArr, int &maxX, int &maxY, Triangle *triArr, double *colors, int &numTri, double *results);

#endif
