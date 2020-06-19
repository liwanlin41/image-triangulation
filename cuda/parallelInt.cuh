#ifndef parallel_int_cuh
#define parallel_int_cuh
#include <nvfunctional>
#include <assert.h>
#include "pixel.cuh"
#include "triangle.cuh"

using namespace std;

// parallelize integral computation

// compute the double integral of func times image over triArr[t] in image represented by pixArr;
// pixArr, results, triArr must already be shared between host and device
// where results is created to avoid unnecessary memory allocation
double doubleIntEval(nvstd::function<double(double, double)> func, Pixel *pixArr, int &maxX, int &maxY, Triangle *triArr, int &t, double *results);

// given the array of pixels and the array of triangles with color colors[i] on triangle i,
// return the energy of the approximation
double constantEnergyEval(Pixel *pixArr, int &maxX, int &maxY, Triangle *triArr, double *colors, int &numTri, double *results);

// compute line integral of func over triangle in image represented by pixArr
double lineIntEval(nvstd::function<double(double, double)> func, Pixel *pixArr, int &maxX, int &maxY, Triangle *triArr, int &t, double *results);

#endif
