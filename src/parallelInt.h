#ifndef parallel_int_h
#define parallel_int_h
#include <cuda.h>
#include <cuda_runtime.h>
#include <functional>
#include <vector>
#include <map>
#include <assert.h>
#include "pixel.hpp"
#include "triangle.hpp"

using namespace std;

// parallelize integral computation

// compute the double integral of func over triangle in image represented by pixArr
double doubleIntEval(function<double(double, double)> func, Pixel *pixArr, int maxX, int maxY, Triangle *triArr, int numTri);
// given the array of pixels and the array of triangles with color map given by colors, return
// the energy of the approximation
double constantEnergyEval(Pixel *pixArr, int maxX, int maxY, Triangle *triArr, map<Triangle*, double> colors);

// compute line integral of func over triangle in image represented by pixVec
double lineIntEval(function<double(double, double)> func, vector<vector<Pixel>> *pixVec, Triangle *triangle);

// for testing
double run();

#endif
