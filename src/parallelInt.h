#ifndef parallel_int_h
#define parallel_int_h
#include <cuda.h>
#include <cuda_runtime.h>
#include <functional>
#include <vector>
#include <assert.h>
#include "pixel.hpp"
#include "triangle.hpp"

using namespace std;

// parallelize integral computation

// compute the double integral of func over triangle in image represented by pixVec
double doubleIntEval(function<double(double, double)> func, vector<vector<Pixel>> *pixVec, Triangle *triangle);

// compute line integral of func over triangle in image represented by pixVec
double lineIntEval(function<double(double, double)> func, vector<vector<Pixel>> *pixVec, Triangle *triangle);

// for testing
double runSum(double *arr, int size);

#endif
