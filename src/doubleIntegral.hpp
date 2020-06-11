#ifndef double_integral_h
#define double_integral_h

#include <vector>
#include "triangle.hpp"
#include "pixel.hpp"

using namespace std;

/**
 * static class for evaluating a double integral over a triangle
 */
class DoubleIntegral {
    private:
        DoubleIntegral();
    public:
        // static method for evaluating integral of func over triangle in image with pixels pixVec
        static double evaluate(double (*func)(double, double), vector<vector<Pixel>> *pixVec, Triangle *triangle);
};

#endif