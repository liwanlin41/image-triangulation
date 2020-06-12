#ifndef line_integral_h
#define line_integral_h

#include <math.h>
#include <vector>
#include <functional>
#include "triangle.hpp"
#include "pixel.hpp"

using namespace std;

/**
 * static class for evaluating line integral over a triangle
 */

class LineIntegral {
    private:
        LineIntegral();
    public:
        // static method for evaluating integral of func over triangle in image with pixels pixVec
        static double evaluate(function<double(double, double)> func, vector<vector<Pixel>> *pixVec, Triangle *triangle);

        // static method for evaluating integral of func over triangle with vertices *a, *b, *c in image with pixels pixVec
        static double evaluate(function<double(double, double)> func, vector<vector<Pixel>> *pixVec, Point *a, Point *b, Point *c);
};

#endif