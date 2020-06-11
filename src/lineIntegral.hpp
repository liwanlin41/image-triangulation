#ifndef line_integral_h
#define line_integral_h

#include <CImg.h>
#include <math.h>
#include "triangle.hpp"
#include "pixel.hpp"

using namespace cimg_library;

/**
 * static class for evaluating line integral over a triangle
 */

class LineIntegral {
    private:
        LineIntegral();
    public:
        // static method for evaluating integral of func over triangle in image img
        static double evaluate(double (*func)(double, double), Triangle *triangle);
        // static double evaluate(double (*func)(double, double), CImg<unsigned char> *img, Triangle *triangle);

        // static method for evaluating integral of func over triangle with vertices *a, *b, *c in image img
        static double evaluate(double (*func)(double, double), Point *a, Point *b, Point *c);
        // static double evaluate(double (*func)(double, double), CImg<unsigned char> *img, Point *a, Point *b, Point *c);
};

#endif