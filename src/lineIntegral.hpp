#ifndef line_integral_h
#define line_integral_h

#include "triangle.hpp"

/**
 * static class for evaluating line integral over a triangle
 */

class LineIntegral {
    private:
        LineIntegral();
    public:
        // static method for evaluating integral of func over triangle
        static double evaluate(double (*func)(double, double), Triangle *triangle);
        // static method for evaluating integral of func over triangle with vertices *a, *b, *c
        static double evaluate(double (*func)(double, double), Point *a, Point *b, Point *c);
        // TODO: static methods specialized for a function over an image (constant per pixel)
};

#endif