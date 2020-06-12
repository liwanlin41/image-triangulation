#ifndef pixel_h
#define pixel_h

#include <vector>
#include <functional>
#include "point.hpp"
#include "segment.hpp"
#include "triangle.hpp"
using namespace std;

/**
 * class representing an image pixel as a square, with pixel (x,y) being a unit square centered at (x,y)
 */

class Pixel {
    private:
        int color; // represent luminance of pixel from 0 to 255; may change to support color images later
        int x, y; // center of pixel
        vector<Point> corners; // corners of pixel in ccw order
    public:
        // create a pixel centered at (x,y) with color color
        Pixel(int x, int y, int color);
        int getColor();
        // compute the length of the intersection of Segment e
        // with this pixel (may be 0); if x, y are given, store the coordinates
        // of the midpoint of the segment in x, y
        double intersectionLength(Segment &e, double *x = NULL, double *y = NULL);
        // return true if this pixel contains point p (p is inside or on the boundary)
        bool containsPoint(Point &p);

        // return approximate line integral contribution of this pixel when multiplied
        // with integrand func (exact for func linear); in other words, integrate
        // color times func over the subsegment of e lying in this pixel
        double lineIntegral(function<double(double, double)> func, Segment &e);

        // compute area of the intersection of this pixel with Triangle t (may be 0);
        // store the boundary of the intersection in polygon if available
        double intersectionArea(Triangle &t, vector<Point> *polygon = NULL);

        // return approximate line integral contribution of this pixel when multiplied
        // with integrand func; in other words, integrate color times func over the part
        // of this pixel that lies inside triangle t;
        // approximation by average value of x, y in intersection of pixel and triangle
        double doubleIntegral(function<double(double, double)> func, Triangle &t);
};

#endif