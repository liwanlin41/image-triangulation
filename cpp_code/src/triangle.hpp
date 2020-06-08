#ifndef triangle_h
#define triangle_h

#include <vector>
#include "point.hpp"
#include "matrix.hpp"

using namespace std;

/**
 * represent a triangle whose vertices are movable
 */

class Triangle {
    private:
        Point a, b, c;
        vector<Point> vertices;
    public:
        Triangle(Point &a, Point &b, Point &c);
        double getArea();
        // get signed area based on the order a, b, c
        // assumed to be counterclockwise
        double getSignedArea();
};

#endif