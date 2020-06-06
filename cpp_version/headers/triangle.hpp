#ifndef triangle_h
#define triangle_h

#include "point.hpp"

/**
 * represent a triangle whose vertices are movable
 */

class Triangle {
    private:
        Point a, b, c;
    public:
        Triangle(Point &a, Point &b, Point &c);
        Point* getVertices();
        double getArea();
        // get signed area based on the order a, b, c
        double getSignedArea();
};

#endif