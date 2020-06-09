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
        vector<Point*> vertices;
    public:
        // construct from three Point pointers and orient in ccw direction
        Triangle(Point *a, Point *b, Point *c);
        double getArea();
        // get signed area based on the order of vertices
        // with ccw direction positive
        double getSignedArea();
        // get the derivative of area when the vertex *p is moving at velocity (vx, vy)
        // where at least one of vx, vy must be 0; corresponds to taking gradient in one direction
        double dA(Point *p, double vx, double vy);
};

#endif