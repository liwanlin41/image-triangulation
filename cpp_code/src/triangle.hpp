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
};

#endif