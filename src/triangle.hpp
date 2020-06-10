#ifndef triangle_h
#define triangle_h

#include <vector>
#include "point.hpp"
#include "segment.hpp"
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
        // get the change in area when the vertex *p is moving at velocity (vx, vy)
        double dA(Point *p, double vx, double vy);
        // get the gradient in the x direction for vertex *p
        double gradX(Point *p);
        // get the gradient in the y direction for vertex *p
        double gradY(Point *p);
        // set input variables to bounding box of triangle
        void boundingBox(double &xMin, double &xMax, double &yMin, double &yMax);
};

#endif