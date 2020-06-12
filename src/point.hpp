#ifndef point_h
#define point_h

#include <math.h>
#include <iostream>
/**
 * represent a point (x,y) on the plane
 * mutable; allows perturbing vertices
 */

using namespace std;

class Point {
    private:
        double x, y;
    public:
        Point(double x, double y);
        double getX();
        double getY();
        double distance(Point &other);
        void move(double deltaX, double deltaY);
        bool operator==(const Point &other) const;
        bool operator!=(const Point &other) const;
        friend ostream& operator<<(ostream &os, const Point &p);
};
#endif