#ifndef point_h
#define point_h

#include <math.h>
#include <iostream>
#include <iomanip>
/**
 * represent a point (x,y) on the plane
 * mutable; allows perturbing vertices
 */

using namespace std;

class Point {
    private:
        double x, y;
        // determine if point is on edge of image and should not be moved
        bool borderX, borderY;
    public:
        Point(double x, double y, bool borderX = false, bool borderY = false);
        double getX() const;
        double getY() const;
        // return true if point was constructed on a vertical image edge
        bool isBorderX() const;
        // return true if point was constructed on a horizontal image edge
        bool isBorderY() const;
        double distance(Point &other);
        void move(double deltaX, double deltaY);
        bool operator==(const Point &other) const;
        bool operator!=(const Point &other) const;
        friend ostream& operator<<(ostream &os, const Point &p);
};
#endif