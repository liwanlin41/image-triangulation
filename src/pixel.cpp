#include <assert.h>
#include "pixel.hpp"

// helper functions

// compute (unsigned) area of the polygon enclosed by points,
// where edges of the polygon are given by points.at(i) -- points.at(i+1)
double shoelace(vector<Point> &points) {
    double area = 0;
    int n = points.size();
    for(int i = 0; i < n; i++) {
        double x0 = points.at(i).getX();
        double y0 = points.at(i).getY();
        double x1 = points.at((i+1) % n).getX();
        double y1 = points.at((i+1) % n).getY();
        area += (x1 * y0 - x0 * y1);
    }
    return area;
}

Pixel::Pixel(int x_, int y_, int c) : x(x_), y(y_), color(c) {
    // confirmed: this is cast correctly 
    corners.push_back(Point(x-0.5,y-0.5));
    corners.push_back(Point(x+0.5,y-0.5));
    corners.push_back(Point(x+0.5,y+0.5));
    corners.push_back(Point(x-0.5,y+0.5));
}

int Pixel::getColor() {
    return color;
}

bool Pixel::containsPoint(Point &p) {
    double px = p.getX();
    double py = p.getY();
    return (-0.5 + x <= px && px <= 0.5 + x) && (-0.5 + y <= py && py <= 0.5 + y);
}

double Pixel::intersectionLength(Segment &e, double *x, double *y) {
    vector<Point> intersections; // hold intersections
    for(int i = 0; i < 4; i++) {
        // retrieve a side of the pixel; at most two will have an intersection
        // unless intersections is at corners
        Segment side(&corners.at(i), &corners.at((i+1) % 4));
        if (side.intersects(e)) {
            Point intersectionPoint = side.getIntersection(e);
            bool isNewPoint = true; // whether this intersection is a new distinct point
            for(Point pt : intersections) {
                if (pt == intersectionPoint) {
                    isNewPoint = false;
                }
            }
            if (isNewPoint) {
                intersections.push_back(intersectionPoint);
            }
        }
    }
    if (intersections.size() < 2) {
        Point start = e.getStart(); // creates a copy
        Point end = e.getEnd();
        if (containsPoint(start)) {
            intersections.push_back(start);
        }
        if (containsPoint(end)) {
            intersections.push_back(end);
        }
    }
    if (intersections.size() < 2) {
        return 0;
    }
    assert(intersections.size() == 2);
    Segment contained(&intersections.at(0), &intersections.at(1));
    // check for null pointers, assign midpoint coords
    if (x && y) {
        *x = (intersections.at(0).getX() + intersections.at(1).getX())/2;
        *y = (intersections.at(0).getY() + intersections.at(1).getY())/2;
    }
    return contained.length();
}

double Pixel::lineIntegral(double (*func)(double, double), Segment &e) {
    double midX, midY;
    double length = intersectionLength(e, &midX, &midY);
    if (length == 0) {
        return 0;
    }
    return (*func)(midX, midY) * length * color;
}

double Pixel::intersectionArea(Triangle &t, double *x, double *y) {

}

double Pixel::doubleIntegral(double (*func)(double, double), Triangle &t) {
    double avgX, avgY;
    double area = intersectionArea(t, &avgX, &avgY);
    if (area == 0) {
        return 0;
    }
    return (*func)(avgX, avgY) * area * color;
}