#include <assert.h>
#include "pixel.hpp"

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
        Segment side(&corners.at(i), &corners.at((i+1) % 4));
        if (side.intersects(e)) {
            intersections.push_back(side.getIntersection(e));
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
    if (intersections.size() == 0) {
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
    double mid_x, mid_y;
    double length = intersectionLength(e, &mid_x, &mid_y);
    if (length == 0) {
        return 0;
    }
    return (*func)(mid_x, mid_y) * length * color;
}