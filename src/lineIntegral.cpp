#include "lineIntegral.hpp"


LineIntegral::LineIntegral() {}

double LineIntegral::evaluate(function<double(double x, double y)> func, vector<vector<Pixel>> *pixVec, Triangle *triangle) {
    int xMin, xMax, yMin, yMax;
    // compute bounding box of triangle, in pixels
    triangle->boundingBox(xMin, xMax, yMin, yMax);
    // iterate over bounding box
    vector<Point> vertices = triangle->copyVertices();
    // preserve ccw direction
    Segment seg01(&vertices.at(0), &vertices.at(1));
    Segment seg12(&vertices.at(1), &vertices.at(2));
    Segment seg20(&vertices.at(2), &vertices.at(0));
    double integral = 0;
    for(int i = xMin; i <= xMax; i++) {
        for(int j = yMin; j <= yMax; j++) {
            Pixel p = (*pixVec).at(i).at(j); // local copy
            integral += p.lineIntegral(func, seg01) + p.lineIntegral(func, seg12) + p.lineIntegral(func, seg20);
        }
    }
    return integral;
}

double LineIntegral::evaluate(function<double(double, double)> func, vector<vector<Pixel>> *pixVec, Point *a, Point *b, Point *c) {
    // first compute bounding box
    double ax = a->getX();
    double ay = a->getY();
    double bx = b->getX();
    double by = b->getY();
    double cx = c->getX();
    double cy = c->getY();
    int xMin = round(min(ax, min(bx, cx)));
    int xMax = round(max(ax, max(bx, cx)));
    int yMin = round(min(ay, min(by, cy)));
    int yMax = round(max(ay, max(by, cy)));
    double integral = 0;
    // ensure order of points within segment is correct
    Segment ab(a, b);
    Segment bc(b, c);
    Segment ca(c, a);
    for(int i = xMin; i <= xMax; i++) {
        for(int j = yMin; j <= yMax; j++) {
            Pixel p = (*pixVec).at(i).at(j); // local copy
            integral += p.lineIntegral(func, ab) + p.lineIntegral(func, bc) + p.lineIntegral(func, ca);
        }
    }
    return integral;
}