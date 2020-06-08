#include "triangle.hpp"

using namespace std;

Triangle::Triangle(Point *a, Point *b, Point *c) {
    vertices.push_back(a);
    vertices.push_back(b);
    vertices.push_back(c);
    if(getSignedArea() < 0) {
        vertices.clear();
    }
    // reverse direction
    vertices.push_back(a);
    vertices.push_back(c);
    vertices.push_back(b);
}

double Triangle::getSignedArea() {
    Point a = *(vertices.at(0));
    Point b = *(vertices.at(1));
    Point c = *(vertices.at(2));
    Matrix matrix(b.getX() - a.getX(), c.getX() - a.getX(), b.getY() - a.getY(), c.getY() - a.getY());
    return matrix.determinant()/2;
}

double Triangle::getArea() {
    double signedArea = getSignedArea();
    if (signedArea < 0) {
        return -signedArea;
    }
    return signedArea;
}
