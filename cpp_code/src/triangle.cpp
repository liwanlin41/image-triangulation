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
    double ax = vertices.at(0)->getX();
    double ay = vertices.at(0)->getY();
    double bx = vertices.at(1)->getX(); 
    double by = vertices.at(1)->getY();
    double cx = vertices.at(2)->getX();
    double cy = vertices.at(2)->getY();
    Matrix matrix(bx - ax, cx - ax, by - ay, cy - ay);
    return matrix.determinant()/2;
}

double Triangle::getArea() {
    double signedArea = getSignedArea();
    if (signedArea < 0) {
        return -signedArea;
    }
    return signedArea;
}
