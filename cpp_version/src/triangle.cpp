#include "../headers/triangle.hpp"

using namespace std;

Triangle::Triangle(Point &a_, Point &b_, Point &c_) : a(a_), b(b_), c(c_) {
    vertices.push_back(a);
    vertices.push_back(b);
    vertices.push_back(c);
}

double Triangle::getSignedArea() {
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