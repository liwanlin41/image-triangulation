#include "../headers/triangle.hpp"

using namespace std;

Triangle::Triangle(Point &a_, Point &b_, Point &c_) : a(a_), b(b_), c(c_) {
    vertices.push_back(a);
    vertices.push_back(b);
    vertices.push_back(c);
}

double Triangle::getSignedArea() {
    
}

double Triangle::getArea() {

}