#include "triangle.hpp"

using namespace std;

// custom rounding function to support needed pixel rounding
int pixelRound(double x) {
    int floor = (int) x;
    if (abs(x - floor) <= 0.5) {
        return floor;
    }
    else if (x > 0) {
        return floor + 1;
    }
    return floor - 1;
}

Triangle::Triangle(Point *a, Point *b, Point *c) {
    vertices.push_back(a);
    vertices.push_back(b);
    vertices.push_back(c);
    if(getSignedArea() < 0) {
        vertices.clear();
        // reverse direction
        vertices.push_back(a);
        vertices.push_back(c);
        vertices.push_back(b);
    }
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

double Triangle::dA(Point *p, double vx, double vy) {
    // first extract the other two endpoints; note order matters
    int index;
    for(int i = 0; i < 3; i++) {
        if(vertices.at(i) == p) {
            index = i;
        }
    }
    Point *edgePoints[2];
    // retrieve in ccw order
    edgePoints[0] = vertices.at((index+1) % 3);
    edgePoints[1] = vertices.at((index+2) % 3);
    // change is -velocity dot edge normal of length |e|/2
    Segment opposite(edgePoints[0], edgePoints[1]);
    Matrix velocity(vx, vy);
    Matrix norm = opposite.scaledNormal();
    // 1 by 1
    Matrix grad = velocity.transpose().multiply(norm);
    return -grad.get(0,0);
}

double Triangle::gradX(Point *p) {
    return dA(p, 1, 0);
}

double Triangle::gradY(Point *p) {
    return dA(p, 0, 1);
}

void Triangle::boundingBox(int &xMin, int &xMax, int &yMin, int &yMax) {
    double ax = vertices.at(0)->getX();
    double ay = vertices.at(0)->getY();
    double bx = vertices.at(1)->getX();
    double by = vertices.at(1)->getY();
    double cx = vertices.at(2)->getX();
    double cy = vertices.at(2)->getY();
    // round to nearest to get pixel coords
    xMin = pixelRound(min(ax, min(bx, cx)));
    xMax = pixelRound(max(ax, max(bx, cx)));
    yMin = pixelRound(min(ay, min(by, cy)));
    yMax = pixelRound(max(ay, max(by, cy)));
}

vector<Point> Triangle::copyVertices() {
    vector<Point> copy;
    for(Point *ptr : vertices) {
        Point copiedPoint = *ptr;
        copy.push_back(copiedPoint);
    }
    return copy;
}

double Triangle::getSignedArea(Point *a, Point *b, Point *c) {
    double ax = a->getX();
    double ay = a->getY();
    double bx = b->getX(); 
    double by = b->getY();
    double cx = c->getX();
    double cy = c->getY();
    Matrix matrix(bx - ax, cx - ax, by - ay, cy - ay);
    return matrix.determinant()/2;
}

bool Triangle::contains(Point &p) {
    // p is inside the triangle iff the orientations of the triangles
    // with vertices (vertices[i], vertices[i+1], p) are all the same
    bool signs[3]; // hold signs of triangles (true if ccw)
    for(int i = 0; i < 3; i++) {
        signs[i] = (Triangle::getSignedArea(vertices.at(i), vertices.at((i+1) % 3), &p) >= 0);
    }
    return (signs[0] == signs[1] && signs[1] == signs[2]);
}

ostream& operator<<(ostream& os, const Triangle &t) {
    os << "Triangle ";
    for(Point *ptr : t.vertices) {
        os << *ptr << " ";
    }
    os << "\n";
    return os;
}