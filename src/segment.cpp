#include <iostream>
#include <math.h>
#include "segment.hpp"
#include <algorithm>

using namespace std;

// helper function for determining if t in [a,b] where order of a,b is unknown
bool isBetween(double t, double a, double b) {
    return (a <= t && t <= b) || (b <= t && t <= a);
}

// helper function for determining whether there exists
// an axis-parallel separator between {a, b} and {c, d}
bool areDisjoint(Point *a, Point *b, Point *c, Point *d) {
    // check x direction
    bool x_separable = (min(a->getX(), b->getX()) > max(c->getX(), d->getX())) || 
        (min(c->getX(), d->getX()) > max(a->getX(), b->getX()));
    bool y_separable = (min(a->getY(), b->getY()) > max(c->getY(), d->getY())) ||
        (min(c->getY(), d->getY()) > max(a->getY(), b->getY()));
    return x_separable || y_separable;
}

Segment::Segment(Point *a, Point *b) : endpoint1(a), endpoint2(b) {}

double Segment::length() {
    double x1 = endpoint1->getX();
    double y1 = endpoint1->getY();
    double x2 = endpoint2->getX();
    double y2 = endpoint2->getY();
    return pow(pow(x1-x2,2) + pow(y1-y2,2),0.5);
}

bool Segment::intersects(Segment other) {
    // parametrize and represent as matrix equation to be solved: t1 * x0 + (1-t1) * x1 = t2 * x2 + (1-t2) * x3
    // (x0-x1) * t1 + (x3-x2) * t2 = x3 - x1
    double matrix[2][2];
    // first column (matches t1)
    matrix[0][0] = endpoint1->getX() - endpoint2->getX();
    matrix[1][0] = endpoint1->getY() - endpoint2->getY();
    // second column (matches t2)
    matrix[0][1] = other.endpoint2->getX() - other.endpoint1->getX();
    matrix[1][1] = other.endpoint2->getY() - other.endpoint1->getY();

    double determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    double adjugate[2][2];
    adjugate[0][0] = matrix[1][1];
    adjugate[1][1] = matrix[0][0];
    adjugate[0][1] = -matrix[0][1];
    adjugate[1][0] = -matrix[1][0];

    // target vector
    double target[2] = {other.endpoint2->getX() - endpoint2->getX(), other.endpoint2->getY() - endpoint2->getY()};

    // compute scaled solution t1, t2
    double t1 = adjugate[0][0] * target[0] + adjugate[0][1] * target[1];
    double t2 = adjugate[1][0] * target[0] + adjugate[1][1] * target[1];

    // first handle degenerate case; checking t1 = t2 = 0 will give the correct result for the parallel case
    // but an additional check needs to be done for the collinear case
    // EDIT: return false for segments intersecting in !=1 point
    if (determinant == 0) {
        return false;
        // return (t1 == 0 && t2 == 0) && !areDisjoint(endpoint1, endpoint2, other.endpoint1, other.endpoint2); 
    }

    return isBetween(t1, 0, determinant) && isBetween(t2, 0, determinant);
}

// only support non-parallel segments for now
// TODO: refactor with intersects
// TODO: precision of result?
Point Segment::getIntersection(Segment other) {
    double matrix[2][2];
    // first column (matches t1)
    matrix[0][0] = endpoint1->getX() - endpoint2->getX();
    matrix[1][0] = endpoint1->getY() - endpoint2->getY();
    // second column (matches t2)
    matrix[0][1] = other.endpoint2->getX() - other.endpoint1->getX();
    matrix[1][1] = other.endpoint2->getY() - other.endpoint1->getY();

    double determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    double adjugate[2][2];
    adjugate[0][0] = matrix[1][1];
    adjugate[1][1] = matrix[0][0];
    adjugate[0][1] = -matrix[0][1];
    adjugate[1][0] = -matrix[1][0];

    // target vector
    double target[2] = {other.endpoint2->getX() - endpoint2->getX(), other.endpoint2->getY() - endpoint2->getY()};

    // compute scaled solution t1, t2
    double t1 = adjugate[0][0] * target[0] + adjugate[0][1] * target[1];
    double t2 = adjugate[1][0] * target[0] + adjugate[1][1] * target[1];

    double x = (endpoint1->getX() * t1 + endpoint2->getX() * (determinant-t1)) / determinant;
    double y = (endpoint1->getY() * t1 + endpoint2->getY() * (determinant-t1)) / determinant;

    return Point(x,y);
}

Matrix Segment::unitNormal() {
    double deltaX = endpoint2->getX() - endpoint1->getX();
    double deltaY = endpoint2->getY() - endpoint1->getY();
    double unitX = deltaX / length();
    double unitY = deltaY / length();
    // rotate pi/2 clockwise
    return Matrix(unitY, -unitX);
}

Matrix Segment::scaledNormal() {
    double deltaX = endpoint2->getX() - endpoint1->getX();
    double deltaY = endpoint2->getY() - endpoint1->getY();
    return Matrix(deltaY/2, -deltaX/2);
}

ostream& Segment::operator<<(ostream& os) {
    os << "(" << endpoint1->getX() << ", " << endpoint1->getY() << ") -- (" << endpoint2->getX() << ", " << endpoint2->getY() << ")\n";
    return os;
}