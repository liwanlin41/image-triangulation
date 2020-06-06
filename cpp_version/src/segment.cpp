#include <iostream>
#include <math.h>
#include "../headers/segment.hpp"

using namespace std;

Segment::Segment(Point &a, Point &b) : endpoint1(a), endpoint2(b) {}

double Segment::length() {
    double x1 = endpoint1.getX();
    double y1 = endpoint1.getY();
    double x2 = endpoint2.getX();
    double y2 = endpoint2.getY();
    return pow(pow(x1-x2,2) + pow(y1-y2,2),0.5);
}

bool Segment::intersects(Segment other) {
    // parametrize and represent as matrix equation to be solved: t1 * x0 + (1-t1) * x1 = t2 * x2 + (1-t2) * x3
    // (x0-x1) * t1 + (x3-x2) * t2 = x3 - x1
    double matrix[2][2];
    // first column (matches t1)
    matrix[0][0] = endpoint1.getX() - endpoint2.getX();
    matrix[1][0] = endpoint1.getY() - endpoint2.getY();
    // second column (matches t2)
    matrix[0][1] = other.endpoint2.getX() - other.endpoint1.getX();
    matrix[1][1] = other.endpoint2.getY() - other.endpoint1.getY();

    double determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    double adjugate[2][2];
    adjugate[0][0] = matrix[1][1];
    adjugate[1][1] = matrix[0][0];
    adjugate[0][1] = -matrix[0][1];
    adjugate[1][0] = -matrix[1][0];

    // target vector
    double target[2] = {other.endpoint2.getX() - endpoint2.getX(), other.endpoint2.getY() - endpoint2.getY()};

    // compute scaled solution t1, t2
    double t1 = adjugate[0][0] * target[0] + adjugate[0][1] * target[1];
    double t2 = adjugate[1][0] * target[0] + adjugate[1][1] * target[1];

    if (determinant > 0) {
        return (0 <= t1 && t1 <= determinant) && (0 <= t2 && t2 <= determinant);
    } 
    // if determinant is 0 (parallel) and all points are collinear, t1 and t2 will be 0 and this line is true;
    // if parallel then t1, t2 are not both zero and this line is false (no intersection)
    // therefore this also handles the singular case
    return (determinant <= t1 && t1 <= 0) && (determinant <= t2 && t2 <= 0);
}

// only support non-parallel segments for now
// TODO: refactor with intersects
// TODO: precision of result?
Point Segment::getIntersection(Segment other) {
    double matrix[2][2];
    // first column (matches t1)
    matrix[0][0] = endpoint1.getX() - endpoint2.getX();
    matrix[1][0] = endpoint1.getY() - endpoint2.getY();
    // second column (matches t2)
    matrix[0][1] = other.endpoint2.getX() - other.endpoint1.getX();
    matrix[1][1] = other.endpoint2.getY() - other.endpoint1.getY();

    double determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    double adjugate[2][2];
    adjugate[0][0] = matrix[1][1];
    adjugate[1][1] = matrix[0][0];
    adjugate[0][1] = -matrix[0][1];
    adjugate[1][0] = -matrix[1][0];

    // target vector
    double target[2] = {other.endpoint2.getX() - endpoint2.getX(), other.endpoint2.getY() - endpoint2.getY()};

    // compute scaled solution t1, t2
    double t1 = adjugate[0][0] * target[0] + adjugate[0][1] * target[1];
    double t2 = adjugate[1][0] * target[0] + adjugate[1][1] * target[1];

    double x = (endpoint1.getX() * t1 + endpoint2.getX() * (determinant-t1)) / determinant;
    double y = (endpoint1.getY() * t1 + endpoint2.getY() * (determinant-t1)) / determinant;

    return Point(x,y);
}