#include "doubleIntegral.hpp"

DoubleIntegral::DoubleIntegral() {}

double DoubleIntegral::evaluate(function<double(double, double)> func, vector<vector<Pixel>> *pixVec, Triangle *triangle) {
    int xMin, xMax, yMin, yMax;
    // compute bounding box of triangle, in pixels
    triangle->boundingBox(xMin, xMax, yMin, yMax);
    // iterate over bounding box
    double integral = 0;
    for(int i = xMin; i <= xMax; i++) {
        for(int j = yMin; j <= yMax; j++) {
            Pixel p = (*pixVec).at(i).at(j); // local copy
            integral += p.doubleIntegral(func, *triangle);
        }
    }
    return integral;
}
