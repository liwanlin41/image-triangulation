#include "linear.h"

const double LinearApprox::matrix[3][3] = {{9,-3,-3}, {-3,9,-3}, {-3,-3,9}};

LinearApprox::LinearApprox(CImg<unsigned char> *img, double step, double ds) : Approx(img, step, ds) {
}

LinearApprox::~LinearApprox() {
    // TODO
    for(int i = 0; i < numTri; i++) {
        delete[] coefficients[i];
        delete[] basisIntegral[i];
    }
    delete[] coefficients;
    delete[] basisIntegral;
}

ApproxType LinearApprox::getApproxType() {
    return APPROXTYPE;
}

void LinearApprox::reallocateSpace(int oldNumTri) {
    for(int i = 0; i < oldNumTri; i++) {
        delete[] coefficients[i];
        delete[] basisIntegral[i];
    }
    delete[] coefficients;
    delete[] basisIntegral;
    coefficients = new double *[numTri];
    basisIntegral = new double *[numTri];
    for(int i = 0; i < numTri; i++) {
        coefficients[i] = new double[3];
        basisIntegral[i] = new double[3];
    }
}

void LinearApprox::initialize(int pixelRate) {
    Approx::initialize(APPROXTYPE, pixelRate);
    coefficients = new double *[numTri];
    basisIntegral = new double *[numTri];
    for(int i = 0; i < numTri; i++) {
        coefficients[i] = new double[3];
        basisIntegral[i] = new double[3];
    }
    updateApprox();
}

void LinearApprox::initialize(vector<Point> &pts, vector<array<int, 3>> &faces) {
    Approx::initialize(APPROXTYPE, pts, faces);
    coefficients = new double *[numTri];
    basisIntegral = new double *[numTri];
    for(int i = 0; i < numTri; i++) {
        coefficients[i] = new double[3];
        basisIntegral[i] = new double[3];
    }
    updateApprox();
}

double LinearApprox::computeEnergy() {

}

void LinearApprox::computeGrad() {

}

void LinearApprox::updateApprox() {

}

void LinearApprox::computeEdgeEnergies(vector<array<double, 3>> *edgeEnergies) {

}

vector<array<double, 3>> LinearApprox::getColors() {
    /*
    vector<array<double, 3>> colors;
    colors.push_back({1,0,0});
    colors.push_back({1,0,0});
    colors.push_back({0.5,0,0});
    
    colors.push_back({0,1,0});
    colors.push_back({1,1,1});
    colors.push_back({0,1,0});
    
    colors.push_back({0,0,1});
    colors.push_back({0,0,1});
    colors.push_back({0,0,0});

    colors.push_back({0,0,0});
    colors.push_back({0,0,0});
    colors.push_back({1,1,1});
    return colors;
    */
}

/*
vector<Point> LinearApprox::testPoints() {
    vector<Point> res;
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 2; j++) {
            res.push_back(Point(i, j));
        }
    }
    return res;
}

vector<array<int, 3>> LinearApprox::testFaces() {
    vector<array<int, 3>> faces;
    faces.push_back({0,1,2});
    faces.push_back({1,2,3});
    faces.push_back({2,3,4});
    faces.push_back({3,4,5});
    return faces;
}
*/