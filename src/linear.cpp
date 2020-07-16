#include "linear.h"

const double LinearApprox::matrix[3][3] = {{9,-3,-3}, {-3,9,-3}, {-3,-3,9}};

LinearApprox::LinearApprox(CImg<unsigned char> *img, double step, double ds) : Approx(img, step, ds) {
}

LinearApprox::~LinearApprox() {
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
    double totalEnergy = 0;
    for(int t = 0; t < numTri; t++) {
        array<int, 3> vertices = faces.at(t);
        //double energy1 = integrator.linearEnergyApprox(points + vertices[0], points + vertices[1], points + vertices[2], coefficients[t], ds);
        totalEnergy += integrator.linearEnergyApprox(triArr + t, coefficients[t], ds);
    }
    return totalEnergy;
}

void LinearApprox::computeGrad() {

}

void LinearApprox::updateApprox() {
    for(int t = 0; t < numTri; t++) {
		// compute image phi_i dA and store it for reference on next iteration
        for(int i = 0; i < 3; i++) {
            basisIntegral[t][i] = integrator.doubleIntEval(triArr + t, ds, GRAY, i);
        }
		double area = triArr[t].getArea();
        // extract coefficients
        for(int i = 0; i < 3; i++) {
            double coeff = 0;
            for(int j = 0; j < 3; j++) {
                coeff += matrix[i][j] * basisIntegral[t][j];
            }
            coefficients[t][i] = min(255.0, coeff / area); // prevent blowup
        }
	}
}

void LinearApprox::computeEdgeEnergies(vector<array<double, 3>> *edgeEnergies) {

}

vector<array<double, 3>> LinearApprox::getColors() {
    vector<array<double, 3>> colors;
    const int scale = 255;
    for(int t = 0; t < numTri; t++) {
        double area = triArr[t].getArea();
        // one vertex at a time and extract rgb separately
        double r[3];
        double g[3];
        double b[3];
        for(int j = 0; j < 3; j++) {
            r[j] = integrator.doubleIntEval(triArr + t, ds, RED, j);
            g[j] = integrator.doubleIntEval(triArr + t, ds, GREEN, j);
            b[j] = integrator.doubleIntEval(triArr + t, ds, BLUE, j);
        }
        // do each vertex
        for(int i = 0; i < 3; i++) {
            array<double, 3> vColor = {0,0,0};
            for(int j = 0; j < 3; j++) {
                vColor[0] += matrix[i][j] * r[j] / (scale * area);
                vColor[1] += matrix[i][j] * g[j] / (scale * area);
                vColor[2] += matrix[i][j] * b[j] / (scale * area);
            }
            colors.push_back(vColor);
        }
    }
    return colors;
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