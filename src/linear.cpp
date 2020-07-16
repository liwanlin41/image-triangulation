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

void LinearApprox::computeCoeffs(Triangle *tri, double *coeffs, ColorChannel channel) {
    double L[3]; // integral of f phi_i dA
    for(int i = 0; i < 3; i++) {
        L[i] = integrator.doubleIntEval(tri, ds, channel, i);
    }
    for(int i = 0; i < 3; i++) {
        double coeff = 0;
        for(int j = 0; j < 3; j++) {
            coeff += matrix[i][j] * L[j];
        }
        coeffs[i] = min(255.0, coeff / tri->getArea());
    }
}

void LinearApprox::computeEdgeEnergies(vector<array<double, 3>> *edgeEnergies) {
    for(auto ii = edgeBelonging.begin(); ii != edgeBelonging.end(); ii++) {
		array<int, 2> edge = ii->first;
		vector<int> triangles = ii->second; // triangles containing edge
		// compute current total energy over these triangles
		double curEnergy = 0;
		for(int t : triangles) {
			curEnergy += integrator.linearEnergyApprox(triArr+t, coefficients[t], ds);
		}
		// find new point that may be added to mesh
		Point endpoint0 = points[edge[0]];
		Point endpoint1 = points[edge[1]];
		double midX = (endpoint0.getX() + endpoint1.getX()) / 2;
		double midY = (endpoint0.getY() + endpoint1.getY()) / 2; 
		Point midpoint(midX, midY);

		double newEnergy = 0;
		for(int t : triangles) {
			Point opposite;
			// get opposite vertex
			for(int v = 0; v < 3; v++) {
				// for accuracy, use raw indices rather than point location
				if(faces.at(t).at(v) != edge[0] && faces.at(t).at(v) != edge[1]) {
					opposite = points[faces.at(t).at(v)];
				}
			}
			// the two triangles formed by cutting this edge
			Triangle t1(&midpoint, &opposite, &endpoint0);
			Triangle t2(&midpoint, &opposite, &endpoint1);
			// equal area of both triangles
			double area = triArr[t].getArea() / 2;
			// get energy on subdivided triangles
            double coeffs1[3];
            double coeffs2[3];
            double L1[3];
            double L2[3];
            // get integral f phi_i dA
            for(int i = 0; i < 3; i++) {
                L1[i] = integrator.doubleIntEval(&t1, ds, GRAY, i);
                L2[i] = integrator.doubleIntEval(&t2, ds, GRAY, i);
            }
			double color1 = integrator.doubleIntEval(&t1, ds) / area;
			double color2 = integrator.doubleIntEval(&t2, ds) / area;
			newEnergy += integrator.constantEnergyEval(&t1, color1, ds) + integrator.constantEnergyEval(&t2, color2, ds);
		}
		// change in energy due to subdivision
		edgeEnergies->push_back({(double) edge[0], (double) edge[1], newEnergy - curEnergy});
	}
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
        computeCoeffs(triArr + t, r, RED);
        computeCoeffs(triArr + t, g, GREEN);
        computeCoeffs(triArr + t, b, BLUE);
        for(int i = 0; i < 3; i++) {
            colors.push_back({r[i]/scale, g[i]/scale, b[i]/scale});
        }
    }
    return colors;
}