#include "constant.h"

const double TOLERANCE = 1e-10;

ConstantApprox::ConstantApprox(CImg<unsigned char> *img, double step, double ds) : Approx(img, step, ds) {
}

void ConstantApprox::reallocateSpace() {
	delete[] imageInt;
	delete[] grays;
	imageInt = new double[numTri];
	grays = new double[numTri];
}

void ConstantApprox::initialize(vector<Point> &pts, vector<array<int, 3>> &inds) {
	Approx::initialize(APPROXTYPE, pts, inds); // call parent method
	imageInt = new double[numTri];
	grays = new double[numTri];

	// create an initial approximation based on this triangulation
	updateApprox();
}

void ConstantApprox::initialize(int pixelRate) {
	Approx::initialize(APPROXTYPE, pixelRate);
	imageInt = new double[numTri];
	grays = new double[numTri];
	// create an initial approximation based on this triangulation
	updateApprox();
}

ConstantApprox::~ConstantApprox() {
	delete[] grays;
	delete[] imageInt;
}

double ConstantApprox::computeEnergy() {
	double totalEnergy = 0;
    for(int t = 0; t < numTri; t++) {
        totalEnergy += integrator.constantEnergyEval(triArr + t, grays[t], ds);
    }
    return totalEnergy;
}

void ConstantApprox::computeGrad() {
	// clear gradients from last iteration
	for(int i = 0; i < numPoints; i++) {
		gradX[points + i] = 0;
		gradY[points + i] = 0;
	}
	for(int i = 0; i < numTri; i++) {
		// integral of fdA, retrieved from last updateApprox iteration
		double imageIntegral = imageInt[i];
		for(int j = 0; j < 3; j++) {
			double changeX, changeY;
			gradient(i, j, imageIntegral, &changeX, &changeY);
			// constrain points on boundary of image
			if(triArr[i].vertices[j]->isBorderX()) {
				changeX = 0;
			}
			if(triArr[i].vertices[j]->isBorderY()) {
				changeY = 0;
			}
			gradX[triArr[i].vertices[j]] += changeX;
			gradY[triArr[i].vertices[j]] += changeY;
		}
	}
}

void ConstantApprox::gradient(int t, int movingPt, double imageIntegral, double *gradX, double *gradY) {
	// to save time, only compute integrals if triangle is non-degenerate;
	// degenerate triangle has 0 energy and is locally optimal, set gradient to 0
	double area = triArr[t].getArea();
	double gradient[2] = {0, 0};
	if (area > TOLERANCE) {
		double dA[2] = {triArr[t].gradX(movingPt), triArr[t].gradY(movingPt)};
		double boundaryChange[2];
		// compute gradient in x and y direction
		for(int i = 0; i < 2; i++) {
			// sample more frequently because both time and space allow (or don't)
			boundaryChange[i] = integrator.lineIntEval(triArr+t, movingPt, (i == 0), ds);
		}
		for(int j = 0; j < 2; j++) {
			gradient[j] = (2 * area * imageIntegral * boundaryChange[j]
				- imageIntegral * imageIntegral * dA[j]) / (-area * area)
				- 100 * dA[j] / area; // add in log barrier gradient
		}
	}
	// check for null pointers
	if (gradX && gradY) {
		*gradX = gradient[0];
		*gradY = gradient[1];
	}
}

void ConstantApprox::updateApprox() {
	for(int t = 0; t < numTri; t++) {
		// compute image dA and store it for reference on next iteration
		double val = integrator.doubleIntEval(triArr+t, ds);
		imageInt[t] = val;
		double area = triArr[t].getArea();
		// take average value
		double approxVal = val / area;
		// handle degeneracy
		if (isnan(approxVal)) {
			assert(area < TOLERANCE);
			approxVal = 255; // TODO: something better than this
		}
		grays[t] = min(255.0, approxVal); // prevent blowup in case of poor approximation
	}
}

void ConstantApprox::computeEdgeEnergies(vector<array<double, 3>> *edgeEnergies) {
	for(auto ii = edgeBelonging.begin(); ii != edgeBelonging.end(); ii++) {
		array<int, 2> edge = ii->first;
		vector<int> triangles = ii->second; // triangles containing edge
		// compute current total energy over these triangles
		double curEnergy = 0;
		for(int t : triangles) {
			curEnergy += integrator.constantEnergyEval(triArr+t, grays[t], ds);
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
			double color1 = integrator.doubleIntEval(&t1, ds) / area;
			double color2 = integrator.doubleIntEval(&t2, ds) / area;
			newEnergy += integrator.constantEnergyEval(&t1, color1, ds) + integrator.constantEnergyEval(&t2, color2, ds);
		}
		// change in energy due to subdivision
		edgeEnergies->push_back({(double) edge[0], (double) edge[1], newEnergy - curEnergy});
	}
}

vector<array<double,3>> ConstantApprox::getColors() {
	vector<array<double, 3>> fullColors;
	for(int t = 0; t < numTri; t++) {
		// scale to fit polyscope colors TODO: check that this is correct
		int scale = 255;
		double area = triArr[t].getArea();
		double r = integrator.doubleIntEval(triArr+t, ds, RED) / (scale * area);
		double g = integrator.doubleIntEval(triArr+t, ds, GREEN) / (scale * area);
		double b = integrator.doubleIntEval(triArr+t, ds, BLUE) / (scale * area);
		fullColors.push_back({r, g, b});
	}
	return fullColors;
}