#include "constant.cuh"

// constants for converting color image to grayscale
const double RED_LUMINANCE = 0.2126;
const double GREEN_LUMINANCE = 0.7152;
const double BLUE_LUMINANCE = 0.0722;

const double TOLERANCE = 1e-10;

// get the luminance at pixel (x, y) by standard luminance transformation
int getLuminance(CImg<unsigned char> *img, int x, int y) {
	int red_val = (int) (*img)(x, y, 0, 0);
	int green_val = (int) (*img)(x, y, 0, 1);
	int blue_val = (int) (*img)(x, y, 0, 2);
	return round(red_val * RED_LUMINANCE + green_val * GREEN_LUMINANCE + blue_val * BLUE_LUMINANCE);
}

// custom rounding function to support needed pixel rounding
int customRound(double x) {
	int floor = (int) x;
	if (abs(x - floor) <= 0.5) {
		return floor;
	} else if (x > 0) {
		return floor + 1;
	}
	return floor - 1;
}

ConstantApprox::ConstantApprox(CImg<unsigned char> *img, vector<Point> *pts, vector<array<int, 3>> &inds, double step) 
: stepSize(step) {
	// create pixel array representation
	maxX = img->width();
	maxY = img->height();
	cout << "image is " << maxX << "x" << maxY << endl;
	// allocate shared space for pixel array
	cudaMallocManaged(&pixArr, maxX * maxY * sizeof(Pixel));
	// create space for results, one slot per pixel
	cudaMallocManaged(&results, maxX * maxY * sizeof(double));
	bool isGrayscale = (img->spectrum() == 1);
	for(int x = 0; x < maxX; x++){
		for(int y = 0; y < maxY; y++) {
			int color = (isGrayscale) ? (*img)(x, y) : getLuminance(img, x, y);
			pixArr[x * maxY + y] = Pixel(x, y, color);
		}
	}

	// load in points of triangulation
	numPoints = pts->size();
	// allocate shared space for points
	cudaMallocManaged(&points, numPoints * sizeof(Point));
	// copy everything in TODO: make this more efficient (get directly from source)
	for(int i = 0; i < numPoints; i++) {
		points[i] = pts->at(i);
	}

	// now load in all the triangles
	numTri = inds.size();
	// allocate shared space for triangles and colors
	cudaMallocManaged(&triArr, numTri * sizeof(Triangle));
	cudaMallocManaged(&colors, numTri * sizeof(double));
	for(int i = 0; i < numTri; i++) {
		array<int, 3> t = inds.at(i); // vertex indices for this triangle
		// constructor takes point addresses
		triArr[i] = Triangle(points + t.at(0), points + t.at(1), points + t.at(2));
	}
	// create shared space for triangle iterations
	cudaMallocManaged(&workingTriangle, 3 * sizeof(Triangle));

	// create an initial approximation based on this triangulation
	updateApprox();
}

ConstantApprox::~ConstantApprox() {
	cudaFree(pixArr);
	cudaFree(results);
	cudaFree(points);
	cudaFree(triArr);
	cudaFree(colors);
	cudaFree(workingTriangle);
}

double ConstantApprox::computeEnergy() {
	return constantEnergyEval(pixArr, maxX, maxY, triArr, colors, numTri, results);
}

void ConstantApprox::computeGrad() {
	// clear gradients from last iteration
	for(int i = 0; i < numPoints; i++) {
		gradX[points + i] = 0;
		gradY[points + i] = 0;
	}
	for(int i = 0; i < numTri; i++) {
		// integral of fdA
		double imageIntegral = doubleIntEval(APPROXTYPE, pixArr, maxX, maxY, triArr, i, results);
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
		// compute gradient in x direction
		boundaryChange[0] = lineIntEval(APPROXTYPE, pixArr, maxX, maxY, triArr, t, movingPt, true, results, workingTriangle);
		// and in y direction
		boundaryChange[1] = lineIntEval(APPROXTYPE, pixArr, maxX, maxY, triArr, t, movingPt, false, results, workingTriangle);
		for(int j = 0; j < 2; j++) {
			gradient[j] = (2 * area * imageIntegral * boundaryChange[j]
				- imageIntegral * imageIntegral * dA[j]) / (-area * area);
		}
	}
	// check for null pointers
	if (gradX && gradY) {
		*gradX = gradient[0];
		*gradY = gradient[1];
	}
}

bool ConstantApprox::gradUpdate() {
	// gradient descent update for each point
	for(int i = 0; i < numPoints; i++) {
		points[i].move(-stepSize * gradX.at(points+i), -stepSize * gradY.at(points+i));
	}
	// check validity of result
	for(int i = 0; i < numTri; i++) {
		if(triArr[i].getSignedArea() < 0) {
			return false;
		}
	}
	return true;
}

void ConstantApprox::undo() {
	for(int i = 0; i < numPoints; i++) {
		points[i].move(stepSize * gradX.at(points+i), stepSize * gradY.at(points+i));
	}
	stepSize /= 2;
}

void ConstantApprox::updateApprox() {
	for(int t = 0; t < numTri; t++) {
		// compute image dA
		double val = doubleIntEval(APPROXTYPE, pixArr, maxX, maxY, triArr, t, results);
		// take average value
		double approxVal = val / triArr[t].getArea();
		// handle degeneracy
		if (isnan(approxVal)) {
			assert(triArr[t].getArea() < TOLERANCE);
			approxVal = 255; // TODO: something better than this
		}
		colors[t] = approxVal;
	}
}

void ConstantApprox::run(int maxIter, double eps) {
	// track change in energy for stopping point
	double newEnergy = computeEnergy();
	// initialize to something higher than newEnergy
	double prevEnergy = newEnergy + 100 * eps;
	int iterCount = 0;
	while(iterCount < maxIter && abs(prevEnergy - newEnergy) > eps) {
		cout << "iteration " << iterCount << endl;
		computeGrad();
		while(!gradUpdate()) {
			undo(); // keep halving stepSize until it works
		}
		updateApprox();
		prevEnergy = newEnergy;
		newEnergy = computeEnergy();
		if(newEnergy > prevEnergy) { // overshot optimum
			do {
				undo();
			} while (!gradUpdate());
			updateApprox();
			newEnergy = computeEnergy();
		}
		cout << "new energy: " << newEnergy << endl;
		cout << "step size: " << stepSize << endl;
		iterCount++;
	}
}

double ConstantApprox::getStep() {
	return stepSize;
}

// inefficient TODO: fix
vector<Point> ConstantApprox::getVertices() {
	vector<Point> vertices;
	for(int i = 0; i < numPoints; i++) {
		vertices.push_back(points[i]);
	}
	return vertices;
}

vector<array<double,3>> ConstantApprox::getColors() {
	vector<array<double, 3>> fullColors;
	for(int t = 0; t < numTri; t++) {
		// scale to fit polyscope colors TODO: check that this is correct
		double approxColor = colors[t]/255;
		fullColors.push_back({approxColor, approxColor, approxColor});
	}
	return fullColors;
}
