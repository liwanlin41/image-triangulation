#include "constant.cuh"

const double TOLERANCE = 1e-10;

ConstantApprox::ConstantApprox(CImg<unsigned char> *img, vector<Point> *pts, vector<array<int, 3>> &inds, double step, double ds_) 
: stepSize(step), ds(ds_) {
	// create pixel array representation
	maxX = img->width();
	maxY = img->height();
	cout << "image is " << maxX << "x" << maxY << endl;
	// allocate shared space for pixel array
	cudaMallocManaged(&pixArr, maxX * maxY * sizeof(Pixel));
	bool isGrayscale = (img->spectrum() == 1);
	for(int x = 0; x < maxX; x++) {
		for(int y = 0; y < maxY; y++) {
			int ind = x * maxY + y; // 1D pixel index
			if(isGrayscale) {
				pixArr[ind] = Pixel(x, y, (*img)(x, y));
			} else {
				int rgb[3];
				for(int i = 0; i < 3; i++) {
					rgb[i] = (*img)(x, y, 0, i);
				}
				int r = (*img)(x, y, 0, 0);
				pixArr[ind] = Pixel(x, y, rgb[0], rgb[1], rgb[2]);
			}
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
	cudaMallocManaged(&grays, numTri * sizeof(double));
	/*
	cudaMallocManaged(&reds, numTri * sizeof(double));
	cudaMallocManaged(&greens, numTri * sizeof(double));
	cudaMallocManaged(&blues, numTri * sizeof(double));
	*/

	double maxLength = 0; // get maximum side length of a triangle for space allocation
	for(int i = 0; i < numTri; i++) {
		array<int, 3> t = inds.at(i); // vertex indices for this triangle
		// constructor takes point addresses
		triArr[i] = Triangle(points + t.at(0), points + t.at(1), points + t.at(2));
		maxLength = max(maxLength, triArr[i].maxLength());
	}
	imageInt = new double[numTri];
	// create shared space for triangle iterations
	cudaMallocManaged(&workingTriangle, 3 * sizeof(Point));

	// create space for results, one slot per gpu worker
	long long maxDivisions = (int) (maxLength/ds + 1); // max num samples per side, rounded up
	// maximum possible number of samples per triangle is loosely upper bounded by 2 * maxDivisions^2
	// assumming edge lengths are bounded above by maxDivisions * 2
	long long resultSlots = max(2 * maxDivisions * maxDivisions, (long long) maxX * maxY); // at least num pixels
	cudaMallocManaged(&results0, resultSlots * sizeof(double));
	cudaMallocManaged(&results1, resultSlots * sizeof(double));
	// the above may cause errors because so much memory is required
	cudaError_t error = cudaGetLastError();
  	if(error != cudaSuccess)
  	{
    	// print the CUDA error message and exit
		printf("CUDA error: %s\n", cudaGetErrorString(error));
		cout << "An approximation of this quality is not possible due to memory limitations." << endl;
    	exit(-1);
  	}

	// create an initial approximation based on this triangulation
	integrator.initialize(pixArr, triArr, maxX, maxY, APPROXTYPE, resultSlots);
	updateApprox();
}

ConstantApprox::~ConstantApprox() {
	cudaFree(pixArr);
	cudaFree(results0);
	cudaFree(results1);
	cudaFree(points);
	cudaFree(triArr);
	cudaFree(grays);
	/*
	cudaFree(reds);
	cudaFree(greens);
	cudaFree(blues);
	*/
	cudaFree(workingTriangle);
	delete[] imageInt;
}

double ConstantApprox::computeEnergy() {
	//return integrator.constantEnergyEval(grays, numTri);
	return integrator.constantEnergyApprox(grays, numTri, ds);
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
			//boundaryChange[i] = lineIntEval(APPROXTYPE, pixArr, maxX, maxY, triArr, t, movingPt, (i == 0), results0, results1, workingTriangle);
			boundaryChange[i] = lineIntApprox(APPROXTYPE, pixArr, maxY, triArr, t, movingPt, (i == 0), results0, results1, ds, workingTriangle);
		}
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
		// compute image dA and store it for reference on next iteration
		//double val = doubleIntEval(APPROXTYPE, pixArr, maxX, maxY, triArr, t, results0);
		double val = doubleIntApprox(APPROXTYPE, pixArr, maxY, triArr + t, results0, results1, ds, workingTriangle);
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
		// get rgb values
		/*
		reds[t] = doubleIntEval(APPROXTYPE, pixArr, maxX, maxY, triArr, t, results, red) / area;
		greens[t] = doubleIntEval(APPROXTYPE, pixArr, maxX, maxY, triArr, t, results, green) / area;
		blues[t] = doubleIntEval(APPROXTYPE, pixArr, maxX, maxY, triArr, t, results, blue) / area;
		*/
	}
}

double ConstantApprox::step(double &prevEnergy, double &newEnergy) {
	double usedStep;
	computeGrad();
    while(!gradUpdate()) {
        undo(); // keep halving stepSize until no triangle is inverted
	}
	updateApprox();
    prevEnergy = newEnergy;
	newEnergy = computeEnergy();
    // TODO: tune this
	if(newEnergy > prevEnergy) { // overshot optimum?
        do {
            undo();
        } while (!gradUpdate());
        updateApprox();
		newEnergy = computeEnergy();
		usedStep = stepSize;
    } else {
		usedStep = stepSize;
		stepSize *= 2; // prevent complete vanishing to zero
	}
    cout << "new energy: " << newEnergy << endl;
	cout << "Step size: " << usedStep << endl;
	return usedStep;
}

void ConstantApprox::run(int maxIter, double eps) {
	// track change in energy for stopping point
	double newEnergy = computeEnergy();
	// initialize to something higher than newEnergy
	double prevEnergy = newEnergy + 100 * eps;
	int iterCount = 0;
	while(iterCount < maxIter && abs(prevEnergy - newEnergy) > eps) {
		cout << "iteration " << iterCount << endl;
		step(prevEnergy, newEnergy);
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
		int scale = 255;
		double area = triArr[t].getArea();
		/*
		double r = doubleIntEval(APPROXTYPE, pixArr, maxX, maxY, triArr, t, results0, RED)/(scale * area);
		double g = doubleIntEval(APPROXTYPE, pixArr, maxX, maxY, triArr, t, results0, GREEN)/(scale * area);
		double b = doubleIntEval(APPROXTYPE, pixArr, maxX, maxY, triArr, t, results0, BLUE)/(scale * area);
		*/
		double r = doubleIntApprox(APPROXTYPE, pixArr, maxY, triArr + t, results0, results1, ds, workingTriangle, RED) / (scale * area);
		double g = doubleIntApprox(APPROXTYPE, pixArr, maxY, triArr + t, results0, results1, ds, workingTriangle, GREEN) / (scale * area);
		double b = doubleIntApprox(APPROXTYPE, pixArr, maxY, triArr + t, results0, results1, ds, workingTriangle, BLUE) / (scale * area);
		fullColors.push_back({r, g, b});
		/*
		if(max(r, max(g, b)) > 1.05) {
			cout << triArr[t];
			double rp = doubleIntApprox(APPROXTYPE, pixArr, maxY, triArr + t, results, ds, workingTriangle, RED);
			double gp = doubleIntApprox(APPROXTYPE, pixArr, maxY, triArr + t, results, ds, workingTriangle, GREEN);
			double bp = doubleIntApprox(APPROXTYPE, pixArr, maxY, triArr + t, results, ds, workingTriangle, BLUE);
			printf("values %f, %f, %f\n", rp, gp, bp); 
			printf("values obtained were %f, %f, %f\n", r * scale * area, g * scale * area, b * scale * area);
			printf("area is %f\n", area);
		}
		*/
	}
	return fullColors;
}
