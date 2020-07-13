#include "constant.cuh"

const double TOLERANCE = 1e-10;

ConstantApprox::ConstantApprox(CImg<unsigned char> *img, double step, double ds_) : originalStep(step), stepSize(step), ds(ds_) {
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
}

void ConstantApprox::initialize(vector<Point> *pts, vector<array<int, 3>> &inds) {
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
	faces = inds;
	for(int i = 0; i < numTri; i++) {
		array<int, 3> t = inds.at(i); // vertex indices for this triangle
		// constructor takes point addresses
		triArr[i] = Triangle(points + t.at(0), points + t.at(1), points + t.at(2));
		maxLength = max(maxLength, triArr[i].maxLength());
	}
	imageInt = new double[numTri];

	// initialize edge dictionary from faces
	for(int i = 0; i < numTri; i++) {
		for(int j = 0; j < 3; j++) {
			array<int, 2> edge = {faces.at(i)[j], faces.at(i)[((j+1)%3)]};
			// ensure edges are ordered
			if(edge[0] > edge[1]) {
				edge[0] = edge[1];
				edge[1] = faces.at(i)[j];
			}
			if(edgeBelonging.find(edge) == edgeBelonging.end()) { // does not exist
				edgeBelonging[edge] = vector<int>();
			}
			edgeBelonging[edge].push_back(i);
		}		
	}

	// initialize integrator

	// find space needed for results, one slot per gpu worker
	long long maxDivisions = (int) (maxLength/ds + 1); // max num samples per side, rounded up
	// maximum possible number of samples per triangle is loosely upper bounded by 2 * maxDivisions^2
	// assumming edge lengths are bounded above by maxDivisions * 2
	long long resultSlots = max(2 * maxDivisions * maxDivisions, (long long) maxX * maxY); // at least num pixels
	integrator.initialize(pixArr, maxX, maxY, APPROXTYPE, resultSlots);

	// create an initial approximation based on this triangulation
	updateApprox();
}

void ConstantApprox::initialize(int pixelRate) {
	// create points
	int numX = ceil(((double) maxX) / pixelRate) + 1; // number of samples in x direction
	int numY = ceil(((double) maxY) / pixelRate) + 1;
	double dx = ((double) maxX) / (numX - 1); // step size in x direction, remembering to get both endpoints
	double dy = ((double) maxY) / (numY - 1);

	// create shared space for points
	numPoints = numX * numY;
	cudaMallocManaged(&points, numPoints * sizeof(Point));

	for(int i = 0; i < numX; i++) {
		bool isBoundX = (i == 0) || (i == numX - 1); // whether point is on vertical boundary
		for(int j = 0; j < numY; j++) {
			bool isBoundY = (j == 0) || (j == numY - 1);
			int index1D = i * numY + j;
			// shift by (-0.5, -0.5) to align to edge of image (lattice points at pixel centers)
			points[index1D] = Point(i * dx - 0.5, j * dy - 0.5, isBoundX, isBoundY);
		}
	}
	cout << "starting grid: " << numX << "x" << numY << endl;

	// create triangles
	numTri = 2 * (numX - 1) * (numY - 1);
	cudaMallocManaged(&triArr, numTri * sizeof(Triangle));
	cudaMallocManaged(&grays, numTri * sizeof(double));
	imageInt = new double[numTri];

	int triInd = 0; // index the triangles
	for(int i = 0; i < numX; i++) {
		for(int j = 0; j < numY; j++) {
			int index1D = i * numY + j;
			// randomly triangulate the square with min x,y corner at this point
			if(i < numX - 1 && j < numY - 1) {
				Point *pt = points + index1D; // easier reference to current point
				if(rand() % 2 == 0) {
					triArr[triInd] = Triangle(pt, pt + numY, pt + numY + 1);
					faces.push_back({index1D, index1D + numY, index1D + numY + 1});
					triArr[triInd+1] = Triangle(pt, pt + numY + 1, pt + 1);
					faces.push_back({index1D, index1D + numY + 1, index1D + 1});
				} else {
					triArr[triInd] = Triangle(pt, pt + 1, pt + numY);
					faces.push_back({index1D, index1D + 1, index1D + numY});
					triArr[triInd+1] = Triangle(pt + numY, pt + 1, pt + numY + 1);
					faces.push_back({index1D + numY, index1D + 1, index1D + numY + 1});
				}
				triInd += 2;
			}
		}
	}
	assert(triInd == numTri);

	// initialize edge dictionary from faces
	for(int i = 0; i < numTri; i++) {
		for(int j = 0; j < 3; j++) {
			array<int, 2> edge = {faces.at(i)[j], faces.at(i)[((j+1)%3)]};
			// ensure edges are ordered
			if(edge[0] > edge[1]) {
				edge[0] = edge[1];
				edge[1] = faces.at(i)[j];
			}
			if(edgeBelonging.find(edge) == edgeBelonging.end()) { // does not exist
				edgeBelonging[edge] = vector<int>();
			}
			edgeBelonging[edge].push_back(i);
		}		
	}

	// initialize integrator
	double maxLength = 2 * max(dx, dy); // generously round up maximum triangle side length

	// find space needed for results, one slot per gpu worker
	long long maxDivisions = (int) (maxLength/ds + 1); // max num samples per side, rounded up
	// maximum possible number of samples per triangle is loosely upper bounded by 2 * maxDivisions^2
	// assumming edge lengths are bounded above by maxDivisions * 2
	long long resultSlots = max(2 * maxDivisions * maxDivisions, (long long) maxX * maxY); // at least num pixels
	integrator.initialize(pixArr, maxX, maxY, APPROXTYPE, resultSlots);

	// create an initial approximation based on this triangulation
	updateApprox();
}

ConstantApprox::~ConstantApprox() {
	cudaFree(pixArr);
	cudaFree(points);
	cudaFree(triArr);
	cudaFree(grays);
	/*
	cudaFree(reds);
	cudaFree(greens);
	cudaFree(blues);
	*/
	delete[] imageInt;
}

double ConstantApprox::computeEnergy() {
	double totalEnergy = 0;
	for(int t = 0; t < numTri; t++) {
		totalEnergy += integrator.constantEnergyEval(triArr+t, grays[t], ds);
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

double ConstantApprox::step(double &prevEnergy, double &newEnergy, bool stringent) {
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
        	do {
            	undo();
        	} while (!gradUpdate());
        	updateApprox();
			newEnergy = computeEnergy();
		} while(stringent && newEnergy > prevEnergy);
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
				// for accuracy, use indices
				if(faces.at(t).at(v) != edge[0] && faces.at(t).at(v) != edge[1]) {
					// note cannot use triArr[t].vertices because the order of vertices
					// may be wrong (Triangle fixes order for ccw orientation)
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

void ConstantApprox::subdivide(int n) {
	vector<array<double, 3>> edgeEnergies;
	computeEdgeEnergies(&edgeEnergies);
	cout << "done computing energies" << endl;
	// sort by energies
	sort(edgeEnergies.begin(), edgeEnergies.end(), [](const array<double, 3> a, const array<double, 3> b) {
		return a[2] < b[2];
	});

	// proceed through edgeEnergies and extract up to n edges to subdivide
	set<int> trianglesToRemove; // hold indices of triangles in triArr
	vector<Point> newPoints; // points to append
	vector<array<int, 3>> newTriangles; // new faces to add

	int numDivided = 0; // number of edges already divided
	int curIndex = 0; // current edge being considered
	while(numDivided < n && curIndex < edgeEnergies.size()) {
		array<int, 2> edge = {(int) edgeEnergies.at(curIndex)[0], (int) edgeEnergies.at(curIndex)[1]};
		vector<int> incidentFaces = edgeBelonging.at(edge);
		// for clean subdivision, don't cut the same face twice
		bool alreadyDivided = false;
		for(int t : incidentFaces) {
			if(trianglesToRemove.find(t) != trianglesToRemove.end()) { 
				alreadyDivided = true;
			}
		}
		if(!alreadyDivided) { // this edge can be used
			// get new point to add to mesh
			double x0 = points[edge[0]].getX();
			double x1 = points[edge[1]].getX();
			double y0 = points[edge[0]].getY();
			double y1 = points[edge[1]].getY();
			bool borderX = points[edge[0]].isBorderX() && x0 == x1;
			bool borderY = points[edge[0]].isBorderY() && y0 == y1;
			Point midpoint((x0 + x1) / 2, (y0 + y1) / 2, borderX, borderY);
			newPoints.push_back(midpoint); // note overall index of this point is numPoints + numDivided
			// handle triangles
			for(int t : incidentFaces) {
				trianglesToRemove.insert(t);
				// get opposite vertex
				int oppositeInd;
				for(int v = 0; v < 3; v++) {
					if(faces.at(t).at(v) != edge[0] && faces.at(t).at(v) != edge[1]) {
						oppositeInd = faces.at(t).at(v);
					}
				}
				newTriangles.push_back({oppositeInd, edge[0], numPoints + numDivided});
				newTriangles.push_back({oppositeInd, edge[1], numPoints + numDivided});
			}
			numDivided++;
		}
		curIndex++;
	}
	cout << "subdivisions extracted" << endl;

	stepSize = originalStep;
	updateMesh(&newPoints, &newTriangles, &trianglesToRemove);
	updateApprox();
}

void ConstantApprox::updateMesh(vector<Point> *newPoints, vector<array<int, 3>> *newFaces, set<int> *discardedFaces) {
	vector<Point> oldPoints = getVertices();
	// free old memory
	cudaFree(points);
	cudaFree(triArr);
	cudaFree(grays);
	delete[] imageInt;
	// reallocate space
	int oldNumPoints = numPoints;
	numPoints += newPoints->size();
	numTri += newFaces->size() - discardedFaces->size();
	cudaMallocManaged(&points, numPoints * sizeof(Point));
	cudaMallocManaged(&triArr, numTri * sizeof(Triangle));
	cudaMallocManaged(&grays, numTri * sizeof(double));
	imageInt = new double[numTri];

	// load points
	for(int i = 0; i < oldNumPoints; i++) {
		points[i] = oldPoints.at(i);
	}
	for(int i = 0; i < newPoints->size(); i++) {
		points[oldNumPoints + i] = newPoints->at(i);
	}

	// handle triangles
	// first remove triangles that were split by going in reverse order, since sets are sorted (?)
	for(auto f = discardedFaces->rbegin(); f != discardedFaces->rend(); f++) {
		faces.erase(faces.begin() + *f);
	}
	// add new triangles
	for(auto f = newFaces->begin(); f != newFaces->end(); f++) {
		faces.push_back(*f);
	}

	// update triArr
	for(int i = 0; i < numTri; i++) {
		array<int, 3> t = faces.at(i);
		triArr[i] = Triangle(points + t[0], points + t[1], points + t[2]);
	}

	// update edges for next subdivision
	// TODO: don't clear the whole array, remove only the necessary parts
	edgeBelonging.clear();
	for(int i = 0; i < numTri; i++) {
		for(int j = 0; j < 3; j++) {
			array<int, 2> edge = {faces.at(i)[j], faces.at(i)[((j+1)%3)]};
			// ensure edges are ordered
			if(edge[0] > edge[1]) {
				edge[0] = edge[1];
				edge[1] = faces.at(i)[j];
			}
			if(edgeBelonging.find(edge) == edgeBelonging.end()) { // does not exist
				edgeBelonging[edge] = vector<int>();
			}
			edgeBelonging[edge].push_back(i);
		}		
	}
	gradX.clear();
	gradY.clear();
}

double ConstantApprox::getStep() {
	return stepSize;
}

vector<Point> ConstantApprox::getVertices() {
	vector<Point> vertices;
	for(int i = 0; i < numPoints; i++) {
		vertices.push_back(points[i]);
	}
	return vertices;
}

vector<array<int, 3>> ConstantApprox::getFaces() {
	return faces;
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
