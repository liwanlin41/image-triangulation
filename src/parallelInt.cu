#include "parallelInt.cuh"

// thread setup
int numThreadsX = 32;
int numThreadsY = 16;
dim3 threadsPerBlock(numThreadsX, numThreadsY);
// NOTE: block grid will have to be set within each function

// compute the sum of an array arr with given size, in parallel
// with 1D thread/blocks, storing the result per block in result
__global__ void sumBlock(double *arr, int size, double *result) {
	extern __shared__ double partial[]; // hold partial results
	int tid = threadIdx.x;
	int ind = blockIdx.x * blockDim.x + tid;
	// load into partial result array
	if(ind < size) {
		partial[tid] = arr[ind];
	} else {
		partial[tid] = 0;
	}
	__syncthreads();

	for(int step = blockDim.x / 2; step > 0; step /= 2) {
		if(tid < step) {
			partial[tid] += partial[tid + step];
		}
		__syncthreads();
	}

	// write output for block to result
	if(tid == 0) {
		result[blockIdx.x] = partial[0];
	}
}

// quickly sum an array with given size in parallel and return the result;
// NOTE: arr, partialRes must already be shared between host and device 
double sumArray(double *arr, int size, double *partialRes) {
	int numThreads = 1024; // threads per block
	// shared memory size for device
	int memSize = numThreads * sizeof(double);
	int numBlocks = (size + numThreads - 1) / numThreads;
	sumBlock<<<numBlocks, numThreads, memSize>>>(arr, size, partialRes);
	// number of elements to sum is now numBlocks
	// number of blocks for next iteration
	int newNumBlocks = (numBlocks + numThreads - 1) / numThreads;
	// repeat until all elements have been summed
	while(numBlocks > 1) {
		sumBlock<<<newNumBlocks,numThreads, memSize>>>(partialRes, numBlocks, partialRes); 
		numBlocks = newNumBlocks;
		newNumBlocks = (newNumBlocks + numThreads - 1) / numThreads;
	}
	// at this point the array has been summed and the result is in partialRes[0]
	cudaDeviceSynchronize();
	double output = partialRes[0];
	return output;
}

// compute double integral of f dA for a single pixel and single triangle triArr[t]
// pixArr is a 1D representation of image, where pixel (x, y) is at x * maxY + y
// reults holds the result for each pixel
__global__ void pixConstantDoubleInt(Pixel *pixArr, int maxX, int maxY, Triangle *triArr, int t, double *results, ColorChannel channel) {
	int x = blockIdx.x * blockDim.x + threadIdx.x;
	int y = blockIdx.y * blockDim.y + threadIdx.y;
	int ind = x * maxY + y; // index in pixArr
	if(x < maxX && y < maxY) { // check bounds
		//double area = pixArr[ind].intersectionArea(triArr[t]);
		double area = pixArr[ind].approxArea(triArr[t]);
		results[ind] = area * pixArr[ind].getColor(channel);
	}
}

double doubleIntEval(ApproxType approx, Pixel *pixArr, int &maxX, int &maxY, Triangle *triArr, int &t, double *results, ColorChannel channel) {
	dim3 numBlocks((maxX + numThreadsX -1) / numThreadsX, (maxY + numThreadsY -1) / numThreadsY);
	// compute integral in parallel based on function to integrate
	switch (approx) {
		case constant: {
			pixConstantDoubleInt<<<numBlocks, threadsPerBlock>>>(pixArr, maxX, maxY, triArr, t, results, channel);
			break;
		}
		case linear: // TODO: fill out
			break;
		case quadratic: // TODO: fill out
			break;
	}
	double answer = sumArray(results, maxX * maxY, results);
	cudaDeviceSynchronize(); // wait for everything to finish
	return answer;
}

// compute the energy of a single pixel on triangle triArr[t]
__global__ void pixConstantEnergyInt(Pixel *pixArr, int maxX, int maxY, Triangle *triArr, double *colors, int t, double *results) {
	int x = blockIdx.x * blockDim.x + threadIdx.x;
	int y = blockIdx.y * blockDim.y + threadIdx.y;
	int ind = x * maxY + y; // index in pixArr;
	if(x < maxX && y < maxY) {
		//double area = pixArr[ind].intersectionArea(triArr[t]);
		double area = pixArr[ind].approxArea(triArr[t]);
		double diff = colors[t] - pixArr[ind].getColor();
		results[ind] = diff * diff * area;
	}
}

// using Point a as vertex point, sample ~samples^2/2 points inside the triangle with an area element of dA
// NOTE: samples does not count endpoints along edge bc as the parallelograms rooted there lie outside the triangle
// maxY is for converting 2D pixel index to 1D index
__global__ void approxConstantEnergySample(Pixel *pixArr, int maxY, Point *a, Point *b, Point *c, double color, double *results, double dA, int samples) {
	int u = blockIdx.x * blockDim.x + threadIdx.x; // component towards b
	int v = blockIdx.y * blockDim.y + threadIdx.y; // component towards c
	int ind = (2 * samples - u + 1) * u / 2 + v; // 1D index in results
	// this is because there are s points in the first column, s-1 in the next, etc. up to s - u + 1
	if(u + v < samples) {
		// get coordinates of this point using appropriate weights
		double x = (a->getX() * (samples - u - v) + b->getX() * u + c->getX() * v) / samples;
		double y = (a->getY() * (samples - u - v) + b->getY() * u + c->getY() * v) / samples;
		// find containing pixel
		int pixX = pixelRound(x);
		int pixY = pixelRound(y);
		double diff = color - pixArr[pixX * maxY + pixY].getColor();
		if(u + v == samples - 1) { // along the opposite edge, dA element is a triangle instead of a parallelogram
			results[ind] = diff * diff * dA / 2;
		} else {
			results[ind] = diff * diff * dA;
		}
	}
}

double constantEnergyApprox(Pixel *pixArr, int &maxY, Triangle *triArr, double *colors, int &numTri, double *results, double ds, Point *workingTri) {
	double totalEnergy = 0;
	for(int t = 0; t < numTri; t++) {
		// compute number of samples needed
		double maxLength = triArr[t].maxLength();
		int samples = ceil(maxLength/ds);
		int i = triArr[t].minVertex(); // vertex opposite shortest side
		// ensure minVertex is copied into location workingTri
		triArr[t].copyVertices(workingTri+((3-i)%3), workingTri+((4-i)%3), workingTri+((5-i)%3));
		// unfortunately half of these threads will not be doing useful work; fix this somehow?
		dim3 numBlocks((samples + numThreadsX - 1) / numThreadsX, (samples + numThreadsY - 1) / numThreadsY);
		double dA = triArr[t].getArea() * 2 / (samples * samples);
		approxConstantEnergySample<<<numBlocks, threadsPerBlock>>>(pixArr, maxY, workingTri, workingTri + 1, workingTri + 2, colors[t], results, dA, samples);
		totalEnergy += sumArray(results, samples * (samples + 1) / 2, results);
	}
	return totalEnergy;
}

double constantEnergyEval(Pixel *pixArr, int &maxX, int &maxY, Triangle *triArr, double *colors, int &numTri, double *results) {
	dim3 numBlocks((maxX + numThreadsX - 1) / numThreadsX, (maxY + numThreadsY - 1) / numThreadsY);
	double totalEnergy = 0;
	for(int t = 0; t < numTri; t++) {
		pixConstantEnergyInt<<<numBlocks, threadsPerBlock>>>(pixArr, maxX, maxY, triArr, colors, t, results);
		totalEnergy += sumArray(results, maxX * maxY, results); // add energy for this triangle
	}
	return totalEnergy;
}

// compute line integral of v dot n f ds for a single pixel and single triangle a, b, c when point b is moving
__global__ void pixConstantLineInt(Pixel *pixArr, int maxX, int maxY, Point *a, Point *b, Point *c, bool isX, double *results) {
	int x = blockIdx.x * blockDim.x + threadIdx.x;
	int y = blockIdx.y * blockDim.y + threadIdx.y;
	int ind = x * maxY + y;
	if (x < maxX && y < maxY) {
		double answer = 0;
		for(int i = 0; i < 2; i++) { // v dot n is nonzero only on a -- b and b -- c
			// extract segment and maintain ccw order for outward normal
			Segment seg = (i == 0) ? Segment(a, b) : Segment(b, c);
			Point *segEnd = (i == 0) ? a : c; // determine endpoint of seg that is not b
			double midX, midY; // to hold midpoint of segment intersection with this pixel
			double length = pixArr[ind].intersectionLength(seg, &midX, &midY);
			if(length != 0) {
				Point midpoint(midX, midY);
				// compute velocity at this point by scaling
				double distanceToVertex = midpoint.distance(*segEnd);
				double scale = distanceToVertex / seg.length(); // 1 if at b, 0 at opposite edge
				// velocity components
				double velX = (isX) ? scale : 0;
				double velY = scale - velX;
				// get unit normal values for this segment
				double nx, ny;
				seg.unitNormal(&nx, &ny);
				double vn = velX * nx + velY * ny; // average value of v * n
				answer += vn * length * pixArr[ind].getColor();
			}
		}
		results[ind] = answer;
	}
}

double lineIntEval(ApproxType approx, Pixel *pixArr, int &maxX, int &maxY, Triangle *triArr, int &t, int &pt, bool isX, double *results, Point *workingTri) {
	dim3 numBlocks((maxX + numThreadsX - 1) / numThreadsX, (maxY + numThreadsY - 1) / numThreadsY);
	triArr[t].copyVertices(workingTri, workingTri+1, workingTri+2);
	// compute integral in parallel based on function to integrate
	switch (approx) {
		case constant: {
			pixConstantLineInt<<<numBlocks, threadsPerBlock>>>(pixArr, maxX, maxY, workingTri+((pt+2)%3), workingTri+pt, workingTri+((pt+1)%3), isX, results);
			break;
		}
		case linear: // TODO
			break;
		case quadratic: // TODO
			break;
	}
	double answer = sumArray(results, maxX * maxY, results);
	return answer;
}