#include "parallelInt.cuh"

// thread setup
int numThreadsX = 32;
int numThreadsY = 32;
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
__global__ void pixConstantDoubleInt(Pixel *pixArr, int &maxX, int &maxY, Triangle *triArr, int &t, double *results) {
	int x = blockIdx.x * blockDim.x + threadIdx.x;
	int y = blockIdx.y * blockDim.y + threadIdx.y;
	int ind = x * maxY + y; // index in pixArr
	if(ind < maxX * maxY) { // check bounds
		double area = pixArr[ind].intersectionArea(triArr[t]);
		results[ind] = area * pixArr[ind].getColor();
	}
}

double doubleIntEval(ApproxType approx, Pixel *pixArr, int &maxX, int &maxY, Triangle *triArr, int &t, double *results) {
	dim3 numBlocks((maxX + numThreadsX -1) / numThreadsX, (maxY + numThreadsY -1) / numThreadsY);
	// compute integral in parallel based on function to integrate
	switch (approx) {
		case constant:
			pixConstantDoubleInt<<<numBlocks, threadsPerBlock>>>(pixArr, maxX, maxY, triArr, t, results);
			break;
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
__global__ void pixConstantEnergyInt(Pixel *pixArr, int &maxX, int &maxY, Triangle *triArr, double *colors, int &t, double *results) {
	int x = blockIdx.x * blockDim.x + threadIdx.x;
	int y = blockIdx.y * blockDim.y + threadIdx.y;
	int ind = x * maxY + y; // index in pixArr;
	if(ind < maxX * maxY) {
		double area = pixArr[ind].intersectionArea(triArr[t]);
		double diff = colors[t] - pixArr[ind].getColor();
		results[ind] = diff * diff * area;
	}
}

double constantEnergyEval(Pixel *pixArr, int &maxX, int &maxY, Triangle *triArr, double *colors, int &numTri, double *results) {
	dim3 numBlocks((maxX + numThreadsX - 1) / numThreadsX, (maxY + numThreadsY - 1) / numThreadsY);
	double totalEnergy = 0;
	for(int t = 0; t < numTri; t++) {
		pixConstantEnergyInt<<<numBlocks, threadsPerBlock>>>(pixArr, maxX, maxY, triArr, colors, t, results);
		totalEnergy += sumArray(results, maxX * maxY, results); // add energy for this triangle
	}
	cudaDeviceSynchronize(); // wait to finish
	return totalEnergy;
}

// compute line integral of v dot n f ds for a single pixel and single triangle a, b, c when point b is moving
__global__ void pixConstantLineInt(Pixel *pixArr, int &maxX, int &maxY, Point *a, Point *b, Point *c, bool isX, double *results) {
	int x = blockIdx.x * blockDim.x + threadIdx.x;
	int y = blockIdx.y * blockDim.y + threadIdx.y;
	int ind = x * maxY + y;
	if (ind < maxX * maxY) {
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
				double velX = (isX) ? scale : 0;
				Matrix v(velX, scale - velX); // velocity vector
				Matrix n = seg.unitNormal();
				double vn = v.transpose().multiply(n).get(0,0); // average value of v * n
				answer += vn * length * pixArr[ind].getColor();
			}
		}
		results[ind] = answer;
	}
}

double lineIntEval(ApproxType approx, Pixel *pixArr, int &maxX, int &maxY, Triangle *triArr, int &t, int &pt, bool isX, double *results) {
	dim3 numBlocks((maxX + numThreadsX - 1) / numThreadsX, (maxY + numThreadsY - 1) / numThreadsY);
	Point vertices[3]; // vertices of triArr[t]
	triArr[t].copyVertices(vertices, vertices+1, vertices+2);
	// compute integral in parallel based on function to integrate
	switch (approx) {
		case constant:
			pixConstantLineInt<<<numBlocks, threadsPerBlock>>>(pixArr, maxX, maxY, vertices + ((pt+2)%3), vertices + pt, vertices + ((pt+1)%3), isX, results);
			break;
		case linear: // TODO
			break;
		case quadratic: // TODO
			break;
	}
	double answer = sumArray(results, maxX * maxY, results);
	cudaDeviceSynchronize();
	return answer;
}