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
	int curSize = size; // current length of array to sum
	int numBlocks = (size + numThreads - 1) / numThreads;
	bool ansArr = true; // whether results are currently held in arr
	while(curSize > 1) {
		if(ansArr) {
			sumBlock<<<numBlocks, numThreads, memSize>>>(arr, curSize, partialRes);
		} else {
			sumBlock<<<numBlocks, numThreads, memSize>>>(partialRes, curSize, arr);
		}
		curSize = numBlocks;
		numBlocks = (numBlocks + numThreads - 1) / numThreads;
		ansArr = !ansArr;
	}
	// at this point the array has been summed
	cudaDeviceSynchronize();
	if(ansArr) { // arr should hold the results
		return arr[0];
	}
	return partialRes[0];
}

// compute double integral of f dA for a single pixel and single triangle triArr[t]
// pixArr is a 1D representation of image, where pixel (x, y) is at x * maxY + y
// reults holds the result for each pixel
__global__ void pixConstantDoubleInt(Pixel *pixArr, int maxX, int maxY, Triangle *triArr, int t, double *results, ColorChannel channel) {
	int x = blockIdx.x * blockDim.x + threadIdx.x;
	int y = blockIdx.y * blockDim.y + threadIdx.y;
	int ind = x * maxY + y; // index in pixArr
	if(x < maxX && y < maxY) { // check bounds
		double area = pixArr[ind].intersectionArea(triArr[t]);
		//double area = pixArr[ind].approxArea(triArr[t]);
		results[ind] = area * pixArr[ind].getColor(channel);
	}
}

double doubleIntEval(ApproxType approx, Pixel *pixArr, int &maxX, int &maxY, Triangle *triArr, int &t, double *res0, double *res1, ColorChannel channel) {
	dim3 numBlocks((maxX + numThreadsX -1) / numThreadsX, (maxY + numThreadsY -1) / numThreadsY);
	// compute integral in parallel based on function to integrate
	switch (approx) {
		case constant: {
			pixConstantDoubleInt<<<numBlocks, threadsPerBlock>>>(pixArr, maxX, maxY, triArr, t, res0, channel);
			break;
		}
		case linear: // TODO: fill out
			break;
		case quadratic: // TODO: fill out
			break;
	}
	double answer = sumArray(res0, maxX * maxY, res1);
	cudaDeviceSynchronize(); // wait for everything to finish
	return answer;
}

// using Point a as vertex point, sample ~samples^2/2 points inside triangle with area element of dA
// for details see approxConstantEnergySample below
__global__ void constDoubleIntSample(Pixel *pixArr, int maxY, Point *a, Point *b, Point *c, double *results, double dA, int samples, ColorChannel channel) {
	int u = blockIdx.x * blockDim.x + threadIdx.x; // component towards b
	int v = blockIdx.y * blockDim.y + threadIdx.y; // component towards c
	int ind = (2 * samples - u + 1) * u / 2 + v; // 1D index in results
	if(u + v < samples) {
		double x = (a->getX() * (samples - u - v) + b->getX() * u + c->getX() * v) / samples;
		double y = (a->getY() * (samples - u - v) + b->getY() * u + c->getY() * v) / samples;
		// find containing pixel
		int pixX = pixelRound(x);
		int pixY = pixelRound(y);
		double areaContrib = (u+v == samples - 1) ? dA : 2 * dA;
		results[ind] = pixArr[pixX * maxY + pixY].getColor(channel) * areaContrib;
	}
}

double doubleIntApprox(ApproxType approx, Pixel *pixArr, int &maxY, Triangle *tri, double *res0, double *res1, double &ds, Point *workingTri, ColorChannel channel) {
	int i = tri->midVertex();
	tri->copyVertices(workingTri+((3-i)%3), workingTri+((4-i)%3), workingTri+((5-i)%3));
	// compute number of samples
	int samples = ceil(workingTri[1].distance(workingTri[2])/ds);
	dim3 numBlocks((samples + numThreadsX - 1) / numThreadsX, (samples + numThreadsY - 1) / numThreadsY);
	double dA = tri->getArea() / (samples * samples);
	switch(approx) {
		case constant: {
			constDoubleIntSample<<<numBlocks, threadsPerBlock>>>(pixArr, maxY, workingTri, workingTri+1, workingTri+2, res0, dA, samples, channel);
			cudaError_t error = cudaGetLastError();
  			if(error != cudaSuccess)
  			{
    			// print the CUDA error message and exit
				printf("CUDA error in double int: %s\n", cudaGetErrorString(error));
				exit(-1);
			}
			break;
		}
		case linear:
			break;
		case quadratic:
			break;
	}
	double answer = sumArray(res0, samples * (samples + 1) / 2, res1);
	return answer;
}

// compute line integral of v dot n f ds where point a is moving; 
// reverse determines if integral should be computed from a to b (false) or opposite
__global__ void constLineIntSample(Pixel *pixArr, int maxY, Point *a, Point *b, bool reverse, bool isX, double *results, double ds, int samples) {
	int k = blockIdx.x * blockDim.x + threadIdx.x; // index along a to b
	if(k < samples) {
		// extract current point and containing pixel
		double x = (a->getX() * (samples - k) + b->getX() * k) / samples;
		double y = (a->getY() * (samples - k) + b->getY() * k) / samples;
		int pixX = pixelRound(x);
		int pixY = pixelRound(y);
		// velocity components
		double scale = ((double) samples - k) / samples; // 1 when k = 0 (evaluate at a) and 0 at b
		double velX = (isX) ? scale : 0;
		double velY = scale - velX;
		// extract unit normal, manually for the sake of speed
		double length = a->distance(*b); // length of whole segment
		// assume going from a to b first, want normal pointing right
		double nx = (b->getY() - a->getY()) / length;
		double ny = (a->getX() - b->getX()) / length;
		double vn = velX * nx + velY * ny; // value of v * n at this point
		// flip vn if normal is actually pointing the other way (integrate from b to a)
		if(reverse) vn *= -1;
		results[k] = vn * ds * pixArr[pixX * maxY + pixY].getColor();
	}
}

double lineIntApprox(ApproxType approx, Pixel *pixArr, int &maxY, Triangle *triArr, int &t, int &pt, bool isX, double *results0, double *results1, double ds, Point *workingTri) {
	// ensure pt is copied into workingTri
	triArr[t].copyVertices(workingTri+((3-pt)%3), workingTri+((4-pt)%3), workingTri+((5-pt)%3));
	int numThreads = numThreadsX * numThreadsY; // this will be a 1D kernel call
	// get number of samples for side pt, pt+1 and side pt, pt+2
	int samples[2];
	int numBlocks[2];
	for(int i = 0; i < 2; i++) {
		// increase the number of samples because space and time both allow
		samples[i] = ceil(10*workingTri->distance(workingTri[i+1])/ds);
		numBlocks[i] = ceil(1.0 * samples[i] / numThreads);
	};
	double answer = 0; // integrate over both moving sides
	switch(approx) {
		case constant: {
			for(int i = 0; i < 2; i++) {
				double totalLength = workingTri->distance(workingTri[i+1]);
				// actual dx being used
				double dx = totalLength / samples[i];
				constLineIntSample<<<numBlocks[i], numThreads>>>(pixArr, maxY, workingTri, workingTri+i+1, (i==1), isX, results0, dx, samples[i]);
				answer += sumArray(results0, samples[i], results1);
			}
		}
		case linear:
			break;
		case quadratic:
			break;
	}
	return answer;
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

double lineIntEval(ApproxType approx, Pixel *pixArr, int &maxX, int &maxY, Triangle *triArr, int &t, int &pt, bool isX, double *res0, double *res1, Point *workingTri) {
	dim3 numBlocks((maxX + numThreadsX - 1) / numThreadsX, (maxY + numThreadsY - 1) / numThreadsY);
	triArr[t].copyVertices(workingTri, workingTri+1, workingTri+2);
	// compute integral in parallel based on function to integrate
	switch (approx) {
		case constant: {
			pixConstantLineInt<<<numBlocks, threadsPerBlock>>>(pixArr, maxX, maxY, workingTri+((pt+2)%3), workingTri+pt, workingTri+((pt+1)%3), isX, res0);
			break;
		}
		case linear: // TODO
			break;
		case quadratic: // TODO
			break;
	}
	double answer = sumArray(res0, maxX * maxY, res1);
	return answer;
}

ParallelIntegrator::ParallelIntegrator() {
	threads2D = dim3(threadsX, threadsY);
}

bool ParallelIntegrator::initialize(Pixel *pix, Triangle *tri, int xMax, int yMax, ApproxType a, long long space, bool exact) {
	// allocate working computation space
	cudaMallocManaged(&arr, space * sizeof(double));
	cudaMallocManaged(&helper, space * sizeof(double));
	// the above may cause errors because so much memory is required
	cudaError_t error = cudaGetLastError();
  	if(error != cudaSuccess) {
		printf("CUDA error: %s\n", cudaGetErrorString(error));
		cout << "An approximation of this quality is not possible due to memory limitations." << endl;
		return false;
	}
	cudaMallocManaged(&curTri, 3 * sizeof(Point));
	if(cudaGetLastError() != cudaSuccess) {
		printf("CUDA error: %s\n", cudaGetErrorString(error));
		return false;
	}
	// steal references for easy access later
	pixArr = pix;
	triArr = tri;
	approx = a;
	maxX = xMax;
	maxY = yMax;
	computeExact = exact;
	return true;
}

ParallelIntegrator::~ParallelIntegrator() {
	cudaFree(arr);
	cudaFree(helper);
	cudaFree(curTri);
}

double ParallelIntegrator::sumArray(int size) {
	// shared memory size for device
	int memSize = threads1D * sizeof(double);
	int curSize = size; // current length of array to sum
	int numBlocks = (size + threads1D - 1) / threads1D;
	bool ansArr = true; // whether results are currently held in arr
	while(curSize > 1) {
		if(ansArr) {
			sumBlock<<<numBlocks, threads1D, memSize>>>(arr, curSize, helper);
		} else {
			sumBlock<<<numBlocks, threads1D, memSize>>>(helper, curSize, arr);
		}
		curSize = numBlocks;
		numBlocks = (numBlocks + threads1D - 1) / threads1D;
		ansArr = !ansArr;
	}
	// at this point the array has been summed
	cudaDeviceSynchronize();
	if(ansArr) { // arr should hold the results
		return arr[0];
	}
	return helper[0];
}

// kernel for constantEnergyEval
// compute the energy of a single pixel on triangle triArr[t]
__global__ void pixConstantEnergyInt(Pixel *pixArr, int maxX, int maxY, Triangle *triArr, double *colors, int t, double *results) {
	int x = blockIdx.x * blockDim.x + threadIdx.x;
	int y = blockIdx.y * blockDim.y + threadIdx.y;
	int ind = x * maxY + y; // index in pixArr;
	if(x < maxX && y < maxY) {
		double area = pixArr[ind].intersectionArea(triArr[t]);
		//double area = pixArr[ind].approxArea(triArr[t]);
		double diff = colors[t] - pixArr[ind].getColor();
		results[ind] = diff * diff * area;
	}
}

double ParallelIntegrator::constantEnergyExact(double *colors, int numTri) {
	dim3 numBlocks((maxX + threadsX - 1) / threadsX, (maxY + threadsY - 1) / threadsY);
	double totalEnergy = 0;
	for(int t = 0; t < numTri; t++) {
		pixConstantEnergyInt<<<numBlocks, threads2D>>>(pixArr, maxX, maxY, triArr, colors, t, arr);
		totalEnergy += sumArray(maxX * maxY);
	}
	return totalEnergy;
}

// kernel for constant energy approx
// using Point a as vertex point, sample ~samples^2/2 points inside the triangle with a triangular area element of dA
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
		// account for points near edge bc having triangle contributions rather than parallelograms,
		// written for fast access and minimal branching
		double areaContrib = (u + v == samples - 1) ? dA : 2 * dA;
		results[ind] = diff * diff * areaContrib;
	}
}

double ParallelIntegrator::constantEnergyApprox(double *colors, int numTri, double ds) {
	double totalEnergy = 0;
	for(int t = 0; t < numTri; t++) {
		int i = triArr[t].midVertex(); // vertex opposite middle side
		// ensure minVertex is copied into location workingTri
		triArr[t].copyVertices(curTri+((3-i)%3), curTri+((4-i)%3), curTri+((5-i)%3));
		// compute number of samples needed, using median number per side as reference
		int samples = ceil(curTri[1].distance(curTri[2])/ds);
		// unfortunately half of these threads will not be doing useful work; no good way to fix this, sqrt is too slow
		dim3 numBlocks((samples + threadsX - 1) / threadsX, (samples + threadsY - 1) / threadsY);
		double dA = triArr[t].getArea() / (samples * samples);
		approxConstantEnergySample<<<numBlocks, threads2D>>>(pixArr, maxY, curTri, curTri + 1, curTri + 2, colors[t], arr, dA, samples);
		/*
		cudaError_t error = cudaGetLastError();
  		if(error != cudaSuccess)
  		{
    		// print the CUDA error message and exit
			printf("CUDA error in energy: %s\n", cudaGetErrorString(error));
			printf("num samples is %d\n", samples);
			exit(-1);
		  }
		*/
		totalEnergy += sumArray(samples * (samples + 1) / 2);
	}
	return totalEnergy;
}

double ParallelIntegrator::constantEnergyEval(double *colors, int numTri, double ds) {
	// switch integration method based on exactness required
	if(computeExact) {
		return constantEnergyExact(colors, numTri);
	}
	return constantEnergyApprox(colors, numTri, ds);
}