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

// compute double integral of func for a single pixel and single triangle triArr[t]
// pixArr is a 1D representation of image, where pixel (x,y) is at x * maxY + y
// results holds the result for each pixel
__global__ void pixDoubleInt(nvstd::function<double(double, double)> &func, Pixel *pixArr, int &maxX, int &maxY, Triangle *triArr, int &t, double *results) {
	int x = blockIdx.x * blockDim.x + threadIdx.x;
	int y = blockIdx.y * blockDim.y + threadIdx.y;
	int ind = x * maxY + y; // index in pixArr
	if(ind < maxX * maxY) { // check bounds
		double integral = pixArr[ind].doubleIntegral(func, triArr[t]);
		results[ind] = integral;
	}
}

double doubleIntEval(nvstd::function<double(double, double)> func, Pixel *pixArr, int &maxX, int &maxY, Triangle *triArr, int &t, double *results) {
	dim3 numBlocks((maxX + numThreadsX -1) / numThreadsX, (maxY + numThreadsY - 1) / numThreadsY);
	// compute integral in parallel
	pixDoubleInt<<<numBlocks, threadsPerBlock>>>(func, pixArr, maxX, maxY, triArr, t, results);
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

// compute the line integral contribution of a single pixel on triangle triArr[t]
__global__ void pixLineInt(nvstd::function<double(double, double)> &func, Pixel *pixArr, int &maxX, int &maxY, Triangle *triArr, int &t, double *results) {
	int x = blockIdx.x * blockDim.x + threadIdx.x;
	int y = blockIdx.y * blockDim.y + threadIdx.y;
	int ind = x * maxY + y;
	if (ind < maxX * maxY) {
		double integral = pixArr[ind].lineIntegral(func, triArr[t]);
		results[ind] = integral;
	}

}

double lineIntEval(nvstd::function<double(double, double)> func, Pixel *pixArr, int &maxX, int &maxY, Triangle *triArr, int &t, double *results) {
	dim3 numBlocks((maxX + numThreadsX - 1) / numThreadsX, (maxY + numThreadsY - 1) / numThreadsY);
	// integrate in parallel
	pixLineInt<<<numBlocks, threadsPerBlock>>>(func, pixArr, maxX, maxY, triArr, t, results);
	double answer = sumArray(results, maxX * maxY, results);
	cudaDeviceSynchronize();
	return answer;
}

/*
double run() {
	Pixel *pixArr;
	Triangle *triArr;
	// dimensions of fake image
	int maxX = 4;
	int maxY = 3;
	cudaMallocManaged(&pixArr, maxX * maxY *sizeof(Pixel));
	cudaMallocManaged(&triArr, sizeof(Triangle));
	for(int i = 0; i < maxX; i++) {
		for(int j = 0; j < maxY; j++) {
			pixArr[i*maxY + j] = Pixel(i, j, 255);
		}
	}
	Point a(0,0);
	Point b(2,0);
	Point c(0,1);
	Triangle t(&a, &b, &c);
	triArr[0] = t;
	
	auto id = [](double x, double y) {
		return 1;
	};
	double area = doubleIntEval(id, pixArr, maxX, maxY, triArr, 1);
	cudaFree(pixArr);
	cudaFree(triArr);
	return area;
}
*/






/*
// compute the contribution of a single pixel to the double integral
__global__ void doubleIntPixEval(function<double(double, double)> func, vector<vector<Pixel>> *pixVec, Triangle *triangle, double *blockResults, int maxX, int maxY) {
	// modified from https://stackoverflow.com/questions/19946286
	extern __shared__ double threads[][]; // hold results per thread
    int x = blockIdx.x * blockDim.x + threadIdx.x;
	int y = blockIdx.y * blockDim.y + threadIdx.y;
	// calculate pixel contribution
	double integral = 0;
	if(x < maxX && y < maxY) {
		integral = pixVec->at(x).at(y).doubleIntegral(func, *triangle);
	}
	threads[x][y] = integral;
	__syncthreads();

	// sum in x direction; divide summed range in half each time,
	// sum number in first half with corresponding number in second half
	for(int offset = blockDim.x / 2; offset > 0; offset /= 2) {
		if(threadId.x < offset) {
			threads[threadIdx.x][threadIdx.y] += threads[threadIdx.x + offset][threadIdx.y];
		}
		__syncthreads();
	}
	// threads[0][y] hold the x direction sums
	
	// now sum in y direction
	for(int offset = blockDim.y / 2; offset > 0; offset /= 2) {
		if(threadIdx.x == 0 && threadIdx.y < offset) {
			threads[0][threadIdx.y] += threads[0][threadIdx.y+offset];
		}
		__syncthreads();
	}

	if(threadIdx.x == 0 && threadIdx.y == 0) {
		blockResults[blockIdx.x * gridDim.y + blockIdx.y] = threads[0][0];
	}
}

double doubleIntEval(function<double(double, double)> func, vector<vector<Pixel>> *pixVec, Triangle *triangle) {
	int maxX = pixVec->size();
	int maxY = pixVec->at(0).size();
	int blocksizeX = 16;
	int blocksizeY = 16;
	dim3 threadsPerBlock(blocksizeX, blocksizeY);
	dim3 numBlocks((maxX + blocksizeX - 1)/blocksizeX, (maxY + blocksizeY - 1)/blocksizeY);
	// create array to hold block results
	double *blockResults = new double[numBlocks.x * numBlocks.y];
	doubleIntPixEval<<<numBlocks,threadsPerBlock>>>(func, pixVec, triangle, blockResults, maxX, maxY);
}

*/
