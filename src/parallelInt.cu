#include "parallelInt.h"

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
// arr must already be shared between host and device 
double sumArray(double *arr, int size) {
	int numThreads = 1024; // threads per block
	// shared memory size for device
	int memSize = numThreads * sizeof(double);
	int numBlocks = (size + numThreads - 1) / numThreads;
	double *partialRes; // block sums
	cudaMallocManaged(&partialRes, numBlocks*sizeof(double));
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
	// remember to free memory
	cudaFree(partialRes);
	return output;
}

double runSum(double *arr, int size) {
	double *sharedArr;
	cudaMallocManaged(&sharedArr, size*sizeof(double));
	for(int i = 0; i < size; i++) {
		sharedArr[i] = arr[i];
	}
	double sum = sumArray(sharedArr, size);
	cudaFree(sharedArr);
	return sum;
}

double doubleIntEval(function<double(double, double)> func, vector<vector<Pixel>> *pixVec, Triangle *triangle) {
	return 0;
}

double lineIntEval(function<double(double, double)> func, vector<vector<Pixel>> *pixVec, Triangle *triangle) {
	return 0;
}






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