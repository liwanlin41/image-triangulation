#ifndef matrix_cuh
#define matrix_cuh

#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#define CUDA_DEV __device__
#else
#define CUDA_HOSTDEV
#define CUDA_DEV
#endif

/**
 * represent an immutable rectangular matrix of arbitrary size,
 * with special support for 2x2 matrices
 */

class Matrix {
	private:
		int numRows, numCols;
		double *matrix;
	public:
		// constructor assuming input vector is rectangular
		//__device__ Matrix(device_vector<double> const &inputVec, int rows, int cols);
		CUDA_HOSTDEV Matrix(double *inputArr, int rows, int cols);
		// special constructor for 2x2 matrix, creating the matrix [[a,b],[c,d]]
		CUDA_HOSTDEV Matrix(double a, double b, double c, double d);
		// special constructor for 2x1 vector matrix
		CUDA_HOSTDEV Matrix(double a, double b);

		// copy constructor
		CUDA_HOSTDEV Matrix(const Matrix &m);
		// copy assignment
		CUDA_HOSTDEV Matrix& operator=(const Matrix &m);
		
		CUDA_HOSTDEV ~Matrix();

		CUDA_DEV int getNumRows();
		CUDA_DEV int getNumCols();
		// return the element in row i, column j of the matrix
		CUDA_HOSTDEV double get(int i, int j);
		// return result of this times other for compatible matrices
		CUDA_HOSTDEV Matrix multiply(Matrix &other);

		// compute determinant of 2x2 matrix
		CUDA_HOSTDEV double determinant();
		// compute adjugate of 2x2 matrix
		CUDA_DEV Matrix adjugate();
		// compute inverse of 2x2 matrix
		CUDA_DEV Matrix inverse();
	
		// transpose a matrix of arbitrary dimension
		CUDA_HOSTDEV Matrix transpose();
};

#endif
