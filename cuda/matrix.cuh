#ifndef matrix_cuh
#define matrix_cuh

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
		__host__ __device__ Matrix(double *inputArr, int rows, int cols);
		// special constructor for 2x2 matrix, creating the matrix [[a,b],[c,d]]
		__host__ __device__ Matrix(double a, double b, double c, double d);
		// special constructor for 2x1 vector matrix
		__host__ __device__ Matrix(double a, double b);

		// copy constructor
		__host__ __device__ Matrix(const Matrix &m);
		// copy assignment
		__host__ __device__ Matrix& operator=(const Matrix &m);
		
		__host__ __device__ ~Matrix();

		__device__ int getNumRows();
		__device__ int getNumCols();
		// return the element in row i, column j of the matrix
		__host__ __device__ double get(int i, int j);
		// return result of this times other for compatible matrices
		__host__ __device__ Matrix multiply(Matrix &other);

		// compute determinant of 2x2 matrix
		__host__ __device__ double determinant();
		// compute adjugate of 2x2 matrix
		__device__ Matrix adjugate();
		// compute inverse of 2x2 matrix
		__device__ Matrix inverse();
	
		// transpose a matrix of arbitrary dimension
		__host__ __device__ Matrix transpose();
};

#endif
