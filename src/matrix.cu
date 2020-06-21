#include "matrix.cuh"

Matrix::Matrix(double *inputArr, int rows, int cols) :
	numRows(rows), numCols(cols) {
	// defensive copy
	matrix = new double[rows * cols];
	for(int i = 0; i < rows * cols; i++) {
		matrix[i] = inputArr[i];
	}
}

Matrix::Matrix(double a, double b, double c, double d) {
	numRows = 2;
	numCols = 2;
	matrix = new double[4];
	matrix[0] = a;
	matrix[1] = b;
	matrix[2] = c;
	matrix[3] = d;
}

Matrix::Matrix(double a, double b) {
	numRows = 2;
	numCols = 1;
	matrix = new double[2];
	matrix[0] = a;
	matrix[1] = b;
}

__device__ int Matrix::getNumRows() {
	return numRows;
}

__device__ int Matrix::getNumCols() {
	return numCols;
}

double Matrix::get(int i, int j) {
	return matrix[i * numCols + j];
}

Matrix Matrix::multiply(Matrix &other) {
	int rows = numRows;
	int cols = other.numCols;
	double *product = new double[rows * cols];
	for(int i = 0; i < rows; i++) {
		for(int j = 0; j < cols; j++) {
			double dotProd = 0;
			for(int k = 0; k < numCols; k++) {
				dotProd += get(i, k) * other.get(k, j);
			}
			product[i * cols + j] = dotProd;
		}
	}
	Matrix result(product, rows, cols);
	delete[] product;
	return result;
}

double Matrix::determinant() {
	return matrix[0] * matrix[3] - matrix[1] * matrix[2];
}

__device__ Matrix Matrix::adjugate() {
	return Matrix(matrix[3], -matrix[1], -matrix[2], matrix[0]);
}
	
__device__ Matrix Matrix::inverse() {
	double det = determinant();
	return Matrix(matrix[3]/det, -matrix[1]/det, -matrix[2]/det, matrix[0]/det);
}

Matrix Matrix::transpose() {
	double *transposeArr = new double[numRows * numCols];
	for(int i = 0; i < numRows; i++) {
		for(int j = 0; j < numCols; j++) {
			transposeArr[j * numRows + i] = get(i, j);
		}
	}
	Matrix result(transposeArr, numCols, numRows);
	delete[] transposeArr;
	return result;
}

Matrix::~Matrix() {
	delete[] matrix;
}

Matrix::Matrix(const Matrix &m) {
	numRows = m.numRows;
	numCols = m.numCols;
	matrix = new double[numRows * numCols];
	for(int i = 0; i < numRows * numCols; i++) {
		matrix[i] = m.matrix[i];
	}
}

Matrix& Matrix::operator=(const Matrix &m) {
	if (this != &m) {
		numRows = m.numRows;
		numCols = m.numCols;
		double* newArray = new double[numRows * numCols];
		for(int i = 0; i < numRows * numCols; i++) {
			newArray[i] = m.matrix[i];
		}
		delete[] matrix;
		matrix = newArray;
	}
	return *this;
}
