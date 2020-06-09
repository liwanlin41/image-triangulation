#include "matrix.hpp"

Matrix::Matrix(vector<vector<double>> const &inputArr) {
    // defensive copy
    numRows = inputArr.size();
    numCols = inputArr.at(0).size();
    // fill matrix with dummy values
    matrix.resize(numRows, vector<double>(numCols, 0));
    for(int i = 0; i < numRows; i++) {
        for(int j = 0; j < numCols; j++) {
            matrix[i][j] = inputArr.at(i).at(j);
        }
    }
}

Matrix::Matrix(double a, double b, double c, double d) {
    numRows = 2;
    numCols = 2;
    matrix.resize(2, vector<double>(2));
    matrix[0][0] = a;
    matrix[0][1] = b;
    matrix[1][0] = c;
    matrix[1][1] = d;
}

Matrix::Matrix(double a, double b) {
    numRows = 2;
    numCols = 1;
    matrix.resize(2, vector<double>(1));
    matrix[0][0] = a;
    matrix[1][0] = b;
}

int Matrix::getNumRows() {
    return numRows;
}

int Matrix::getNumCols() {
    return numCols;
}

double Matrix::get(int i, int j) {
    return matrix.at(i).at(j);
}

double Matrix::determinant() {
    if (numRows != 2 || numCols != 2) {
        throw domain_error("unsupported dimensions");
    }
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
}

Matrix Matrix::adjugate() {
    if (numRows != 2 || numCols != 2) {
        throw domain_error("unsupported dimensions");
    }
    return Matrix(matrix[1][1], -matrix[0][1], -matrix[1][0], matrix[0][0]);
}

Matrix Matrix::inverse() {
    if (numRows != 2 || numCols != 2) {
        throw domain_error("unsupported dimensions");
    }
    double det = determinant();
    if (det == 0) {
        throw domain_error("singular matrix");
    }
    return Matrix(matrix[1][1]/det, -matrix[0][1]/det, -matrix[1][0]/det, matrix[0][0]/det);
}

Matrix Matrix::multiply(Matrix &other) {
    if(numCols != other.getNumRows()) {
        throw domain_error("incompatible dimensions");
    }
    int rows = numRows;
    int cols = other.getNumCols();
    vector<vector<double>> product;
    product.resize(rows, vector<double>(cols));
    for(int i = 0; i < rows; i++) {
        for(int j = 0; j < cols; j++) {
            double dotProd = 0;
            for(int k = 0; k < numCols; k++) {
                dotProd += matrix.at(i).at(k) * other.matrix.at(k).at(j);
            }
            product[i][j] = dotProd;
        }
    }
    return Matrix(product);
}

Matrix Matrix::transpose() {
    vector<vector<double>> transposeArr(numCols, vector<double>(numRows));
    for(int i = 0; i < numRows; i++) {
        for(int j = 0; j < numCols; j++) {
            transposeArr[j][i] = matrix.at(i).at(j);
        }
    }
    return Matrix(transposeArr);
}

bool Matrix::operator==(const Matrix &other) const {
    if(numRows == other.numRows && numCols == other.numCols) {
        for(int i = 0; i < numRows; i++) {
            for(int j = 0; j < numCols; j++) {
                if(matrix.at(i).at(j) != other.matrix.at(i).at(j)) {
                    return false;
                }
            }
        }
        return true;
    }
    return false;
}

bool Matrix::operator!=(const Matrix &other) const {
    return !(*this == other);
}
