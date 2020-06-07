#ifndef matrix_h
#define matrix_h

#include <iostream>
#include <vector>
using namespace std;

/**
 * represent a rectangular matrix of arbitrary size,
 * with special support for 2x2 matrices
 */

class Matrix {
    private:
        int numRows, numCols;
        vector<vector<double>> matrix; // hold matrix values
    public:
        // constructor assuming the input vector is rectangular
        Matrix(vector<vector<double>> const &);
        // special constructor for 2x2 matrix, creating the matrix [[a,b],[c,d]]
        Matrix(double a, double b, double c, double d);

        int getNumRows();
        int getNumCols();
        // return the element in row i, column j
        double get(int i, int j);
        // return result of this @ other for compatible matrices
        Matrix multiply(Matrix &other);

        // compute determinant of 2x2 matrix
        double determinant();
        // compute adjugate of 2x2 matrix
        Matrix adjugate();
        // compute inverse of 2x2 matrix
        Matrix inverse();

        bool operator==(const Matrix &other) const;
        bool operator!=(const Matrix &other) const;
};

#endif