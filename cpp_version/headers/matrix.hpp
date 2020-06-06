#ifndef matrix_h
#define matrix_h

/**
 * represent a 2 x 2 square matrix
 */

class Matrix {
    private:
        double matrix[2][2];
    public:
        // create matrix a b
        //               c d
        Matrix(double a, double b, double c, double d);
        // create a matrix from a 2 by 2 array
        Matrix(double[2][2]);
        // create a matrix by row input
        Matrix(double row1[2], double row2[2]);
        
        // get the matrix element in row i, column j
        double get(int i, int j);

        double determinant();
        Matrix adjugate();
        Matrix inverse();
};

#endif