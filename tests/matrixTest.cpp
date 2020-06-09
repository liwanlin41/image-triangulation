#include <gtest/gtest.h>
#include "../src/matrix.hpp"

// test determinant of singular matrix
TEST(DeterminantTest, Zero) {
    Matrix matrix(1,2,2,4);
    ASSERT_EQ(0, matrix.determinant());
}

// test determinant of full rank matrix
TEST(DeterminantTest, FullRank) {
    vector<vector<double>> arr = {{1.0/3, 1}, {0,3}};
    Matrix matrix(arr);
    ASSERT_EQ(1, matrix.determinant());
}

// from https://stackoverflow.com/questions/23270078
TEST(DeterminantTest, Error) {
    vector<vector<double>> arr = {{1,2,3},{4,5,6}};
    Matrix matrix(arr);
    EXPECT_THROW({
        try
        {
            matrix.determinant();
        }
        catch( const domain_error& e )
        {
            // and this tests that it has the correct message
            EXPECT_STREQ("unsupported dimensions", e.what() );
            throw;
        }
        }, domain_error);
}

TEST(AdjugateTest, Basic) {
    vector<vector<double>> arr = {{3,4}, {100,-3.4}};
    Matrix matrix(arr);
    double expected[2][2] = {{-3.4,-4}, {-100, 3}};
    Matrix result = matrix.adjugate();
    for(int i = 0; i < 2; i++) {
        for(int j = 0; j < 2; j++) {
            ASSERT_EQ(expected[i][j], result.get(i,j));
        }
    }
}

TEST(SizeTest, NonSquare) {
    vector<vector<double>> arr = {{1,2,3}, {0,0,0}, {4,5,6},{7,8,9}};
    Matrix matrix(arr);
    int numRows = matrix.getNumRows();
    int numCols = matrix.getNumCols();
    ASSERT_EQ(4, numRows);
    ASSERT_EQ(3, numCols);
}

TEST(MultiplyTest, Square) {
    // choose two non-commuting matrices
    Matrix matrix1(1,2,3,4);
    Matrix matrix2(1,0,-1,3);
    Matrix expected(-1,6,-1,12);
    Matrix product = matrix1.multiply(matrix2);
    ASSERT_EQ(expected, product);
}

// multiply two larger non-square matrices
TEST(MultiplyTest, NonSquare) {
    vector<vector<double>> arr1 = {{1, 1, 1}, {0, 1, 2}};
    vector<vector<double>> arr2 = {{1,0,0,3},{0,1,0,-1},{0,0,1,2}};
    vector<vector<double>> expectedArr = {{1,1,1,4},{0,1,2,3}};
    Matrix matrix1(arr1);
    Matrix matrix2(arr2);
    Matrix expected(expectedArr);
    Matrix product = matrix1.multiply(matrix2);
    ASSERT_EQ(2, product.getNumRows());
    ASSERT_EQ(4, product.getNumCols());
    ASSERT_EQ(expected, product);
}

TEST(MultiplyTest, Incompatible) {
    vector<vector<double>> arr1 = {{1,1,1},{0,1,2}};
    Matrix matrix1(arr1);
    Matrix matrix2(1,2,3,4);
    EXPECT_THROW({
        try {
            matrix1.multiply(matrix2);
        } catch (domain_error &e) {
            EXPECT_STREQ("incompatible dimensions", e.what());
            throw;
        }
    }, domain_error);
}

TEST(InverseTest, SpecialLinear) {
    Matrix identity(1,0,0,1);
    Matrix matrix(2,5,0,-0.5);
    Matrix inverse = matrix.inverse();
    ASSERT_EQ(identity, matrix.multiply(inverse));
    ASSERT_EQ(identity, inverse.multiply(matrix));
}

TEST(InverseTest, GeneralLinear) {
    Matrix identity(1,0,0,1);
    Matrix matrix(1,3,4,1);
    Matrix inverse = matrix.inverse();
    Matrix product = matrix.multiply(inverse);
    for(int i = 0; i < 2; i++) {
        for(int j = 0; j < 2; j++) {
            ASSERT_DOUBLE_EQ(identity.get(i,j), product.get(i,j));
        }
    }
}

TEST(InverseTest, Singular) {
    Matrix matrix(1,2,2,4);
    EXPECT_THROW({
        try {
            matrix.inverse();
        } catch (domain_error &e) {
            EXPECT_STREQ("singular matrix", e.what());
            throw;
        }
    }, domain_error);
}

TEST(EqualityTest, DimensionMismatch) {
    Matrix square(1,1,1,1);
    vector<vector<double>> arr = {{1,1,1,1}};
    Matrix row(arr);
    ASSERT_NE(row, square);
}

TEST(TransposeTest, NonSquare) {
    Matrix vec(1,0);
    Matrix transpose = vec.transpose();
    vector<vector<double>> arr = {{1,0}};
    Matrix expected(arr);
    ASSERT_EQ(expected, transpose);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
