#include <gtest/gtest.h>
#include "../headers/matrix.hpp"

// test determinant of singular matrix
TEST(DeterminantTest, Zero) {
    Matrix matrix(1,2,2,4);
    ASSERT_EQ(0, matrix.determinant());
}

// test determinant of full rank matrix
TEST(DeterminantTest, FullRank) {
    double row1[] = {1/3, 1};
    double row2[] = {0, 3};
    Matrix matrix(row1, row2);
    ASSERT_EQ(1, matrix.determinant());
}

TEST(AdjugateTest, ArrConstructor) {
    double arr[2][2] = {{3,4}, {100,-3.4}};
    Matrix matrix(arr);
    double expected[2][2] = {{-3.4,-4}, {-100, 3}};
    Matrix result = matrix.adjugate();
    for(int i = 0; i < 2; i++) {
        for(int j = 0; j < 2; j++) {
            ASSERT_EQ(expected[i][j], result.get(i,j));
        }
    }
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}