#include <gtest/gtest.h>
#include "../src/parallelInt.h"

TEST(CudaTest, SumTest) {
    int size = 16;
    double arr[size]; // try a large array
    double sum = 0; // compute sum serially
    for(int i = 0; i < size; i++) {
        sum += i;
        arr[i] = i;
    }
    double res = runSum(arr, size);
    ASSERT_EQ(sum, res);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}