#include <gtest/gtest.h>
#include "../src/parallelInt.h"

using namespace std;

TEST(CudaTest, SumTest) {
    int size = pow(2, 18);
    double arr[size]; // try a large array; there is an upper bound
    double sum = 0; // compute sum serially
    for(int i = 0; i < size; i++) {
        sum += i/2.0;
        arr[i] = i/2.0;
    }
    double res = runSum(arr, size);
    ASSERT_EQ(sum, res);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}