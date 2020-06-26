#include <gtest/gtest.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "../src/parallelInt.cuh"

using namespace std;

/*
TEST(CudaTest, SumTest) {
    int size = pow(2, 18);
    double arr[size]; // try a large array; there is an upper bound
    double sum = 0; // compute sum serially
    for(int i = 0; i < size; i++) {
        sum += i/2.0;
        arr[i] = i/2.0;
    }
    //double res = runSum(arr, size);
    //ASSERT_EQ(sum, res);
}
*/

TEST(LineTest, Triangle) {
    ApproxType APPROXTYPE = constant;
    int maxX = 5;
    int maxY = 4;
    Pixel *pixArr;
    double *results;
    cudaMallocManaged(&pixArr, 20 * sizeof(Pixel));
    cudaMallocManaged(&results, 20 * sizeof(double));
    for(int i = 0; i < 5; i++) {
        for(int j = 0; j < 4; j++) {
            pixArr[4*i+j] = Pixel(i,j,1);
        }
    }
    Point *points;
    cudaMallocManaged(&points, 3*sizeof(Point));
    points[0] = Point(0,0,true, true);
    points[1] = Point(4,0,true, true);
    points[2] = Point(0,3,true, true);
    Triangle *tri;
    Point *working;
    cudaMallocManaged(&tri, sizeof(Triangle));
    cudaMallocManaged(&working, 3*sizeof(Point));
    tri[0] = Triangle(points, points+1, points+2);
    int t = 0;
    int pt = 2;
    double ans = lineIntEval(APPROXTYPE, pixArr, maxX, maxY, tri, t, pt, true, results, working);
    EXPECT_EQ(0, ans);
    cudaFree(pixArr);
    cudaFree(results);
    cudaFree(points);
    cudaFree(tri);
    cudaFree(working);
}

TEST(AreaTest, Simple) {
    ApproxType testApprox = constant;
    int maxX = 5;
    int maxY = 4;
    Pixel *pixArr;
    Triangle *triArr;
    Point *points;
    int t = 0;
    cudaMallocManaged(&pixArr, maxX * maxY*sizeof(Pixel));
    cudaMallocManaged(&triArr, sizeof(Triangle));
    cudaMallocManaged(&points, 3*sizeof(Point));
    for(int i = 0; i < maxX; i++) {
        for(int j = 0; j < maxY; j++) {
            pixArr[i * maxY + j] = Pixel(i, j, 1);
        }
    }
    points[0] = Point(0,0);
    points[1] = Point(4,0);
    points[2] = Point(0,3);
    triArr[0] = Triangle(points, points+1, points+2);
    double *results;
    cudaMallocManaged(&results, maxX * maxY * sizeof(double));
    double test = doubleIntEval(testApprox, pixArr, maxX, maxY, triArr, t, results);
    cout << test << endl;
    cudaFree(pixArr);
    cudaFree(triArr);
    cudaFree(points);
    cudaFree(results);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}