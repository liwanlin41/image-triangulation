// test Integrator

#include <gtest/gtest.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <vector>
#include "parallelInt.cuh"

using namespace std;

static const int white = 255;
static const int black = 0;

// dimensions of test square
static const int dimX = 100;
static const int dimY = 100;

// degree of tolerance (multiplicative)
static const double TOLERANCE = 0.01;
static const double ds = 0.01;
// for additional testing, uncomment the cout lines
// and check that the values of exact and approximate
// integrals get closer as ds -> 0
// note that for area integrals, ds cannot be too small
// or else floating point error becomes a factor
// (particularly for linear approximation - scales cubically)

// given an uninitialized array pointer pixArr,
// make pixArr represent a 100x100 solid square with
// a given color
void setSolidColor(Pixel* &pixArr, int color) {
    cudaMallocManaged(&pixArr, dimX * dimY * sizeof(Pixel));
    for(int i = 0; i < dimX; i++) {
        for(int j = 0; j < dimY; j++) {
            pixArr[i * dimY + j] = Pixel(i, j, color);
        }
    }
}

// generate random 100x100 (grayscale) image
void setRandomColor(Pixel* &pixArr) {
    cudaMallocManaged(&pixArr, dimX * dimY * sizeof(Pixel));
    for(int i = 0; i < dimX; i++) {
        for(int j = 0; j < dimY; j++) {
            int color = rand() % 256;
            pixArr[i * dimY + j] = Pixel(i, j, color);
        }
    }
}

// color (x,y) by x
void setGradientColor(Pixel* &pixArr) {
    cudaMallocManaged(&pixArr, dimX * dimY * sizeof(Pixel));
    for(int i = 0; i < dimX; i++) {
        for(int j = 0; j < dimY; j++) {
            pixArr[i * dimY + j] = Pixel(i, j, i);
        }
    }
}

// evaluate the exact integral of f dA
TEST(ExactIntegralTest, DoubleInt) {
    Pixel* pixArr;
    ParallelIntegrator integrator;
    setSolidColor(pixArr, white);
    integrator.initialize(pixArr, dimX, dimY, constant, dimX * dimY, true);

    Point *points;
    cudaMallocManaged(&points, 4 * sizeof(Point));

    points[0] = Point(0,0);
    points[1] = Point(4,0);
    points[2] = Point(0,3);
    points[3] = Point(5,1);
    Triangle rightTri(points, points + 1, points + 2);
    Triangle tri(points + 1, points + 2, points + 3);
    ASSERT_FLOAT_EQ(1530, integrator.doubleIntExact(&rightTri, GRAY));
    ASSERT_FLOAT_EQ(892.5, integrator.doubleIntExact(&tri, GRAY));
    cudaFree(pixArr);
    cudaFree(points);
}

// evaluate exact energy integral
TEST(ExactIntegralTest, EnergyInt) {
    Pixel* pixArr;
    ParallelIntegrator integrator;
    setSolidColor(pixArr, black);
    integrator.initialize(pixArr, dimX, dimY, constant, dimX * dimY, true);
    
    Point* points;
    cudaMallocManaged(&points, 3 * sizeof(Point));
    points[0] = Point(10,10);
    points[1] = Point(15,15);
    points[2] = Point(10,20);
    Triangle tri(points, points + 1, points + 2);
    double area = tri.getArea();
    ASSERT_EQ(25, area);
    int color = 10; // approximation color to test
    ASSERT_EQ(area * color * color - LOG_AREA_MULTIPLIER * log(area), integrator.constantEnergyExact(&tri, color));
    cudaFree(pixArr);
    cudaFree(points);
}

// evaluate exact line integral of v * n f dA
TEST(ExactIntegralTest, LineInt) {
    Pixel *pixArr;
    ParallelIntegrator integrator;
    setSolidColor(pixArr, white); 
    integrator.initialize(pixArr, dimX, dimY, constant, dimX * dimY, true);

    Point *points;
    cudaMallocManaged(&points, 3 * sizeof(Point));
    points[0] = Point(0,0);
    points[1] = Point(3,0);
    points[2] = Point(0,4);
    Triangle tri(points, points + 1, points + 2);
    ASSERT_EQ(-2 * white, integrator.lineIntExact(&tri, 0, true));
    ASSERT_DOUBLE_EQ(-1.5 * white, integrator.lineIntExact(&tri, 0, false));
    cudaFree(pixArr);
    cudaFree(points);
}

// test approximations
TEST(ApproximationTest, DoubleInt) {
    Pixel *pixArr;
    ParallelIntegrator integrator;
    setRandomColor(pixArr);
    long long space = dimX / ds * dimY / ds + 1;
    integrator.initialize(pixArr, dimX, dimY, constant, space);

    Point *points;
    cudaMallocManaged(&points, 3 * sizeof(Point));
    points[0] = Point(0,0);
    points[1] = Point(4,2);
    points[2] = Point(1,3);
    Triangle tri(points, points + 1, points + 2);

    double exact = integrator.doubleIntExact(&tri, GRAY);
    double approx = integrator.doubleIntApprox(&tri, ds, GRAY, 0);
    // cout << exact << ", " << approx << endl;
    ASSERT_TRUE(abs(exact - approx) < TOLERANCE * abs(exact));
    cudaFree(pixArr);
    cudaFree(points);
}

TEST(ApproximationTest, EnergyInt) {
    Pixel *pixArr;
    ParallelIntegrator integrator;
    setRandomColor(pixArr);
    long long space = dimX / ds * dimY / ds + 1;
    integrator.initialize(pixArr, dimX, dimY, constant, space);
    Point *points;
    cudaMallocManaged(&points, 3 * sizeof(Point));
    points[0] = Point(1,1);
    points[1] = Point(5,3);
    points[2] = Point(2,4);
    Triangle tri(points, points + 1, points + 2);

    const int color = integrator.doubleIntEval(&tri, ds, GRAY) / tri.getArea();
    double exact = integrator.constantEnergyExact(&tri, color);
    double approx = integrator.constantEnergyApprox(&tri, color, ds);
    // cout << exact << ", " << approx << endl;
    ASSERT_TRUE(abs(exact - approx) < TOLERANCE * abs(exact));
    cudaFree(pixArr);
    cudaFree(points);
}

TEST(ApproximationTest, LineInt) {
    Pixel *pixArr;
    ParallelIntegrator integrator;
    setRandomColor(pixArr);
    long long space = dimX / ds * dimY / ds + 1;
    integrator.initialize(pixArr, dimX, dimY, constant, space);
    Point *points;
    cudaMallocManaged(&points, 3 * sizeof(Point));
    points[0] = Point(1,1);
    points[1] = Point(5,3);
    points[2] = Point(2,4);
    Triangle tri(points, points + 1, points + 2);

    // compute in both x and y directions
    for(int i = 0; i < 2; i++) {
        double exact = integrator.lineIntExact(&tri, 0, (i == 0));
        double approx = integrator.lineIntApprox(&tri, 0, (i == 0), ds, 0);
        //cout << exact << ", " << approx << endl;
        ASSERT_TRUE(abs(exact - approx) < TOLERANCE * abs(exact));
    }
    cudaFree(pixArr);
    cudaFree(points);
}

// test integrals for linear color approximation
TEST(LinearIntegrationTest, DoubleIntConstant) {
    Pixel *pixArr;
    ParallelIntegrator integrator;
    setSolidColor(pixArr, white);
    long long space = dimX / ds * dimY / ds + 1;
    integrator.initialize(pixArr, dimX, dimY, linear, space);
    Point *points;
    cudaMallocManaged(&points, 3 * sizeof(Point));
    points[0] = Point(0,0);
    points[1] = Point(4,0);
    points[2] = Point(0,3);
    Triangle tri(points, points + 1, points + 2);

    double expected[3] = {2 * white, 2 * white, 2 * white};
    for(int i = 0; i < 3; i++) {
        double approx = integrator.doubleIntEval(&tri, ds, GRAY, i);
        //cout << expected[i] << ", " << approx << endl;
        ASSERT_TRUE(abs(expected[i] - approx) < TOLERANCE * abs(expected[i]));
    }
    cudaFree(pixArr);
    cudaFree(points);
}

// test double integral with gradient coloring
TEST(LinearIntegrationTest, DoubleIntShaded) {
    Pixel *pixArr;
    ParallelIntegrator integrator;
    setGradientColor(pixArr);
    double ds = 0.001; // higher degree of accuracy needed
    long long space = dimX / ds * dimY / ds + 1;
    integrator.initialize(pixArr, dimX, dimY, linear, space);
    Point *points;
    cudaMallocManaged(&points, 3 * sizeof(Point));
    points[0] = Point(0,0);
    points[1] = Point(8,0);
    points[2] = Point(0,6);
    Triangle tri(points, points + 1, points + 2);

    double expected[3] = {16, 32, 16};
    for(int i = 0; i < 3; i++) {
        double approx = integrator.doubleIntEval(&tri, ds, GRAY, i);
        //cout << expected[i] << ", " << approx << endl;
        ASSERT_TRUE(abs(expected[i] - approx) < TOLERANCE * abs(expected[i]));
    }
    cudaFree(pixArr);
    cudaFree(points);
}

TEST(LinearIntegrationTest, EnergyInt) {
    Pixel *pixArr;
    ParallelIntegrator integrator;
    setGradientColor(pixArr);
    long long space = dimX / ds * dimY / ds + 1;
    integrator.initialize(pixArr, dimX, dimY, linear, space);
    Point *points;
    cudaMallocManaged(&points, 3 * sizeof(Point));
    points[0] = Point(0,0);
    points[1] = Point(12,0);
    points[2] = Point(0,9);
    Triangle tri(points, points + 1, points + 2);
    
    // coefficients to test; somewhat arbitrary
    double coeffs[3] = {0, 2, 0.5};
    double expected = 857.25;
    double computed = integrator.linearEnergyApprox(&tri, coeffs, ds) + LOG_AREA_MULTIPLIER * log(tri.getArea());
    // cout << expected << ", " << computed << endl;
    ASSERT_TRUE(abs(expected - computed) < TOLERANCE * abs(expected));
    cudaFree(pixArr);
    cudaFree(points);
}

TEST(LinearIntegrationTest, LowEnergy) {
    Pixel *pixArr;
    ParallelIntegrator integrator;
    setGradientColor(pixArr);
    long long space = dimX / ds * dimY / ds + 1;
    integrator.initialize(pixArr, dimX, dimY, linear, space);
    Point *points;
    cudaMallocManaged(&points, 3 * sizeof(Point));
    points[0] = Point(0,0);
    points[1] = Point(48,0);
    points[2] = Point(0,36);
    Triangle tri(points, points + 1, points + 2);

    // exact value
    double coeffs[3] = {0, 48, 0};
    // due to image discretization, the energy integral will end up
    // resembling the integral of (x-0.5)^2 which is 1/12 on a square
    double expected = tri.getArea() / 12;
    double computed = integrator.linearEnergyApprox(&tri, coeffs, ds) + LOG_AREA_MULTIPLIER * log(tri.getArea());
    // cout << expected << ", " << computed << endl;
    ASSERT_TRUE(abs(expected - computed) < TOLERANCE * abs(expected));
    cudaFree(pixArr);
    cudaFree(points);
}

TEST(LinearIntegrationTest, LineInt) {
    Pixel *pixArr;
    ParallelIntegrator integrator;
    setGradientColor(pixArr);
    long long space = dimX / ds * dimY / ds + 1;
    integrator.initialize(pixArr, dimX, dimY, linear, space);
    Point *points;
    cudaMallocManaged(&points, 3 * sizeof(Point));
    points[0] = Point(0,0);
    points[1] = Point(4,0);
    points[2] = Point(0,3);
    Triangle tri(points, points + 1, points + 2);

    double expectedX[3] = {0,3,1};
    // first handle x movement
    for(int i = 0; i < 3; i++) {
        double computed = integrator.lineIntEval(&tri, 1, true, ds, i);
        //cout << expectedX[i] << ", " << computed << endl;
        // add slight extra tolerance to handle expected value being 0;
        // additionally loosen bound for small segment size
        ASSERT_TRUE(abs(expectedX[i] - computed) < 2 * TOLERANCE * max(abs(expectedX[i]), 1.0));
    }
    // then y movement
    double expectedY[3] = {-4.0/3, 0, 4.0/3};
    for(int i = 0; i < 3; i++) {
        double computed = integrator.lineIntEval(&tri, 1, false, ds, i);
        //cout << expectedY[i] << ", " << computed << endl;
        ASSERT_TRUE(abs(expectedY[i] - computed) < 2 * TOLERANCE * max(abs(expectedY[i]), 1.0));
    }
    cudaFree(pixArr);
    cudaFree(points);
}

TEST(LinearIntegrationTest, ImageGradient) {
    Pixel *pixArr;
    ParallelIntegrator integrator;
    setGradientColor(pixArr);
    long long space = dimX / ds * dimY / ds + 1;
    integrator.initialize(pixArr, dimX, dimY, linear, space);
    Point *points;
    cudaMallocManaged(&points, 3 * sizeof(Point));
    points[0] = Point(0,0);
    points[1] = Point(4,0);
    points[2] = Point(0,3);
    Triangle tri(points, points + 1, points + 2);

    // animate points[1] in the x direction
    double expectedX[3] = {1, -1, 0};
    for(int i = 0; i < 3; i++) {
        double computed = integrator.linearImageGradient(&tri, 1, true, ds, i);
        //cout << expectedX[i] << ", " << computed << endl;
        ASSERT_TRUE(abs(expectedX[i] - computed) < TOLERANCE * max(abs(expectedX[i]), 1.0));
    }
    double expectedY[3] = {4.0/3, 0, -4.0/3};
    for(int i = 0; i < 3; i++) {
        double computed = integrator.linearImageGradient(&tri, 1, false, ds, i);
        //cout << expectedY[i] << ", " << computed << endl;
        ASSERT_TRUE(abs(expectedY[i] - computed) < TOLERANCE * max(abs(expectedY[i]), 1.0));
    }
    cudaFree(pixArr);
    cudaFree(points);
}