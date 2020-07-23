// test gradient computations

#include <gtest/gtest.h>
#include "CImg.h"
#include "approx.h"
#include "constant.h"
#include "linear.h"

using namespace std;
using namespace cimg_library;

static const int TEST_ITERATIONS = 1000;
static const double eps = 0.001; // central finite difference step
static const double TOLERANCE = 0.02; // approximation ratio

static const int white = 255;
static const int black = 0;

static const int dimX = 100;
static const int dimY = 100;

static const double step = 0.5; // starting step value
static const double ds = 0.01;

// edge connections for a single triangle
static vector<array<int, 3>> faces = {{0,1,2}};

// randomly color an image in black and white
CImg<unsigned char> randomColor() {
    CImg<unsigned char> image(dimX, dimY, 1, 1, 0);
    for(int i = 0; i < dimX; i++) {
        for(int j = 0; j < dimY; j++) {
            image(i, j, 0, 0) = (rand() % 256);
        }
    }
    return image;
}

// color pixel (i,j) in color i (grayscale)
CImg<unsigned char> gradientColor() {
    CImg<unsigned char> image(dimX, dimY, 1, 1, 0);
    for(int i = 0; i < dimX; i++) {
        for(int j = 0; j < dimY; j++) {
            image(i, j, 0, 0) = i;
        }
    }
    return image;
}

CImg<unsigned char> constantColor(int c) {
    CImg<unsigned char> image(dimX, dimY, 1, 1, c);
    return image;
}

// convert CImg to pixel array
void convertToPixel(Pixel* &pixArr, CImg<unsigned char> *img) {
    cudaMallocManaged(&pixArr, dimX * dimY * sizeof(Pixel));
    for(int i = 0; i < dimX; i++) {
        for(int j = 0; j < dimY; j++) {
            pixArr[i * dimY + j] = Pixel(i, j, (*img)(i, j, 0, 0));
        }
    }
}

// due to approximation error, check if expected and computed
// are within a satisfactory approximation ratio
bool approxEqual(double expected, double computed, double tolerance = TOLERANCE) {
    return (abs(expected - computed) < tolerance * max(abs(expected), 1.0));
}

// however, note that the approximation error is highly variable
// and realistically anything within 10% is a decent approximation


TEST(GradientTest, ConstantBug) {
    CImg<unsigned char> image = constantColor(white);
    Pixel *pixArr;
    convertToPixel(pixArr, &image);
    ParallelIntegrator integrator; // external integrator
    double ds = 0.005;
    long long space = dimX / ds * dimY / ds + 1;
    integrator.initialize(pixArr, dimX, dimY, constant, space);
    Point *pts;
    cudaMallocManaged(&pts, 3 * sizeof(Point));
    Triangle tri(pts, pts + 1, pts + 2);

    ConstantApprox approx(&image, step, ds); // note step value doesn't really matter here

    // only test over a single triangle
    vector<Point> points;
    // generate random coordinates
    double coords[6] = {31.98845171, 88.80946765, 46.13423158, 9.284511394, 23.071969, 87.76759507}; 
    for(int i = 0; i < 3; i++) {
        points.push_back(Point(coords[2 * i], coords[2 * i + 1]));
        pts[i] = points.at(i); // for external integrator
    }
    approx.initialize(points, faces);
    // pick a random number from -1 to 1 for velocity directions
    double vx = 0.239132; 
    double vy = -0.800479;
    double gradX;
    double gradY;
    approx.gradient(0, 0, &gradX, &gradY);
    double gradApprox = gradX * vx + gradY * vy;

    // compute central finite difference
    // note that approx.gradient only computes -d/dt int fg dA,
    // discarding the int f^2 term because the integral is the same over all triangles;
    // however this is not the case for a single triangle
    // Therefore compute energy as int (f-g)^2 dA - int f^2 dA
    // to account for the actual gradient we need to check
    pts[0].move(eps * vx, eps * vy);
    double futureColor = integrator.doubleIntEval(&tri, ds) / tri.getArea();
    double futureEnergy = integrator.constantEnergyEval(&tri, futureColor, ds)
        - (integrator.constantEnergyEval(&tri, 0, ds) + LOG_AREA_MULTIPLIER * log(tri.getArea() - AREA_THRESHOLD));
        
    pts[0].move(-2 * eps * vx, -2 * eps * vy);
    double pastColor = integrator.doubleIntEval(&tri, ds) / tri.getArea();
    double pastEnergy = integrator.constantEnergyEval(&tri, pastColor, ds)
        - (integrator.constantEnergyEval(&tri, 0, ds) + LOG_AREA_MULTIPLIER * log(tri.getArea() - AREA_THRESHOLD));

    double finiteApprox = (futureEnergy - pastEnergy) / (2 * eps);
    if (!approxEqual(finiteApprox, gradApprox, 0.1)) {
        cout << finiteApprox << ", " << gradApprox << endl;
    }
    EXPECT_TRUE(approxEqual(finiteApprox, gradApprox, 0.1));
    cudaFree(pts);
    cudaFree(pixArr);
}

TEST(GradientTest, ConstantCase) {
    CImg<unsigned char> image = constantColor(1);
    Pixel *pixArr;
    convertToPixel(pixArr, &image);
    ParallelIntegrator integrator; // external integrator
    long long space = dimX / ds * dimY / ds + 1;
    integrator.initialize(pixArr, dimX, dimY, constant, space);
    Point *pts;
    cudaMallocManaged(&pts, 3 * sizeof(Point));
    Triangle tri(pts, pts + 1, pts + 2);

    for(int ii = 0; ii < TEST_ITERATIONS; ii++) {
        cout << "iteration " << ii << endl;
        ConstantApprox approx(&image, step, ds); // note step value doesn't really matter here

        // only test over a single triangle
        vector<Point> points;
        // generate random coordinates
        double coords[6] = {1,1,5,1,1,4};
        for(int i = 0; i < 6; i++) {
            // transform so that point lies well within image
            coords[i] = (double) rand() / RAND_MAX * (dimX - 3) + 1;
        }
        
        for(int i = 0; i < 3; i++) {
            points.push_back(Point(coords[2 * i], coords[2 * i + 1]));
            pts[i] = points.at(i); // for external integrator
        }
        approx.initialize(points, faces);
        // pick a random number from -1 to 1 for velocity directions
        double vx = 2 * ((double) rand() / RAND_MAX) - 1;
        double vy = 2 * ((double) rand() / RAND_MAX) - 1;
        double gradX;
        double gradY;
        approx.gradient(0, 0, &gradX, &gradY);
        double gradApprox = gradX * vx + gradY * vy;

        // compute central finite difference
        // note that approx.gradient only computes -d/dt int fg dA,
        // discarding the int f^2 term because the integral is the same over all triangles;
        // however this is not the case for a single triangle
        // Therefore compute energy as int (f-g)^2 dA - int f^2 dA
        // to account for the actual gradient we need to check
        pts[0].move(eps * vx, eps * vy);
        double futureColor = integrator.doubleIntEval(&tri, ds) / tri.getArea();
        double futureEnergy = integrator.constantEnergyEval(&tri, futureColor, ds)
            - (integrator.constantEnergyEval(&tri, 0, ds) + LOG_AREA_MULTIPLIER * log(tri.getArea() - AREA_THRESHOLD));
        
        pts[0].move(-2 * eps * vx, -2 * eps * vy);
        double pastColor = integrator.doubleIntEval(&tri, ds) / tri.getArea();
        double pastEnergy = integrator.constantEnergyEval(&tri, pastColor, ds)
            - (integrator.constantEnergyEval(&tri, 0, ds) + LOG_AREA_MULTIPLIER * log(tri.getArea() - AREA_THRESHOLD));

        double finiteApprox = (futureEnergy - pastEnergy) / (2 * eps);
        if (!approxEqual(finiteApprox, gradApprox)) {
            cout << finiteApprox << ", " << gradApprox << endl;
            pts[0].move(eps * vx, eps * vy);
            cout << tri;
            printf("vx, vy: %f, %f\n", vx, vy);
        }
        EXPECT_TRUE(approxEqual(finiteApprox, gradApprox, 0.1));
    }
    cudaFree(pts);
    cudaFree(pixArr);
}

// with smaller step size (decrease from 0.01 to 0.005), this passes
TEST(GradientTest, LinearBug1) {
    CImg<unsigned char> image = gradientColor();
    Pixel *pixArr;
    convertToPixel(pixArr, &image);
    ParallelIntegrator integrator; // external integrator
    double ds = 0.005;
    long long space = dimX / ds * dimY / ds + 1;
    integrator.initialize(pixArr, dimX, dimY, linear, space);
    Point *pts;
    cudaMallocManaged(&pts, 3 * sizeof(Point));
    Triangle tri(pts, pts + 1, pts + 2);

    LinearApprox approx(&image, step, ds); // note step value doesn't really matter here

    // only test over a single triangle
    vector<Point> points;
    // generate random coordinates

    double coords[6] = {87.58614578, 10.26248323, 89.98234349, 92.82157702, 51.00707523, 83.3485548}; 
    for(int i = 0; i < 3; i++) {
        points.push_back(Point(coords[2 * i], coords[2 * i + 1]));
        pts[i] = points.at(i); // for external integrator
    }
    // account for orientation
    if(Triangle::getSignedArea(pts, pts + 1, pts + 2) < 0) {
        cout << "backwards" << endl;
        pts[1] = points.at(2);
        pts[2] = points.at(1);
        points[1] = pts[1];
        points[2] = pts[2];
    }
    approx.initialize(points, faces);
    double vx = -0.205279;
    double vy = -0.201196;
    double gradX;
    double gradY;
    approx.gradient(0, 0, &gradX, &gradY);
    double gradApprox = gradX * vx + gradY * vy;

    // compute central finite difference
    // note that approx.gradient only computes -d/dt int fg dA,
    // discarding the int f^2 term because the integral is the same over all triangles;
    // however this is not the case for a single triangle
    // Therefore compute energy as int (f-g)^2 dA - int f^2 dA
    // to account for the actual gradient we need to check
    pts[0].move(eps * vx, eps * vy);
    double fCoeff[3]; // future coefficients
    approx.computeCoeffs(&tri, fCoeff);
    double futureEnergy = integrator.linearEnergyApprox(&tri, fCoeff, ds)
        - (integrator.constantEnergyEval(&tri, 0, ds) + LOG_AREA_MULTIPLIER * log(tri.getArea() - AREA_THRESHOLD));
        
    pts[0].move(-2 * eps * vx, -2 * eps * vy);
    double pCoeff[3]; // past coefficients
    approx.computeCoeffs(&tri, pCoeff);
    double pastEnergy = integrator.linearEnergyApprox(&tri, pCoeff, ds)
        - (integrator.constantEnergyEval(&tri, 0, ds) + LOG_AREA_MULTIPLIER * log(tri.getArea() - AREA_THRESHOLD));

    double finiteApprox = (futureEnergy - pastEnergy) / (2 * eps);
    EXPECT_TRUE(approxEqual(finiteApprox, gradApprox, 0.1));
    cudaFree(pts);
    cudaFree(pixArr);
}

// same here; decreasing ds from 0.01 to 0.003 passes test
TEST(GradientTest, LinearBug2) {
    CImg<unsigned char> image = gradientColor();
    Pixel *pixArr;
    convertToPixel(pixArr, &image);
    ParallelIntegrator integrator; // external integrator
    double ds = 0.003;
    long long space = dimX / ds * dimY / ds + 1;
    integrator.initialize(pixArr, dimX, dimY, linear, space);
    Point *pts;
    cudaMallocManaged(&pts, 3 * sizeof(Point));
    Triangle tri(pts, pts + 1, pts + 2);

    LinearApprox approx(&image, step, ds); // note step value doesn't really matter here

    // only test over a single triangle
    vector<Point> points;
    // generate random coordinates

    double coords[6] = {20.94301549, 96.00313661, 83.7031769, 29.76096179, 97.34614442, 52.41566814}; 
    for(int i = 0; i < 3; i++) {
        points.push_back(Point(coords[2 * i], coords[2 * i + 1]));
        pts[i] = points.at(i); // for external integrator
    }
    // account for orientation
    if(Triangle::getSignedArea(pts, pts + 1, pts + 2) < 0) {
        cout << "backwards" << endl;
        pts[1] = points.at(2);
        pts[2] = points.at(1);
        points[1] = pts[1];
        points[2] = pts[2];
    }
    approx.initialize(points, faces);
    double vx = -0.812471; 
    double vy = -0.270102;
    double gradX;
    double gradY;
    approx.gradient(0, 0, &gradX, &gradY);
    double gradApprox = gradX * vx + gradY * vy;

    // compute central finite difference
    // note that approx.gradient only computes -d/dt int fg dA,
    // discarding the int f^2 term because the integral is the same over all triangles;
    // however this is not the case for a single triangle
    // Therefore compute energy as int (f-g)^2 dA - int f^2 dA
    // to account for the actual gradient we need to check
    pts[0].move(eps * vx, eps * vy);
    double fCoeff[3]; // future coefficients
    approx.computeCoeffs(&tri, fCoeff);
    double futureEnergy = integrator.linearEnergyApprox(&tri, fCoeff, ds)
        - (integrator.constantEnergyEval(&tri, 0, ds) + LOG_AREA_MULTIPLIER * log(tri.getArea() - AREA_THRESHOLD));
        
    pts[0].move(-2 * eps * vx, -2 * eps * vy);
    double pCoeff[3]; // past coefficients
    approx.computeCoeffs(&tri, pCoeff);
    double pastEnergy = integrator.linearEnergyApprox(&tri, pCoeff, ds)
        - (integrator.constantEnergyEval(&tri, 0, ds) + LOG_AREA_MULTIPLIER * log(tri.getArea() - AREA_THRESHOLD));

    double finiteApprox = (futureEnergy - pastEnergy) / (2 * eps);
    EXPECT_TRUE(approxEqual(finiteApprox, gradApprox, 0.2));
    cudaFree(pts);
    cudaFree(pixArr);
}

TEST(GradientTest, LinearCase) {
    CImg<unsigned char> image = gradientColor();
    Pixel *pixArr;
    convertToPixel(pixArr, &image);
    ParallelIntegrator integrator; // external integrator
    long long space = dimX / ds * dimY / ds + 1;
    integrator.initialize(pixArr, dimX, dimY, linear, space);
    Point *pts;
    cudaMallocManaged(&pts, 3 * sizeof(Point));
    Triangle tri(pts, pts + 1, pts + 2);

    // discrepancy between finite difference and gradient seems
    // to be due to approximation error; in particular increasing sampling
    // rate seems to solve the problem for the extracted failing cases
    // and finite difference is highly susceptible to approximation error
    // (gradient is fairly stable)
    // to this end, check that the number of failing cases is relatively small
    int numBadApprox = 0;

    for(int ii = 0; ii < TEST_ITERATIONS; ii++) {
        cout << "iteration " << ii << endl;
        LinearApprox approx(&image, step, ds); // note step value doesn't really matter here

        // only test over a single triangle
        vector<Point> points;
        // generate random coordinates
        double coords[6];
        for(int i = 0; i < 6; i++) {
            // transform so that point lies well within image
            coords[i] = (double) rand() / RAND_MAX * (dimX - 3) + 1;
        }
        
        for(int i = 0; i < 3; i++) {
            points.push_back(Point(coords[2 * i], coords[2 * i + 1]));
            pts[i] = points.at(i); // for external integrator
        }
        // account for orientation
        if(Triangle::getSignedArea(pts, pts + 1, pts + 2) < 0) {
            pts[1] = points.at(2);
            pts[2] = points.at(1);
            points[1] = pts[1];
            points[2] = pts[2];
        }
        approx.initialize(points, faces);
        // pick a random number from -1 to 1 for velocity directions
        double vx = 2 * ((double) rand() / RAND_MAX) - 1;
        double vy = 2 * ((double) rand() / RAND_MAX) - 1;
        double gradX;
        double gradY;
        approx.gradient(0, 0, &gradX, &gradY);
        double gradApprox = gradX * vx + gradY * vy;

        // compute central finite difference
        // note that approx.gradient only computes -d/dt int fg dA,
        // discarding the int f^2 term because the integral is the same over all triangles;
        // however this is not the case for a single triangle
        // Therefore compute energy as int (f-g)^2 dA - int f^2 dA
        // to account for the actual gradient we need to check
        pts[0].move(eps * vx, eps * vy);
        double fCoeff[3]; // future coefficients
        approx.computeCoeffs(&tri, fCoeff);
        double futureEnergy = integrator.linearEnergyApprox(&tri, fCoeff, ds)
            - (integrator.constantEnergyEval(&tri, 0, ds) + LOG_AREA_MULTIPLIER * log(tri.getArea() - AREA_THRESHOLD));
        
        pts[0].move(-2 * eps * vx, -2 * eps * vy);
        double pCoeff[3]; // past coefficients
        approx.computeCoeffs(&tri, pCoeff);
        double pastEnergy = integrator.linearEnergyApprox(&tri, pCoeff, ds)
            - (integrator.constantEnergyEval(&tri, 0, ds) + LOG_AREA_MULTIPLIER * log(tri.getArea() - AREA_THRESHOLD));

        double finiteApprox = (futureEnergy - pastEnergy) / (2 * eps);
        /*
        if (!approxEqual(finiteApprox, gradApprox)) {
            cout << finiteApprox << ", " << gradApprox << endl;
            pts[0].move(eps * vx, eps * vy);
            cout << tri;
            printf("vx, vy: %f, %f\n", vx, vy);
        }
        */
        if(!approxEqual(finiteApprox, gradApprox, 5 * TOLERANCE)) {
            numBadApprox++;
        }
    }
    cout << "number of poor approximations: " << numBadApprox << endl;
    ASSERT_TRUE(numBadApprox < TOLERANCE * TEST_ITERATIONS);
    cudaFree(pts);
    cudaFree(pixArr);
}