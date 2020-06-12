#include <gtest/gtest.h>
#include <math.h>
#include "../src/triangle.hpp"
#include "../src/lineIntegral.hpp"
#include "../src/doubleIntegral.hpp"
#include "../src/matrix.hpp"
using namespace std;

const double eps = 0.001;
const int num_iter = 1000;
const double tolerance = 1e-6;

// for testing purposes, create an 100 x 100 image all of whose pixels have color 1
vector<vector<Pixel>> generateFakeImage() {
    vector<vector<Pixel>> image;
    for(int i = 0; i < 100; i++) {
        vector<Pixel> holder;
        for(int j = 0; j < 100; j++) {
            holder.push_back(Pixel(i, j, 1));
        }
        image.push_back(holder);
    }
    return image;
}

TEST(DerivativeTest, FixedTriangle) {
    double vx = 2 * ((double) rand() / RAND_MAX) - 1;
    // pick unit vector
    double vy = pow((1 - vx * vx), 0.5);
    // randomly pick sign
    if ((double) rand() / RAND_MAX < 0.5) {
        vy *= -1;
    }
    Point a(0,0);
    Point b(5,0);
    Point c(3,4);
    Triangle triangle(&a, &b, &c);
    double gradX = triangle.gradX(&c);
    double gradY = triangle.gradY(&c);
    // simulate eps timestep
    double curArea = triangle.getSignedArea();
    c.move(eps * vx, eps * vy);
    double futureArea = triangle.getSignedArea();
    c.move(-2 * eps * vx, -2 * eps * vy);
    double pastArea = triangle.getSignedArea();
    double change1 = (futureArea - curArea) / eps;
    double change2 = (curArea - pastArea) / eps;
    ASSERT_TRUE(abs(change1 - change2) < tolerance);
    double finiteDiff = (futureArea - pastArea) / (2 * eps);
    double gradApprox = (gradX * vx + gradY * vy);
    ASSERT_TRUE(abs(gradApprox - finiteDiff) < tolerance);
}

TEST(DerivativeTest, TriangleArea) {
    for(int i = 0; i < num_iter; i++) {
        // pick a random number from -1 to 1 for velocity
        double vx = 2 * ((double) rand() / RAND_MAX) - 1;
        double vy = 2 * ((double) rand() / RAND_MAX) - 1;
        // generate some random triangle with coordinates < 100
        double random_coords[6];
        for(int j = 0; j < 6; j++) {
            random_coords[j] = (double) rand() / RAND_MAX * 100;
        }
        Point a(random_coords[0], random_coords[1]);
        Point b(random_coords[2], random_coords[3]);
        Point c(random_coords[4], random_coords[5]);
        Segment e(&a, &b);
        Triangle triangle(&a, &b, &c);
        double gradX = triangle.gradX(&a);
        double gradY = triangle.gradY(&a);
        // simulate area(x + eps * v), will have to undo later
        a.move(eps * vx, eps * vy);
        double area1 = triangle.getSignedArea();
        a.move(- 2 * eps * vx, - 2 * eps * vy);
        double area2 = triangle.getSignedArea();
        double areaDiff = (area1 - area2) / (2 * eps);
        double gradApprox = gradX * vx + gradY * vy;
        ASSERT_TRUE(abs(areaDiff - gradApprox) < tolerance);
    }
}

// test derivative for energy function with constant approximation
// over a single triangle
TEST(DerivativeTest, ConstantApprox) {
    // for just integrating the given function
    auto identity = [](double x, double y) {return 1.0;};
    vector<vector<Pixel>> image = generateFakeImage();
    for(int i = 0; i < num_iter; i++) {
        cout << "iteration " << i << endl;
        // pick a random number from -1 to 1 for velocity directions
        // for testing, point a will be moved
        double vx = 2 * ((double) rand() / RAND_MAX) - 1;
        double vy = 2 * ((double) rand() / RAND_MAX) - 1;
        // generate some random triangle with coordinates < 100
        double random_coords[6];
        for(int j = 0; j < 6; j++) {
            random_coords[j] = (double) rand() / RAND_MAX * 100;
        }
        Point a(random_coords[0], random_coords[1]);
        Point b(random_coords[2], random_coords[3]);
        Point c(random_coords[4], random_coords[5]);
        Triangle triangle(&a, &b, &c);
        // for integrating v dot n in the line integral when a is moving
        // at a velocity of (1, 0);
        // (x, y) will be on either ab or ca
        auto vn = [&a, &b, &c](double x, double y, bool dirX) {
            Point current(x, y);
            // determine which segment current is on
            Point segmentEnd = (Triangle::getSignedArea(&a, &b, &current) == 0) ? b : c;
            // compute velocity at this point by scaling
            double distanceToVertex = current.distance(segmentEnd);
            cout << "distanceToVertex: " << distanceToVertex << endl;
            Segment gamma(&a, &segmentEnd);
            double scale = distanceToVertex / gamma.length(); // 1 if at a, 0 if at opposite edge
            cout << "scale factor " << scale << endl;
            // determine whether x or y gradient is being computed
            double velX = (dirX) ? scale : 0;
            Matrix v(velX, scale - velX);
            Matrix n = gamma.unitNormal();
            return (v.transpose().multiply(n)).get(0, 0);
        };
        // break down into x and y components
        auto vnx = [&vn](double x, double y) {
            return vn(x, y, true);
        };
        auto vny = [&vn](double x, double y) {
            return vn(x, y, false);
        };
        // integral of fdA
        double imageIntegral = DoubleIntegral::evaluate(identity, &image, &triangle);
        cout << "imageIntegral " << imageIntegral << endl;
        double area = triangle.getArea();
        cout << "area: " << area << endl;
        double dA[2] = {triangle.gradX(&a), triangle.gradY(&a)};
        cout << "dA: " << dA[0] << ", " << dA[1] << endl;
        double boundaryChange[2];
        // compute gradient in x direction
        boundaryChange[0] = LineIntegral::evaluate(vnx, &image, &triangle);
        // compute gradient in y direction
        boundaryChange[1] = LineIntegral::evaluate(vny, &image, &triangle);
        cout << "boundaryChange values " << boundaryChange[0] << ", " << boundaryChange[1] << endl;
        double gradient[2];
        for(int j = 0; j < 2; j++) {
            gradient[j] = (2 * area * imageIntegral * boundaryChange[j]
                - imageIntegral * imageIntegral * dA[j]) / (area * area);
        }
        cout << "gradient values: " << gradient[0] << ", " << gradient[1] << endl;
        double gradApprox = gradient[0] * vx + gradient[1] * vy;

        // now compute central finite difference
        a.move(eps * vx, eps * vy);
        double futureImgInt = DoubleIntegral::evaluate(identity, &image, &triangle);
        double futureArea = triangle.getArea();
        double futureEnergy = futureImgInt * futureImgInt / futureArea;

        a.move(-2 * eps * vx, -2 * eps * vy);
        double pastImgInt = DoubleIntegral::evaluate(identity, &image, &triangle);
        double pastArea = triangle.getArea();
        double pastEnergy = pastImgInt * pastImgInt / pastArea;

        double finiteApprox = (futureEnergy - pastEnergy) / (2 * eps);
        cout << "finiteApprox " << finiteApprox << endl;
        cout << "gradApprox " << gradApprox << endl;
        ASSERT_TRUE(abs(gradApprox - finiteApprox) < tolerance);
    }
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
