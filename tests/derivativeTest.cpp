#include <gtest/gtest.h>
#include <math.h>
#include "../src/triangle.hpp"
using namespace std;

const double eps = 0.001;
const int num_iter = 1000;
const double tolerance = 1e-6;

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

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
