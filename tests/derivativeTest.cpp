#include <gtest/gtest.h>
#include "../src/triangle.hpp"
using namespace std;

const double eps = 0.001;
const int num_iter = 1000;
const double tolerance = 1e-6;

TEST(DerivativeTest, TriangleArea) {
    for(int i = 0; i < num_iter; i++) {
        // pick a random number from -1 to 1 for velocity
        double vx = 2 * ((double) rand() / RAND_MAX) - 1;
        double vy = 2 * ((double) rand() / RAND_MAX) - 1;
        // genreate some random triangle
        double random_coords[6];
        for(int j = 0; j < 6; j++) {
            random_coords[j] = rand();
        }
        Point a(random_coords[0], random_coords[1]);
        Point b(random_coords[2], random_coords[3]);
        Point c(random_coords[4], random_coords[5]);
        Triangle triangle(&a, &b, &c);
        double gradx = triangle.dA(&a, vx, 0);
        double grady = triangle.dA(&a, 0, vy);
        // simulate area(x + eps * v), will have to undo later
        a.move(eps * vx, eps * vy);
        double area1 = triangle.getSignedArea();
        a.move(- 2 * eps * vx, - 2 * eps * vy);
        double area2 = triangle.getSignedArea();
        double areaDiff = (area1 - area2) / (2 * eps);
        double gradApprox = gradx * vx + grady * vy;
        ASSERT_TRUE(abs(areaDiff - gradApprox) < tolerance);
    }
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}