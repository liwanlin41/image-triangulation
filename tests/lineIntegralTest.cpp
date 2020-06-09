#include <gtest/gtest.h>
#include "../src/lineIntegral.hpp"
#include "../src/triangle.hpp"
#include "../src/point.hpp"

using namespace std;

// some test functions

double z(double x, double y) {
    return 0;
}

// constant function
double constant(double x, double y) {
    return 1;
}

// some non-constant function
double linear(double x, double y) {
    return x + y;
}

TEST(IntegralTest, Zero) {
    Point a(0,0);
    Point b(10,10);
    Point c(20,30);
    Triangle triangle(&a, &b, &c);
    ASSERT_EQ(0, LineIntegral::evaluate(z, &triangle));
    ASSERT_EQ(0, LineIntegral::evaluate(z, &a, &b, &c));
}

TEST(IntegralTest, Constant) {
    Point a(0,0);
    Point b(3,0);
    Point c(0,4);
    Triangle triangle(&a, &b, &c);
    ASSERT_EQ(12, LineIntegral::evaluate(constant, &triangle));
    ASSERT_EQ(12, LineIntegral::evaluate(constant, &a, &b, &c));
}

TEST(IntegralTest, Linear) {
    Point a(0,0);
    Point b(3,0);
    Point c(0,4);
    Triangle triangle(&a, &b, &c);
    ASSERT_EQ(30, LineIntegral::evaluate(linear, &triangle));
    ASSERT_EQ(30, LineIntegral::evaluate(linear, &a, &b, &c));
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}