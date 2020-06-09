#include <gtest/gtest.h>
#include "../src/triangle.hpp"

using namespace std;

// test area of degenerate triangle
TEST(AreaTest, Degenerate) {
    Point a(0,0);
    Point b(0,1);
    Point c(0,2);
    Triangle triangle(&a, &b, &c);
    ASSERT_EQ(0, triangle.getArea());
}

TEST(AreaTest, RightTriangle) {
    Point a(0,0);
    Point b(0,3);
    Point c(4,0);
    Triangle triangle(&a, &b, &c);
    ASSERT_EQ(6, triangle.getArea());
}

// test that constructor re-orders in ccw direction
TEST(SignedAreaTest, Clockwise) {
    Point a(0,0);
    Point b(2,3);
    Point c(5,0);
    Triangle triangle(&a, &b, &c);
    ASSERT_EQ(7.5, triangle.getSignedArea());
}

TEST(SignedAreaTest, CounterClockwise) {
    Point a(0,0);
    Point b(5,0);
    Point c(3,4);
    Triangle triangle(&a, &b, &c);
    ASSERT_EQ(10, triangle.getSignedArea());
}

// test change in area when a triangle vertex moves
TEST(SignedAreaTest, Movement) {
    Point a(0,0);
    Point b(5,0);
    Point c(3,4);
    Triangle triangle(&a, &b, &c);
    c.move(0,-8);
    ASSERT_EQ(-10, triangle.getSignedArea());
    ASSERT_EQ(10, triangle.getArea());
}

// vertex moves parallel to opposite edge
TEST(DerivativeTest, Parallel) {
    Point a(0,0);
    Point b(5,0);
    Point c(3,4);
    Triangle triangle(&a, &b, &c);
    double change = triangle.dA(&c, 5, 0);
    ASSERT_EQ(0, change);
}

// vertex moves perpendicular to opposite dge
TEST(DerivativeTest, Perpendicular) {
    Point a(5,0);
    Point b(3,4);
    Point c(0,0);
    Triangle triangle(&a, &b, &c);
    double change = triangle.dA(&b, 0, 2);
    ASSERT_EQ(5, change);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
