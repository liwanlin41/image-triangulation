#include <gtest/gtest.h>
#include "../headers/triangle.hpp"

using namespace std;

// test area of degenerate triangle
TEST(AreaTest, Degenerate) {
    Point a(0,0);
    Point b(0,1);
    Point c(0,2);
    Triangle triangle(a, b, c);
    ASSERT_EQ(0, triangle.getArea());
}

TEST(AreaTest, RightTriangle) {
    Point a(0,0);
    Point b(0,3);
    Point c(4,0);
    Triangle triangle(a, b, c);
    ASSERT_EQ(6, triangle.getArea());
}

TEST(SignedAreaTest, Clockwise) {
    Point a(0,0);
    Point b(2,3);
    Point c(5,0);
    Triangle triangle(a, b, c);
    ASSERT_EQ(-7.5, triangle.getSignedArea());
}

TEST(SignedAreaTest, CounterClockwise) {
    Point a(0,0);
    Point b(5,0);
    Point c(3,4);
    Triangle triangle(a, b, c);
    ASSERT_EQ(10, triangle.getSignedArea());
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}