// just test compilation for gtest
#include <gtest/gtest.h>
#include "../src/point.cuh"
#include "../src/triangle.cuh"

using namespace std;

TEST(CompileTest, TriangleCompile) {
    Point a(0,0);
    Point b(4,0);
    Point c(0,3);
    Triangle tri(&a, &b, &c);
    ASSERT_EQ(6, tri.getArea());
}