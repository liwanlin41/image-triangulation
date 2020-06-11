#include <gtest/gtest.h>
#include <math.h>
#include "../src/pixel.hpp"

using namespace std;

TEST(LengthTest, NoIntersect){
    Pixel p(0,0, 255);
    Point a(1,1);
    Point b(2,2);
    Segment e(&a, &b);
    ASSERT_EQ(0, p.intersectionLength(e));
}

// segment crosses entirety of pixel
TEST(LengthTest, CrossesPixel) {
    Pixel p(0,0,255);
    Point a(-0.5,-1);
    Point b(0.5,1);
    Segment e(&a, &b);
    // intersections should be (-0.25,-0.5) -- (0.25,0.5)
    ASSERT_DOUBLE_EQ(pow(1.25,0.5), p.intersectionLength(e));
}

// one endpoint of segment is inside pixel
TEST(LengthTest, PointInsidePixel) {
    Pixel p(0,0,255);
    Point a(0,0);
    Point b(0.5,1);
    Segment e(&a, &b);
    ASSERT_DOUBLE_EQ(pow(0.25 * 0.25 + 0.25, 0.5), p.intersectionLength(e));
}

// both endpoints inside pixel
TEST(LengthTest, Contained) {
    Pixel p(0,0,255);
    Point a(-0.25,0);
    Point b(0.25,0);
    Segment e(&a, &b);
    ASSERT_DOUBLE_EQ(0.5, p.intersectionLength(e));
}

TEST(ContainsTest, Inside) {
    Pixel p(0,0,255);
    Point a(0.2,-0.3);
    ASSERT_TRUE(p.containsPoint(a));
}

TEST(ContainsTest, Boundary) {
    Pixel p(0,0,255);
    Point a(0.5,-0.3);
    ASSERT_TRUE(p.containsPoint(a));
}

TEST(ContainsTest, Corner) {
    Pixel p(0,0,255);
    Point a(0.5,0.5);
    ASSERT_TRUE(p.containsPoint(a));
}

TEST(ContainsTest, Outside) {
    Pixel p(0,0,255);
    Point a(1,1);
    ASSERT_FALSE(p.containsPoint(a));
}

// pixel entirely contained within triangle
TEST(AreaTest, PixelContained) {
    Pixel p(1,1,255);
    Point a(0,0);
    Point b(3,0);
    Point c(0,4);
    Triangle t(&a, &b, &c);
    ASSERT_EQ(1,p.intersectionArea(t));
}

// one vertex of triangle inside pixel
TEST(AreaTest, VertexInside) {
    Pixel p(0,0,255);
    Point a(0,0);
    Point b(3,0);
    Point c(0,4);
    Triangle t(&a, &b, &c);
    ASSERT_EQ(0.25, p.intersectionArea(t));
}

// two vertices of triangle inside pixel
TEST(AreaTest, EdgeInside) {
    Pixel p(1,1,255);
    Point a(0.75, 1);
    Point b(1.25, 1);
    Point c(1, 2);
    Triangle t(&a, &b, &c);
    ASSERT_EQ(0.25 * 0.75, p.intersectionArea(t));
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}