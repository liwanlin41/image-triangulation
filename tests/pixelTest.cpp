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

TEST(COntainsTest, Outside) {
    Pixel p(0,0,255);
    Point a(1,1);
    ASSERT_FALSE(p.containsPoint(a));
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}