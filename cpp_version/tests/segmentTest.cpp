#include <iostream>
#include <assert.h>
#include <gtest/gtest.h>
#include "../headers/segment.hpp"

using namespace std;

TEST(LengthTest, Positive) {
    Point a(0,0);
    Point b(3,4);
    Segment e(a,b);
    ASSERT_EQ(5, e.length());
}

TEST(LengthTest, Zero) {
    Point a(0,0);
    Segment e(a,a);
    ASSERT_EQ(0, e.length());
}

TEST(IntersectTest, NormalIntersect) {
    Point a(0,0);
    Point b(3,4);
    Point c(0,-2);
    Point d(3,8);
    Segment e(a,b);
    Segment f(c,d);
    ASSERT_TRUE(e.intersects(f));
    Point intersection = e.getIntersection(f);
    Point expected(1,4.0/3);
    ASSERT_EQ(expected, intersection);
}

TEST(IntersectTest, Parallel) {
    Point a(0,0);
    Point b(3,4);
    Point shiftA(0,2);
    Point shiftB(3,6);
    Segment e(a,b);
    Segment parallel(shiftA, shiftB);
    ASSERT_FALSE(e.intersects(parallel));
}

// test collinear segments with positive overlap length
TEST(IntersectTest, Overlap) {
    Point a(0,0);
    Point b(2,0);
    Point c(1,0);
    Point d(3,0);
    Segment original(a,b);
    Segment overlap(c,d);
    ASSERT_TRUE(original.intersects(overlap));
}

// test no intersection 
TEST(IntersectTest, NoIntersect) {
    Point a(0,0);
    Point b(2,1);
    Point c(3,3);
    Point d(5,0);
    Segment e(a,b);
    Segment f(c,d);
    ASSERT_FALSE(e.intersects(f));
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}