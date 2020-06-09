#include <iostream>
#include <assert.h>
#include <gtest/gtest.h>
#include "../src/segment.hpp"

using namespace std;

TEST(LengthTest, Positive) {
    Point a(0,0);
    Point b(3,4);
    Segment e(&a,&b);
    ASSERT_EQ(5, e.length());
}

TEST(LengthTest, Zero) {
    Point a(0,0);
    Segment e(&a,&a);
    ASSERT_EQ(0, e.length());
}

TEST(IntersectTest, NormalIntersect) {
    Point a(0,0);
    Point b(3,4);
    Point c(0,-2);
    Point d(3,8);
    Segment e(&a,&b);
    Segment f(&c,&d);
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
    Segment e(&a,&b);
    Segment parallel(&shiftA, &shiftB);
    ASSERT_FALSE(e.intersects(parallel));
}

// test collinear segments with positive overlap length
TEST(IntersectTest, Overlap) {
    Point a(0,0);
    Point b(2,0);
    Point c(1,0);
    Point d(3,0);
    Segment original(&a,&b);
    Segment overlap(&c,&d);
    ASSERT_TRUE(original.intersects(overlap));
}

// test collinear segments where one segment is entirely contained within the other
TEST(IntersectTest, Containment) {
    Point a(0,0);
    Point b(3,0);
    Point c(1,0);
    Point d(2,0);
    Segment e(&a,&b);
    Segment f(&c,&d);
    ASSERT_TRUE(e.intersects(f));
    ASSERT_TRUE(f.intersects(e));
}

// test collinear segments with no overlap
TEST(IntersectTest, DisjointCollinear) {
    Point a(0,0);
    Point b(1,0);
    Point c(2,0);
    Point d(3,0);
    Segment e(&a,&b);
    Segment f(&c,&d);
    ASSERT_FALSE(e.intersects(f));
}

// test no intersection 
TEST(IntersectTest, NoIntersect) {
    Point a(0,0);
    Point b(2,1);
    Point c(3,3);
    Point d(5,0);
    Segment e(&a,&b);
    Segment f(&c,&d);
    ASSERT_FALSE(e.intersects(f));
}

TEST(NormalTest, VerticalSegment) {
    Point a(0,0);
    Point b(0,1);
    Segment e(&a, &b);
    Matrix unit = e.unitNormal();
    Matrix expected(1,0);
    ASSERT_EQ(expected, unit);
}

// test orientation of a slanted segment
TEST(NormalTest, Slant) {
    Point a(6,8);
    Point b(0,0);
    Segment e(&a, &b);
    Matrix unit = e.unitNormal();
    Matrix scaled = e.scaledNormal();
    Matrix expectedUnit(-0.8,0.6);
    Matrix expectedScale(-4, 3);
    ASSERT_EQ(expectedUnit, unit);
    ASSERT_EQ(expectedScale, scaled);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
