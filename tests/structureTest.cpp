// test objects: point, segment, triangle, pixel

#include <gtest/gtest.h>
#include "triangle.cuh"

using namespace std;

// note that device side functions cannot be directly called
// and have all been commented out;
// they are indirectly tested in integration calls

// tests for Point

TEST(PointTest, Movement) {
    Point a(0,0);
    a.move(3, 4.5);
    //ASSERT_EQ(a, a);
    ASSERT_EQ(3, a.getX());
    ASSERT_EQ(4.5, a.getY());
}

/*
TEST(PointTest, Equality) {
    Point a(0,0);
    Point b(0,0);
    ASSERT_EQ(a, b);
    a.move(1,1);
    ASSERT_NE(a, b);
    b.move(1,1);
    ASSERT_EQ(a, b);
}
*/

TEST(PointTest, Distance) {
    Point a(0,0);
    Point b(3,4);
    ASSERT_EQ(5, a.distance(b));
}

// tests for segment

TEST(SegmentTest, PositiveLength) {
    Point a(0,0);
    Point b(3,4);
    Segment e(&a,&b);
    ASSERT_EQ(5, e.length());
}

TEST(SegmentTest, ZeroLength) {
    Point a(0,0);
    Segment e(&a,&a);
    ASSERT_EQ(0, e.length());
}

// test segment intersection; note these tests
// cannot be directly called because they are device functions,
// but they were thoroughly tested for accuracy before moving to GPU
// instead, they are tested indirectly with integration

/*
TEST(SegmentTest, IntersectionExists) {
    Point a(0,0);
    Point b(3,4);
    Point c(0,-2);
    Point d(3,8);
    Segment e(&a,&b);
    Segment f(&c,&d);
    Point intersection;
    ASSERT_TRUE(e.intersection(f, &intersection));
    Point expected(1,4.0/3);
    ASSERT_EQ(expected, intersection);
}

TEST(SegmentTest, IntersectionParallel) {
    Point a(0,0);
    Point b(3,4);
    Point shiftA(0,2);
    Point shiftB(3,6);
    Segment e(&a,&b);
    Segment parallel(&shiftA, &shiftB);
    ASSERT_FALSE(e.intersection(parallel));
}

TEST(SegmentTest, IntersectionOverlap) {
    Point a(0,0);
    Point b(2,0);
    Point c(1,0);
    Point d(3,0);
    Segment original(&a,&b);
    Segment overlap(&c,&d);
    ASSERT_FALSE(original.intersection(overlap));
}

// test collinear segments where one segment is entirely contained within the other
TEST(SegmentTest, IntersectionContainment) {
    Point a(0,0);
    Point b(3,0);
    Point c(1,0);
    Point d(2,0);
    Segment e(&a,&b);
    Segment f(&c,&d);
    ASSERT_FALSE(e.intersection(f));
    ASSERT_FALSE(f.intersection(e));
}

// test collinear segments with no overlap
TEST(SegmentTest, IntersectDisjointCollinear) {
    Point a(0,0);
    Point b(1,0);
    Point c(2,0);
    Point d(3,0);
    Segment e(&a,&b);
    Segment f(&c,&d);
    ASSERT_FALSE(e.intersection(f));
}

// test no intersection 
TEST(SegmentTest, IntersectNone) {
    Point a(0,0);
    Point b(2,1);
    Point c(3,3);
    Point d(5,0);
    Segment e(&a,&b);
    Segment f(&c,&d);
    ASSERT_FALSE(e.intersection(f));
}
*/

// test unit normal to vertical segment
TEST(SegmentTest, HorizontalNormal) {
    Point a(0,0);
    Point b(0,1);
    Segment e(&a, &b);
    double nx, ny;
    e.scaledNormal(&nx, &ny);
    ASSERT_EQ(0.5, nx);
    ASSERT_EQ(0, ny);
}

// test orientation of a slanted segment
TEST(SegmentTest, Slant) {
    Point a(6,8);
    Point b(0,0);
    Segment e(&a, &b);
    double nx, ny;
    e.scaledNormal(&nx, &ny);
    ASSERT_EQ(-4, nx);
    ASSERT_EQ(3, ny);
}

// tests for Triangle

// test area of degenerate triangle
TEST(TriangleTest, DegenerateArea) {
    Point a(0,0);
    Point b(0,1);
    Point c(0,2);
    Triangle triangle(&a, &b, &c);
    ASSERT_EQ(0, triangle.getArea());
}

TEST(TriangleTest, RightTriangleArea) {
    Point a(0,0);
    Point b(0,3);
    Point c(4,0);
    Triangle triangle(&a, &b, &c);
    ASSERT_EQ(6, triangle.getArea());
}

// test that constructor re-orders in ccw direction
TEST(TriangleTest, Clockwise) {
    Point a(0,0);
    Point b(2,3);
    Point c(5,0);
    Triangle triangle(&a, &b, &c);
    ASSERT_EQ(7.5, triangle.getSignedArea());
}

TEST(TriangleTest, CCWSignedArea) {
    Point a(0,0);
    Point b(5,0);
    Point c(3,4);
    Triangle triangle(&a, &b, &c);
    ASSERT_EQ(10, triangle.getSignedArea());
}

// test change in area when a triangle vertex moves
TEST(TriangleTest, VertexMovement) {
    Point a(0,0);
    Point b(5,0);
    Point c(3,4);
    Triangle triangle(&a, &b, &c);
    c.move(0,-8);
    ASSERT_EQ(-10, triangle.getSignedArea());
    ASSERT_EQ(10, triangle.getArea());
}

// vertex moves parallel to opposite edge
TEST(TriangleTest, DerivativeParallel) {
    Point a(0,0);
    Point b(5,0);
    Point c(3,4);
    Triangle triangle(&a, &b, &c);
    int mvPt = 2;
    double change = triangle.dA(mvPt, 5, 0);
    ASSERT_EQ(0, change);
}

// vertex moves perpendicular to opposite dge
TEST(TriangleTest, DerivativePerpendicular) {
    Point a(5,0);
    Point b(3,4);
    Point c(0,0);
    Triangle triangle(&a, &b, &c);
    int mvPt = 1;
    double change = triangle.dA(mvPt, 0, 2);
    ASSERT_EQ(5, change);
}

/*
TEST(TriangleTest, ContainedInside) {
    Point a(0,0);
    Point b(3,0);
    Point c(0,4);
    Point p(1,1);
    Triangle triangle(&a, &b, &c);
    ASSERT_TRUE(triangle.contains(p));
}

TEST(TriangleTest, ContainedBoundary) {
    Point a(0,0);
    Point b(3,0);
    Point c(0,4);
    Point p(2,0);
    Triangle triangle(&a, &b, &c);
    ASSERT_TRUE(triangle.contains(p));
}

// test point outside triangle
TEST(TriangleTest, ContainedOutside) {
    Point a(0,0);
    Point b(3,0);
    Point c(0,4);
    Point p(-0.1,5);
    Triangle triangle(&a, &b, &c);
    ASSERT_FALSE(triangle.contains(p));
}
*/