#include <gtest/gtest.h>
#include <math.h>
#include "../src/pixel.hpp"

using namespace std;

vector<vector<Pixel>> generateFakeImage() {
    vector<vector<Pixel>> image;
    for(int i = 0; i < 100; i++) {
        vector<Pixel> column;
        for(int j = 0; j < 50; j++) {
            int color = (i < 50 && j < 25) ? 0 : 255;
            Pixel p(i, j, color);
            column.push_back(p);
        }
        image.push_back(column);
    }
    return image;
}

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

// one vertex inside pixel, intersection resembles a "stripe"
// down the middle
TEST(AreaTest, Bug) {
    Pixel p(10,10,255);
    Point a(10.25, 10.75);
    Point b(8.25, 8.75);
    Point c(10.25, 9.75);
    Triangle t(&a, &b, &c);
    ASSERT_EQ(0.5625, p.intersectionArea(t));
}

// reason for failure: very small parts of the triangle lie outside
// the image and are not counted; in practice hopefully this
// will not happen due to boundary points preventing movement outside
TEST(AreaTest, Bug2) {
    Point a(-0.5, -0.5);
    Point b(58.2550236, -0.5000000000001);
    Point c(40.53696033, 13.74525581);
    vector<vector<Pixel>> img = generateFakeImage();
    Triangle t(&a, &b, &c);
    double area = 0;
    for(int i = 0; i < img.size(); i++) {
        for(int j = 0; j < img.at(0).size(); j++) {
            Pixel p = img.at(i).at(j);
            double A = p.intersectionArea(t);            
            area += A;
        }
    }
    ASSERT_FLOAT_EQ(t.getArea(), area);
}

TEST(AreaTest, Bug3) {
    Point a(99.5, -0.5);
    Point b(79.5, 14.57797241);
    Point c(58.2550236, -0.5);
    vector<vector<Pixel>> img = generateFakeImage();
    Triangle t(&a, &b, &c);
    double area = 0;
    for(int i = 0; i < img.size(); i++) {
        for(int j = 0; j < img.at(0).size(); j++) {
            Pixel p = img.at(i).at(j);
            double A = p.intersectionArea(t);
            area += A;
        }
    }
    ASSERT_FLOAT_EQ(t.getArea(), area);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}