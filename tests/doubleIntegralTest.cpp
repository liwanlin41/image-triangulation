#include <gtest/gtest.h>
#include "../src/doubleIntegral.hpp"
#include "../src/triangle.hpp"
#include "../src/point.hpp"

using namespace std;

// for testing purposes, create an 1000 x 500 image all of whose pixels have color 1
vector<vector<Pixel>> generateFakeImage() {
    vector<vector<Pixel>> image;
    for(int i = 0; i < 1000; i++) {
        vector<Pixel> holder;
        for(int j = 0; j < 500; j++) {
            holder.push_back(Pixel(i, j, 1));
        }
        image.push_back(holder);
    }
    return image;
}

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
    vector<vector<Pixel>> image = generateFakeImage();
    Point a(0,0);
    Point b(10,10);
    Point c(20,30);
    Triangle triangle(&a, &b, &c);
    ASSERT_EQ(0, DoubleIntegral::evaluate(z, &image, &triangle));
    ASSERT_EQ(0, DoubleIntegral::evaluate(z, &image, &a, &b, &c));
}

TEST(IntegralTest, Constant) {
    vector<vector<Pixel>> image = generateFakeImage();
    Point a(0,0);
    Point b(3,0);
    Point c(0,4);
    Triangle triangle(&a, &b, &c);
    ASSERT_EQ(6, DoubleIntegral::evaluate(constant, &image, &triangle));
    ASSERT_EQ(6, DoubleIntegral::evaluate(constant, &image, &a, &b, &c));
}

TEST(IntegralTest, Linear) {
    vector<vector<Pixel>> image = generateFakeImage();
    Point a(0,0);
    Point b(3,0);
    Point c(0,4);
    Triangle triangle(&a, &b, &c);
    ASSERT_EQ(14, DoubleIntegral::evaluate(linear, &image, &triangle));
    ASSERT_EQ(14, DoubleIntegral::evaluate(linear, &image, &a, &b, &c));
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}