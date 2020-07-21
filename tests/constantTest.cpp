#include <gtest/gtest.h>
#include "../src/constant.h"

using namespace std;

// generate 100 x 50 bw image with black rectangle in upper left corner
CImg<unsigned char> generateFakeImage() {
    CImg<unsigned char> result(100, 50, 1, 1, 0);
    unsigned char white[] = {255, 255, 255};
    unsigned char black[] = {0,0,0};
    result.draw_rectangle(0,0,100,50, white, 1);
    result.draw_rectangle(0, 0, 50, 25, black, 1);
    return result;
}

// generate 100 x 50 bw image with a black rectangle in the upper left corner
/*
vector<vector<Pixel>> generateFakeImage() {
    vector<vector<Pixel>> image;
    for(int i = 0; i < 100; i++) {
        vector<Pixel> column;
        for(int j = 0; j < 50; j++) {
            int color = (i < 50 && j < 25) ? 0 : 1;
            Pixel p(i, j, color);
            column.push_back(p);
        }
        image.push_back(column);
    }
    return image;
}
*/

//vector<vector<Pixel>> image = generateFakeImage();
CImg<unsigned char> image = generateFakeImage();
ConstantApprox approx(&image, 0.5);

TEST(EnergyTest, Initial) {
    // top and left triangles are the only two
    // with positive energy contribution; difference is 0.5
    double expectedEnergy = 0.25 * (25 * 50 + 100 * 25) / 2;
    ASSERT_EQ(expectedEnergy, approx.computeEnergy());
}

// check that all the assertions hold
TEST(GradientTest, Syntax) {
    approx.computeGrad();
    approx.gradUpdate();
}

TEST(RunTest, Syntax) {
    approx.run();
}