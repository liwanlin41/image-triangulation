#include <gtest/gtest.h>
#include "../src/constant.hpp"

using namespace std;

// generate 100 x 50 bw image with a black rectangle in the upper left corner
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

vector<vector<Pixel>> image = generateFakeImage();
ConstantApprox approx(&image, 3, 0.5);

TEST(EnergyTest, Initial) {
    // top and left triangles are the only two
    // with positive energy contribution; difference is 0.5
    double expectedEnergy = 0.25 * (25 * 50 + 100 * 25) / 2;
    ASSERT_EQ(expectedEnergy, approx.computeEnergy());
}

TEST(GradientTest, Syntax) {
    approx.computeGrad();
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}