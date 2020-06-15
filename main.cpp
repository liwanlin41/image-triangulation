#include <iostream>
//#include <fstream>
#define cimg_use_png 1
#define cimg_use_jpg 1
#include "CImg.h"

#include "constant.hpp"

using namespace cimg_library;

const double RED_LUMINANCE = 0.2126;
const double GREEN_LUMINANCE = 0.7152;
const double BLUE_LUMINANCE = 0.0722;

// get the luminance at pixel (x, y) by standard luminance transformation
int getLuminance(CImg<unsigned char> *img, int x, int y) {
	int red_val = (int) (*img)(x, y, 0, 0);	
	int green_val = (int) (*img)(x, y, 0, 1);
	int blue_val = (int) (*img)(x, y, 0, 2);
	return round(red_val * RED_LUMINANCE + green_val * GREEN_LUMINANCE + blue_val * BLUE_LUMINANCE);
}

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

// eventually input will be the path to an image file?
int main(int argc, char* argv[]) {
	/*
	CImg<unsigned char> image("../images/black_white.png");
	for(int x = 0; x < image.width(); x+= 20) {
		for(int y = 0; y < image.height(); y+= 20) {
			std::cout << x << ", " << y << ": " << getLuminance(&image, x, y) << std::endl;
		}
	}
	image.display("Image");
	return 0;
	*/
	vector<vector<Pixel>> image = generateFakeImage();
	ConstantApprox approx(&image, 4, 0.05);
	CImg<unsigned char> before = approx.show();
	before.display("Before");
	// step by step runthrough
	/*
	double newEnergy = approx.computeEnergy();
    // initialize to something higher than newEnergy
	int maxIter = 1000;
	double eps = 0.0001;
    double prevEnergy = newEnergy + 100 * eps;
    int iterCount = 0;
    while(iterCount < maxIter && prevEnergy-newEnergy > eps) {
        cout << "iteration " << iterCount << endl;
        approx.computeGrad();
        while(!approx.gradUpdate()) {
            approx.undo(); // keep halving stepSize until it works
        }
        approx.updateApprox();
        prevEnergy = newEnergy;
        newEnergy = approx.computeEnergy();
        while(newEnergy > prevEnergy) { // overshot optimum?
            do {
                approx.undo();
            } while (!approx.gradUpdate());
            approx.updateApprox();
            newEnergy = approx.computeEnergy();
        }
        cout << "new energy: " << newEnergy << endl;
        iterCount++;
		CImg<unsigned char> result = approx.show();
		result.display("Step");
    }
	*/
	approx.run();
	CImg<unsigned char> result = approx.show();
	result.display("Result");
	return 0;
}
