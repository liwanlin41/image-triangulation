#include <iostream>
#include <fstream>
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
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

// let polyscope read values from point
double adaptorF_custom_accessVector2Value(const Point& p, unsigned int ind) {
    // reflect everything so that it displays correctly
    if (ind == 0) return -p.getX(); 
    if (ind == 1) return -p.getY();
    throw std::logic_error("bad access");
    return -1.;
}

/*
int maxIter = 1000;
double eps = 0.001;
double prevEnergy = 100 * eps; // placeholder values
double newEnergy = 0;
int iterCount = 0;
*/

/*
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
*/

// initialize needed objects
/*
vector<vector<Pixel>> image = generateFakeImage();
ConstantApprox approx(&image, 4, 0.5);
*/

/*
void initialize() {
    auto triangulation = polyscope::registerSurfaceMesh2D("Triangulation", approx.getVertices(), approx.getFaces());
    auto colors = triangulation->addFaceColorQuantity("approximate colors", approx.getColors());
    // allow colors by default
    colors->setEnabled(true);
    // setup gradient descent
    newEnergy = approx.computeEnergy();
    // initialize to something higher than newEnergy
    prevEnergy = newEnergy + 100 * eps;
    iterCount = 0;
}
*/

/*
// run the triangulation process
void step() {
    if(iterCount < maxIter && prevEnergy - newEnergy > eps) {
        cout << "iteration " << iterCount << endl;
        approx.computeGrad();
        while(!approx.gradUpdate()) {
            approx.undo(); // keep halving stepSize until it works
        }
        approx.updateApprox();
        prevEnergy = newEnergy;
        newEnergy = approx.computeEnergy();
        cout << "new energy: " << newEnergy << endl;
        cout << "Step size: " << approx.getStep() << endl;
        iterCount++;
        auto triangulation = polyscope::getSurfaceMesh("Triangulation");
        triangulation->updateVertexPositions2D(approx.getVertices());
        triangulation->addFaceColorQuantity("approximate colors", approx.getColors());
    } else {
        polyscope::warning("done");
    }
}
*/

/*
// callback function to play mesh updates on GUI
void callback() {
    ImGui::InputInt("max iterations", &maxIter); 
    ImGui::InputDouble("stopping condition", &eps);  

    if (ImGui::Button("Start")) {
        initialize();
    }
    // run the triangulation
    if (ImGui::Button("Step")) {
        step();
    }
}
*/

// eventually input will be the path to an image file?
int main(int argc, char* argv[]) {
    int maxIter = 10;
    double eps = 0.001;
    double prevEnergy = 100 * eps; // placeholder values
    double newEnergy = 0;
    int iterCount = 0;

    // take image input
    vector<vector<Pixel>> pixVec; // hold image pixels
    // read image pixels
    CImg<unsigned char> image("../images/black_white.png");
	for(int x = 0; x < image.width(); x++) {
        vector<Pixel> imCol; // column of image, but will be stored as row vector
		for(int y = 0; y < image.height(); y++) {
            Pixel pix(x, y, getLuminance(&image, x, y));
            imCol.push_back(pix);
		}
        pixVec.push_back(imCol);
	}
    ConstantApprox approx(&pixVec, 4, 0.5);

    // lambda for initializing triangulation
    auto initialize = [&approx, &newEnergy, &prevEnergy, &eps, &iterCount]() {
        auto triangulation = polyscope::registerSurfaceMesh2D("Triangulation", approx.getVertices(), approx.getFaces());
        auto colors = triangulation->addFaceColorQuantity("approximate colors", approx.getColors());
        // allow colors by default
        colors->setEnabled(true);
        // setup gradient descent
        newEnergy = approx.computeEnergy();
        // initialize to something higher than newEnergy
        prevEnergy = newEnergy + 100 * eps;
        iterCount = 0;
    };

    // lambda for running a step
    auto step = [&]() {
        if(iterCount < maxIter && prevEnergy - newEnergy > eps) {
            cout << "iteration " << iterCount << endl;
            approx.computeGrad();
            while(!approx.gradUpdate()) {
                approx.undo(); // keep halving stepSize until it works
            }
            approx.updateApprox();
            prevEnergy = newEnergy;
            newEnergy = approx.computeEnergy();
            cout << "new energy: " << newEnergy << endl;
            cout << "Step size: " << approx.getStep() << endl;
            iterCount++;
            auto triangulation = polyscope::getSurfaceMesh("Triangulation");
            triangulation->updateVertexPositions2D(approx.getVertices());
            triangulation->addFaceColorQuantity("approximate colors", approx.getColors());
        } else {
            polyscope::warning("done");
        }
    };

    // lambda to handle GUI updates
    auto callback = [&]() {
        ImGui::InputInt("max iterations", &maxIter); 
        ImGui::InputDouble("stopping condition", &eps);  

        if (ImGui::Button("Start")) {
            initialize();
        }
        // run the triangulation
        if (ImGui::Button("Step")) {
            step();
        }
    };

    polyscope::init();
    polyscope::view::style = polyscope::view::NavigateStyle::Planar;
    polyscope::state::userCallback = callback;
    polyscope::show();
	return 0;
}
