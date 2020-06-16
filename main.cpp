#include <iostream>
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "constant.hpp"


const double RED_LUMINANCE = 0.2126;
const double GREEN_LUMINANCE = 0.7152;
const double BLUE_LUMINANCE = 0.0722;

// get the luminance at pixel (x, y) by standard luminance transformation
/*
int getLuminance(CImg<unsigned char> *img, int x, int y) {
	int red_val = (int) (*img)(x, y, 0, 0);	
	int green_val = (int) (*img)(x, y, 0, 1);
	int blue_val = (int) (*img)(x, y, 0, 2);
	return round(red_val * RED_LUMINANCE + green_val * GREEN_LUMINANCE + blue_val * BLUE_LUMINANCE);
}
*/

// let polyscope read values from point
double adaptorF_custom_accessVector2Value(const Point& p, unsigned int ind) {
    // reflect everything so that it displays correctly
    if (ind == 0) return -p.getX(); 
    if (ind == 1) return -p.getY();
    throw std::logic_error("bad access");
    return -1.;
}

const int maxIter = 1000;
const double eps = 0.001;

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
    polyscope::init();
    polyscope::view::style = polyscope::view::NavigateStyle::Planar;

	vector<vector<Pixel>> image = generateFakeImage();
	ConstantApprox approx(&image, 4, 0.5);
    auto triangulation = polyscope::registerSurfaceMesh2D("Triangulation", approx.getVertices(), approx.getFaces());
    //auto colors = polyscope::getSurfaceMesh("Triangulation")->addFaceColorQuantity("approxiate colors", approx.getColors());
    // auto colors = triangulation->addFaceColorQuantity("approximate colors", approx.getColors());
    // allow colors by default
    // colors->setEnabled(true);
    // polyscope::show();

    /*
    // track change in energy for stopping point
    double newEnergy = approx.computeEnergy();
    // initialize to something higher than newEnergy
    double prevEnergy = newEnergy + 100 * eps;
    int iterCount = 0;
    while(iterCount < maxIter && abs(prevEnergy-newEnergy) > eps) {
        cout << "iteration " << iterCount << endl;
        approx.computeGrad();
        while(!approx.gradUpdate()) {
            approx.undo(); // keep halving stepSize until it works
        }
        approx.updateApprox();
        prevEnergy = newEnergy;
        newEnergy = approx.computeEnergy();
        cout << "new energy: " << newEnergy << endl;
        cout << "step size: " << approx.getStep() << endl;
        iterCount++;
        polyscope::getSurfaceMesh("Triangulation")->updateVertexPositions2D(approx.getVertices());
        polyscope::getSurfaceMesh("Triangulation")->addFaceColorQuantity("approxiate colors", approx.getColors());
        polyscope::draw();
    }
    */
    approx.run();
    triangulation->updateVertexPositions2D(approx.getVertices());
    auto colors = triangulation->addFaceColorQuantity("approximate colors", approx.getColors());
    colors->setEnabled(true);
    // polyscope::getSurfaceMesh("Triangulation")->addFaceColorQuantity("approximate colors", approx.getColors());
    polyscope::show();
	return 0;
}
