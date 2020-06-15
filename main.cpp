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
  if (ind == 0) return p.getX();
  if (ind == 1) return p.getY();
  throw std::logic_error("bad access");
  return -1.;
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
    polyscope::init();
    polyscope::view::style = polyscope::view::NavigateStyle::Planar;

	vector<vector<Pixel>> image = generateFakeImage();
	ConstantApprox approx(&image, 4, 0.5);
    polyscope::registerSurfaceMesh2D("Triangulation", approx.getVertices(), approx.getFaces());
    polyscope::getSurfaceMesh("Triangulation")->addFaceColorQuantity("approximate colors", approx.getColors());
    polyscope::show();
    approx.run();
    // polyscope::registerSurfaceMesh2D("Triangulation", approx.getVertices(), approx.getFaces());
    // polyscope::getSurfaceMesh("Triangulation")->addFaceColorQuantity("approximate colors", approx.getColors());
    // polyscope::show();
	return 0;
}
