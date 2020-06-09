//#include <iostream>
//#include <fstream>
#define cimg_use_png 1
#define cimg_use_jpg 1
#include "CImg.h"

/*
#include "point.hpp"
#include "segment.hpp"
#include "matrix.hpp"
#include "triangle.hpp"
*/

using namespace cimg_library;

// eventually input will be the path to an image file?
int main(int argc, char* argv[]) {
	// hard code for now
	CImg<unsigned char> image("../images/piecewise_constant.png");
	image.display("Image");
	return 0;
/*
	CImg<unsigned char> img(640,400,1,3);
	img.fill(0);
	unsigned char purple[] = { 255, 0, 255 };
	img.draw_text(100,100,"Hello", purple);
	img.display("CImg");
	return 0;
*/
}
