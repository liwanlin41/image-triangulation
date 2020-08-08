#include <iostream>
#include <fstream>
#define cimg_use_png 1
#define cimg_use_jpg 1
#include "CImg.h"

using namespace std;
using namespace cimg_library;

int main(int argc, char* argv[]) {
    // default image path and density
    const char *imgPath;

    string inputPath; // to ensure non-null pointer later; find image directory
    if (argc >= 2) {
        inputPath = argv[1];
        imgPath = inputPath.c_str();
    } else {
        cout << "Please choose an image to crop." << endl;
    }

    CImg<unsigned char> image(imgPath);
    CImg<unsigned char> cropped = image.autocrop();
    cropped.save_jpeg("cropped.jpeg");
    return 0;
}