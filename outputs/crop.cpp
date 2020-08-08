#include <iostream>
#include <fstream>
#define cimg_use_png 1
#define cimg_use_jpg 1
#include "CImg.h"

using namespace std;
using namespace cimg_library;

int main(int argc, char* argv[]) {
    const char *imgPath;
    const char *outputName;

    string inputPath; 
    string outputPath = "cropped.jpeg"; // default output file name
    if(argc < 2) {
        cout << "Please choose an image to crop." << endl;
        return 0;
    }
    inputPath = argv[1];
    imgPath = inputPath.c_str();
    if (argc == 3) {
        outputPath = argv[2];
    }
    outputName = outputPath.c_str();

    CImg<unsigned char> image(imgPath);
    CImg<unsigned char> cropped = image.autocrop();
    cropped.save_jpeg(outputName);
    return 0;
}