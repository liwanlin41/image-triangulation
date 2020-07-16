#include <iostream>
#include <fstream>
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"
#include "polyscope/view.h"
#define cimg_use_png 1
#define cimg_use_jpg 1
#include "CImg.h"
#include <cuda.h>
#include <cuda_runtime.h>

#include "src/simulator.h"

using namespace cimg_library;
using namespace matlab::engine;

static const double DENSITY_DEFAULT = 0.05; // experimentally a decent starting density
static const int DX_DEFAULT = 50; // sampling default

// let polyscope read values from point
double adaptorF_custom_accessVector2Value(const Point& p, unsigned int ind) {
    // reflect everything so that it displays correctly
    if (ind == 0) return -p.getX(); 
    if (ind == 1) return -p.getY();
    throw std::logic_error("bad access");
    return -1.;
}

// to aid in mesh registration, put constant triangle indexing here
const array<int, 3> triInds = {0,1,2};
const vector<array<int, 3>> singleTriangle = {triInds};

void registerMesh(Approx *approx, bool first) {
    if(approx->getApproxType() == constant) {
        auto triangulation = polyscope::registerSurfaceMesh2D("Triangulation", approx->getVertices(), approx->getFaces());
        auto colors = triangulation->addFaceColorQuantity("approximate colors", approx->getColors());
	    if(first) {
		    // allow colors by default
    	    colors->setEnabled(true);
    	    // set material to flat to get more accurate rgb values
    	    triangulation->setMaterial("flat");
	    }
    } else if (approx->getApproxType() == linear) {
        //vector<Point> pts = ((LinearApprox*) approx)->testPoints();
        vector<Point> pts = approx->getVertices();
        vector<array<int, 3>> faces = approx->getFaces();
        //vector<array<int, 3>> faces = ((LinearApprox*) approx)->testFaces();
        vector<array<double, 3>> colors = approx->getColors();
        // register a separate mesh for each triangle
        for(int t = 0; t < faces.size(); t++) {
            vector<Point> thisTriangle; // hold vertices of faces.at(t)
            vector<array<double, 3>> vertexColors;
            for(int i = 0; i < 3; i++) {
                thisTriangle.push_back(pts.at(faces.at(t).at(i)));
                vertexColors.push_back(colors.at(3*t+i));
            }
            auto triangle = polyscope::registerSurfaceMesh2D(to_string(t), thisTriangle, singleTriangle);
            auto colorPtr = triangle->addVertexColorQuantity("linear approx", vertexColors);
            if(first) {
                colorPtr->setEnabled(true);
                triangle->setMaterial("flat");
            }
        }
    }
}

void updateMesh(Approx *approx) {
    if(approx->getApproxType() == constant) {
        auto triangulation = polyscope::getSurfaceMesh("Triangulation");
        triangulation->updateVertexPositions2D(approx->getVertices());
        triangulation->addFaceColorQuantity("approximate colors", approx->getColors());
    } else if(approx->getApproxType() == linear) {
        vector<Point> pts = approx->getVertices();
        vector<array<int, 3>> faces = approx->getFaces();
        vector<array<double, 3>> colors = approx->getColors();
        for(int t = 0; t < faces.size(); t++) {
            vector<Point> thisTriangle;
            vector<array<double, 3>> vertexColors;
            for(int i = 0; i < 3; i++) {
                thisTriangle.push_back(pts.at(faces.at(t).at(i)));
                vertexColors.push_back(colors.at(3*t+i));
            }
            auto triangle = polyscope::getSurfaceMesh(to_string(t));
            triangle->updateVertexPositions2D(thisTriangle);
            triangle->addVertexColorQuantity("linear approx", vertexColors);
        }
    }
}

// eventually input will be the path to an image file?
int main(int argc, char* argv[]) {
    // default image path and density
    const char *imgPath = "../images/flower.jpg";

    string inputPath = "../images/"; // to ensure non-null pointer later; find image directory
    if (argc >= 2) {
        inputPath += argv[1];
        imgPath = inputPath.c_str();
    }

    CImg<unsigned char> image(imgPath);

    // wrapper for running approximation
    Simulator sim(imgPath, &image, linear);

    // set default values
    int maxIter = 100;
    int subdivisions = 50;
    double eps = 0.001;

    // lambda for displaying the starting triangulation
    auto initialize = [&]() {
       sim.initialize();
        // center mesh
        polyscope::view::resetCameraToHomeView();
        polyscope::resetScreenshotIndex();
        // screenshot
        polyscope::screenshot(false);
        polyscope::screenshot("../outputs/initial.tga", false);
    };

    // lambda for making a single step of gradient flow
    auto step = [&]() {
       if(sim.step(maxIter, eps)) {
            polyscope::screenshot(false);
       } else {
           polyscope::warning("done");
       }
    };

    // lambda for running the entire gradient flow
    auto runGradient = [&]() {
       while(sim.step(maxIter, eps)) {
           polyscope::screenshot(false);
       }
    };

    // lambda for retriangulating by subdivision
    auto retriangulate = [&]() {
        sim.retriangulate(subdivisions);
        polyscope::screenshot(false);
    };

    // lambda to handle GUI updates
    vector<string> angryButtonVec = {"MoRe TrIaNgLeS", "MORE TRIANGLES", "Are you done yet?", 
        "MORE", "I SAID MORE", "GIVE ME MORE", "We're done here.", "MORE",
        "Why are you so demanding?", "MORE", "I'm trying to be nice.", "MORE",
        "But you're really trying me here.", "MORE", "Okay, have it your way..", "More Triangles"};
    int numPresses = 0;
    string angryButton = "More Triangles";
    auto callback = [&]() {
        ImGui::InputInt("max iterations", &maxIter); 
        ImGui::InputInt("# edges to divide", &subdivisions);
        ImGui::InputDouble("stopping condition", &eps);  

        if (ImGui::Button("Start")) {
            initialize();
        }

        // step by step
        if (ImGui::Button("Step")) {
            step();
        }
        ImGui::SameLine();
        // run the triangulation
        if (ImGui::Button("Gradient Flow")) {
            runGradient();
        }

        // allow retriangulation
        if (ImGui::Button(angryButton.c_str())) {
            retriangulate();
            if(numPresses >= 3 && numPresses < 3 + angryButtonVec.size()) {
                angryButton = angryButtonVec.at(numPresses - 3);
            } else if (numPresses >= 3 + angryButtonVec.size()) {
                angryButton.append("!");
            }
            numPresses++;
        }
    };

    polyscope::init();
    polyscope::view::style = polyscope::view::NavigateStyle::Planar;
    polyscope::state::userCallback = callback;
    polyscope::show();
    polyscope::screenshot("../outputs/triangulation.tga", false);

    sim.cleanup();
	return 0;
}
