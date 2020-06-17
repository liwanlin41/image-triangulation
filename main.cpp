#include <iostream>
#include <fstream>
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#define cimg_use_png 1
#define cimg_use_jpg 1
#include "CImg.h"
#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"

#include "constant.hpp"

using namespace cimg_library;
using namespace matlab::engine;


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

// eventually input will be the path to an image file?
int main(int argc, char* argv[]) {
    const char *imgPath = "../images/clouds.jpg";
    double density = 0.005;
    if (argc >= 2) {
        density = atof(argv[1]);
    }
    cout << density << endl;
    /*
    // string imgPath;
    if (argc >= 2) {
        imgPath = argv[1];
    } else {
        imgPath = "../images/tim.png"; // default image path
    }
    */
    // set default values
    int maxIter = 20;
    double eps = 0.001;
    double prevEnergy = 100 * eps; // placeholder values
    double newEnergy = 0;
    int iterCount = 0;
    // double density = 0.005; // TODO: figure out how to set this        

    // take image input
    vector<vector<Pixel>> pixVec; // hold image pixels
    // read image pixels
    CImg<unsigned char> image(imgPath);
    bool isGrayscale = (image.spectrum() == 1);
    cout << image.width() << " x " << image.height() << endl;
	for(int x = 0; x < image.width(); x++) {
        vector<Pixel> imCol; // column of image, but will be stored as row vector
		for(int y = 0; y < image.height(); y++) {
            int color;
            if(isGrayscale) {
                color = image(x, y);
            } else {
                color = getLuminance(&image, x, y);
            }
            Pixel pix(x, y, color);
            imCol.push_back(pix);
		}
        pixVec.push_back(imCol);
	}

    // start matlab to get initial triangulation
    std::unique_ptr<MATLABEngine> matlabPtr = startMATLAB();
    matlab::data::ArrayFactory factory;
    // add path to TRIM code
    vector<matlab::data::Array> genPathArgs({
        factory.createCharArray("../deps/trim")
    });
    auto generatedPath = matlabPtr->feval(u"genpath",genPathArgs);
    matlabPtr->feval(u"addpath", generatedPath);

    // read image
    vector<matlab::data::Array> pathToImage({
        factory.createCharArray(imgPath)
    });
    auto img = matlabPtr->feval(u"imread", pathToImage);

    // generate triangulation
    vector<matlab::data::Array> tArgs({
        img,
        factory.createScalar<double>(density) // density 
    });
    cout << "Triangulating...\n";
    vector<matlab::data::Array> output = matlabPtr->feval(u"imtriangulate", 3, tArgs);
    cout << "done\n";
    // vertices of triangulation
    matlab::data::Array vertices = output.at(0);
    matlab::data::Array triangleConnections = output.at(1);
    int n = vertices.getDimensions().at(0); // number of points in triangulation

    // initialize points of mesh
    cout << "Getting " << n << " points...\n";
    vector<Point> points;
    for(int i = 0; i < n; i++) {
        // note these are 1-indexed pixel values; will need
        // to convert to usable points 
        int x = vertices[i][0];
        int y = vertices[i][1];
        // determine whether this point lies on edge of image
        bool isBorderX = (x == 1 || x == pixVec.size());
        bool isBorderY = (y == 1 || y == pixVec.at(0).size());
        Point p(x-1.5, y-1.5, isBorderX, isBorderY); // translate to coordinate system in this code
        points.push_back(p);
    }
    cout << "done\n";

    // convert connections to vector<vector<int>>
    cout << "Getting edges...\n";
    vector<array<int, 3>> edges;
    int f = triangleConnections.getDimensions().at(0); // number of triangles
    cout << "number of triangles: " << f << endl;
    for(int i = 0; i < f; i++) {
        array<int, 3> vertexInds;
        for(int j = 0; j < 3; j++) {
            // matlab is 1 indexed for some bizarre reason;
            // change back to zero indexing
            int ind = triangleConnections[i][j];
            vertexInds[j] = ind - 1;
        }
        edges.push_back(vertexInds);
    }
    cout << "done\n";

    cout << "Initializing mesh...\n";
    ConstantApprox approx(&pixVec, &points, edges, 0.5);
    cout << "ready\n";

    // lambda for initializing triangulation
    auto initialize = [&approx, &newEnergy, &prevEnergy, &eps, &iterCount]() {
        cout << "creating mesh\n";
        auto triangulation = polyscope::registerSurfaceMesh2D("Triangulation", approx.getVertices(), approx.getFaces());
        cout << "getting colors\n";
        auto colors = triangulation->addFaceColorQuantity("approximate colors", approx.getColors());
        // allow colors by default
        colors->setEnabled(true);
        // setup gradient descent
        cout << "finding energy..." << endl;
        newEnergy = approx.computeEnergy();
        cout << "done\n";
        // initialize to something higher than newEnergy
        prevEnergy = newEnergy + 100 * eps;
        iterCount = 0;
    };

    // lambda for running a step
    auto step = [&]() {
        if(iterCount < maxIter && abs(prevEnergy - newEnergy) > eps) {
            cout << "iteration " << iterCount << endl;
            approx.computeGrad();
            while(!approx.gradUpdate()) {
                approx.undo(); // keep halving stepSize until it works
            }
            approx.updateApprox();
            prevEnergy = newEnergy;
            newEnergy = approx.computeEnergy();
            if(newEnergy > prevEnergy) { // overshot optimum?
                do {
                    approx.undo();
                } while (!approx.gradUpdate());
                approx.updateApprox();
                newEnergy = approx.computeEnergy();
            }
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
