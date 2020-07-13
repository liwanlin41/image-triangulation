#include <iostream>
#include <fstream>
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"
#define cimg_use_png 1
#define cimg_use_jpg 1
#include "CImg.h"
#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"
#include <cuda.h>
#include <cuda_runtime.h>

#include "src/constant.cuh"

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

// eventually input will be the path to an image file?
int main(int argc, char* argv[]) {
    // default image path and density
    const char *imgPath = "../images/flower.jpg";
    double density; // density input for TRIM
    int dx; // take one sample every dx pixels

    string inputPath = "../images/"; // to ensure non-null pointer later; find image directory
    if (argc >= 2) {
        inputPath += argv[1];
        imgPath = inputPath.c_str();
    }
    // set default values
    int maxIter = 100;
    int subdivisions = 5;
    double eps = 0.001;
    double prevEnergy = 100 * eps; // placeholder values
    double newEnergy = 0;
    int iterCount = 0;
    int totalIters = 0; // total iterations over multiple retriangulations

    string trimString;
    bool useTRIM = false; // default to uniform initialization
    cout << "Use TRIM initialization? y/N: ";
    cin >> trimString;
    vector<string> yesAnswers = {"y", "Y"};
    // anything other than y/Y will be false
    for(int i = 0; i < 2; i++) {
        if(trimString == yesAnswers.at(i)) {
            useTRIM = true;
        }
    }

    // create image to read pixels
    CImg<unsigned char> image(imgPath);

    // create approximation instance
    ConstantApprox approx(&image, 0.05);

    // start matlab to get initial triangulation
    std::unique_ptr<MATLABEngine> matlabPtr = startMATLAB();
    matlab::data::ArrayFactory factory;

    // initialize triangulation
    if(useTRIM) {
        // prompt for density
        cout << "Density argument: ";
        cin >> density;
        if(cin.fail()) {
            cin.clear();
            cout << "defaulting to " << DENSITY_DEFAULT << endl;
            density = DENSITY_DEFAULT;
        }
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
        // appears to affect only large images:
        // adjust image boundaries to TRIM result (may crop a line of pixels)
        int minX = 1000;
        int minY = 1000;
        int maxX = 0;
        int maxY = 0;
        for(int i = 0; i < n; i++) {
            int x = vertices[i][0];
            int y = vertices[i][1];
            maxX = max(x, maxX);
            maxY = max(y, maxY);
            minX = min(x, minX);
            minY = min(y, minY);
        }
        for(int i = 0; i < n; i++) {
            // note these are 1-indexed pixel values; will need
            // to convert to usable points 
            int x = vertices[i][0];
            int y = vertices[i][1];
            // determine whether this point lies on edge of image
            bool isBorderX = (x == minX || x == maxX);
            bool isBorderY = (y == minY || y == maxY);
            Point p(x-0.5 - minX, y-0.5 - minY, isBorderX, isBorderY); // translate to coordinate system in this code
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
        approx.initialize(&points, edges);
    } else {
        // prompt for dx
        cout << "Sample once every __ pixels? ";
        cin >> dx;
        if(cin.fail()) {
            cin.clear();
            dx = DX_DEFAULT;
            cout << "defaulting to " << DX_DEFAULT << endl;
        }
        approx.initialize(dx);
    }
    cout << "ready\n";

    vector<double> elapsedTimeVec; // hold cumulative step size
    vector<double> energyVec; // hold energy per iteration
    double totalStep = 0; // running total step

    // lambda for initializing triangulation
    auto initialize = [&]() {
        cout << "creating mesh\n";
        auto triangulation = polyscope::registerSurfaceMesh2D("Triangulation", approx.getVertices(), approx.getFaces());
        cout << "getting colors\n";
        auto colors = triangulation->addFaceColorQuantity("approximate colors", approx.getColors());
        // allow colors by default
        colors->setEnabled(true);
        // set material to flat to get more accurate rgb values
        triangulation->setMaterial("flat");
        // setup gradient descent
        cout << "finding energy..." << endl;
        newEnergy = approx.computeEnergy();
        cout << "done, energy is " << newEnergy << endl;
        // initialize to something higher than newEnergy
        prevEnergy = newEnergy * 2;
        iterCount = 0;
        totalIters = 0;
        elapsedTimeVec.push_back(totalStep); // initial values
        energyVec.push_back(newEnergy);
        polyscope::screenshot(false);
        polyscope::screenshot("../outputs/initial.tga", false);
    };

    // lambda for running a step
    // track number of times for which change in energy has been "small"
    int numSmallChanges = 0;
    const int maxSmallChanges = 3; // maximum number of consecutive "small" changes allowed
    auto step = [&]() {
        if(iterCount < maxIter && numSmallChanges < maxSmallChanges) {
        //if(iterCount < maxIter && abs(prevEnergy - newEnergy) > eps * prevEnergy) {
            cout << "iteration " << iterCount << " (" << totalIters << " total)" << endl;
            // allow energy increases in the first 10 iterations but not after that
            // in order to handle poor initialization
            totalStep += approx.step(prevEnergy, newEnergy, totalIters >= 10);
            // data collection
            elapsedTimeVec.push_back(totalStep);
            energyVec.push_back(newEnergy);
            iterCount++;
            totalIters++;
            if(abs(prevEnergy - newEnergy) > eps * prevEnergy) {
                numSmallChanges = 0;
            } else {
                numSmallChanges++;
            }
            // handle display
            auto triangulation = polyscope::getSurfaceMesh("Triangulation");
            triangulation->updateVertexPositions2D(approx.getVertices());
            triangulation->addFaceColorQuantity("approximate colors", approx.getColors());
            polyscope::screenshot(false);
        } else {
            polyscope::warning("done");
        }
    };

    // lambda for retriangulating by subdivision
    auto retriangulate = [&]() {
        approx.subdivide(subdivisions);
        // re-initialize new mesh
        auto triangulation = polyscope::registerSurfaceMesh2D("Triangulation", approx.getVertices(), approx.getFaces());
        auto colors = triangulation->addFaceColorQuantity("approximate colors", approx.getColors());
        // reset values
        iterCount = 0;
        //totalStep = 0; 
        numSmallChanges = 0;
        totalIters++;
        newEnergy = approx.computeEnergy();
        elapsedTimeVec.push_back(totalStep);
        energyVec.push_back(newEnergy);
        cout << "energy after subdivision: " << newEnergy << endl;
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
        ImGui::InputInt("# triangles on subdivision", &subdivisions);
        ImGui::InputDouble("stopping condition", &eps);  

        if (ImGui::Button("Start")) {
            initialize();
        }
        // run the triangulation
        if (ImGui::Button("Step")) {
            step();
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

    polyscope::options::autocenterStructures = true;
    polyscope::options::autoscaleStructures = true;
    polyscope::init();
    polyscope::view::style = polyscope::view::NavigateStyle::Planar;
    polyscope::state::userCallback = callback;
    polyscope::show();
    polyscope::screenshot("../outputs/triangulation.tga", false);

    // create suitable matlab arrays for data display purposes
    matlab::data::TypedArray<int> iters = factory.createArray<int>({1, (unsigned long) totalIters + 1});
    matlab::data::TypedArray<double> elapsedTime = factory.createArray<double>({1, (unsigned long) totalIters + 1});
    matlab::data::TypedArray<double> energy = factory.createArray<double>({1, (unsigned long) totalIters + 1});
    for(int i = 0; i <= totalIters; i++) {
        iters[i] = i;
        elapsedTime[i] = elapsedTimeVec.at(i);
        energy[i] = energyVec.at(i);
    }
    matlabPtr->setVariable(u"x", iters);
    matlabPtr->setVariable(u"t", elapsedTime);
    matlabPtr->setVariable(u"E", energy);
    // create and save matlab plots
    matlabPtr->eval(u"f=figure('visible', 'off'); plot(x, t); title('Elapsed Time')");
    matlabPtr->eval(u"exportgraphics(f, '../outputs/data_time.png')");
    matlabPtr->eval(u"g=figure('visible', 'off'); plot(x, E); axis([0 inf 0 inf]); title('Approximation Error')");
    matlabPtr->eval(u"exportgraphics(g, '../outputs/data_energy.png')");

    // convert screenshot sequences to video
    system("ffmpeg -hide_banner -loglevel warning -framerate 2 -i screenshot_%06d.tga -vcodec mpeg4 ../outputs/video_output.mp4");
    system("rm screenshot_*");
	return 0;
}
