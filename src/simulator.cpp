#include "simulator.h"

Simulator::Simulator(const char *imgPath, CImg<unsigned char> *img, ApproxType approxtype) {
    if(approxtype == constant) {
        approx = new ConstantApprox(img, 0.05);
    } else if (approxtype == linear) {
        approx = new LinearApprox(img, 0.05);
    }

    // determine initialization method
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

    matlabPtr = startMATLAB();

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
        approx->initialize(points, edges);
    } else {
        // prompt for dx
        cout << "Sample once every __ pixels? ";
        cin >> dx;
        if(cin.fail()) {
            cin.clear();
            dx = DX_DEFAULT;
            cout << "defaulting to " << DX_DEFAULT << endl;
        }
        approx->initialize(dx);
    }
    cout << "ready\n";
}

Simulator::~Simulator() {
    delete approx;
}

void Simulator::initialize() {
    // complete restart
    elapsedTimeVec.clear();
    energyVec.clear();
    //approx->registerMesh();
    registerMesh(approx);
    // setup gradient descent
    cout << "finding energy..." << endl;
    newEnergy = approx->computeEnergy();
    cout << "done, energy is " << newEnergy << endl;
    // initialize to something higher than newEnergy
    prevEnergy = newEnergy * 2;
    iterCount = 0;
    totalIters = 0;
    elapsedTimeVec.push_back(totalStep); // initial values
    energyVec.push_back(newEnergy);
    // center mesh
    /*
    polyscope::view::resetCameraToHomeView();
    polyscope::resetScreenshotIndex();
    // screenshot
    polyscope::screenshot(false);
    polyscope::screenshot("../outputs/initial.tga", false);
    */
}

void Simulator::step(double eps) {
    cout << "iteration " << iterCount << " (" << totalIters << " total)" << endl;
    // allow a fixed number of energy increases to avoid getting stuck
    totalStep += approx->step(prevEnergy, newEnergy, (totalIters >= 10) && (iterCount >= 5));
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
    updateMesh(approx);
    //approx->updateMesh();
    //polyscope::screenshot(false);
}

bool Simulator::step(int maxIter, double eps) {
    if(iterCount < maxIter && numSmallChanges < maxSmallChanges) {
        step(eps);
        return true;
    } 
    return false;
}

void Simulator::flow(int maxIter, double eps) {
    while(iterCount < maxIter && numSmallChanges < maxSmallChanges) {
        step(eps);
    } 
}

void Simulator::retriangulate(int num) {
    approx->subdivide(num);
    // re-initialize new mesh
    //approx->registerMesh();
    registerMesh(approx);
    // reset values
    iterCount = 0;
    numSmallChanges = 0;
    totalIters++;
    newEnergy = approx->computeEnergy();
    elapsedTimeVec.push_back(totalStep);
    energyVec.push_back(newEnergy);
    cout << "energy after subdivision: " << newEnergy << endl;
    //polyscope::screenshot(false);
}

void Simulator::cleanup() {
    /*
    polyscope::screenshot("../outputs/triangulation.tga", false);
    */
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
    system("rm *.tga");
}