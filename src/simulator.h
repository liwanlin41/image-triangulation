#ifndef simulator_h
#define simulator_h

#include <iostream>
#include "CImg.h"
#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"
/*
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"
#include "polyscope/view.h"
*/
#include "constant.h"
#include "linear.h"

/**
 * run gradient flow and subdivisions on a triangular mesh approximation
 */

using namespace cimg_library;
using namespace matlab::engine;

// let polyscope read values from point
/*
double adaptorF_custom_accessVector2Value(const Point& p, unsigned int ind) {
    // reflect everything so that it displays correctly
    if (ind == 0) return -p.getX(); 
    if (ind == 1) return -p.getY();
    throw std::logic_error("bad access");
    return -1.;
}
*/

// rendering functions that are currently implemented in main file

// render mesh with colors
// first indicates whether this is the first registration of the mesh
void registerMesh(Approx *approx);
void updateMesh(Approx *approx);

class Simulator {
    private:
        static constexpr double DENSITY_DEFAULT = 0.05; // experimentally a decent starting density
        static const int DX_DEFAULT = 50; // sampling default

        std::unique_ptr<MATLABEngine> matlabPtr;
        matlab::data::ArrayFactory factory;

        Approx *approx; // actual approximation instance

        // initialization values
        double density; // density input for TRIM
        int dx; // take one sample every dx pixels (uniform)

        // for rendering
        double prevEnergy = 0;
        double newEnergy = 0;
        int iterCount = 0;
        int totalIters = 0; // total iterations
        vector<double> elapsedTimeVec; // hold cumulative step size
        vector<double> energyVec; // hold energy per iteration
        vector<double> errorVec; // hold approximation error per iteration
        double totalStep = 0; // running total step

        // for controlling convergence 
        int numSmallChanges = 0;
        static const int maxSmallChanges = 3;

    public:
        // create an approximation instance from an image path
        Simulator(const char *imgPath, CImg<unsigned char> *img, ApproxType approxtype);
        ~Simulator();

        // display the starting triangulation
        void initialize();
        // run one step of gradient flow, incrementing numSmallChanges
        // if a change within a multiplicative factor of eps is detected
        void step(double eps);
        // run one step of gradient flow with stopping conditions;
        // return true if a step was made
        bool step(int maxIter = 100, double eps = 0.001);
        // run gradient flow to convergence
        void flow(int maxIter = 100, double eps = 0.001);

        // handle retriangulation of num edges
        void retriangulate(int num);
        
        // handle output
        void cleanup();
};

#endif