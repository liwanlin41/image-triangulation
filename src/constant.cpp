#include "constant.hpp"

using namespace std;

double tolerance = 1e-10;

// custom rounding function to support needed pixel rounding
int customRound(double x) {
    int floor = (int) x;
    if (abs(x - floor) <= 0.5) {
        return floor;
    }
    else if (x > 0) {
        return floor + 1;
    }
    return floor - 1;
}

// identity lambda function
auto identity = [](double x, double y) {return 1.0;};

ConstantApprox::ConstantApprox(vector<vector<Pixel>> *img, vector<Point> *pts, vector<array<int, 3>> &inds, double step)
    : image(img), points(pts), triangleInd(inds), stepSize(step) {
    maxX = img->size();
    maxY = img->at(0).size();
    // create triangles
    for(int i = 0; i < inds.size(); i++) {
        Triangle t(&points->at(inds.at(i).at(0)), &points->at(inds.at(i).at(1)), &points->at(inds.at(i).at(2)));
        triangles.push_back(t);
    }
    /*
    cout << "POINT ADDRESSES\n";
    for(Point &p : *points) {
        cout << &p << endl;
    }
    cout << "END POINT ADDRESSES\nSTART TRIANGLE ADDRESSES\n";
    for(Triangle &t : triangles) {
        cout << t.vertices.at(0) << " " << t.vertices.at(1) << " " << t.vertices.at(2) << endl;
    }
    cout << "END\n";
    */
    updateApprox();
}

double ConstantApprox::computeEnergy() {
    double energy = 0;
    for(Triangle &t: triangles) {
        double approxVal = approx[&t];
        // compute by iterating over pixels
        for(int x = 0; x < maxX; x++) {
            for(int y = 0; y < maxY; y++) {
                //point to pixel being referenced
                Pixel *p = &(image->at(x).at(y));
                double area = p->intersectionArea(t);
                /*
                if (isnan(area)) {
                    cout << "pixel " << x << ", " << y << endl;
                    cout << "area: " << area << endl;
                }
                */
                double diff = approxVal - (p->getColor());
                /*
                if (isnan(diff)) {
                    cout << "pixel " << x << ", " << y << endl;
                    cout << "diff: " << diff << endl;
                }
                */
                energy += (diff * diff) * area;
            }
        }
    }
    return energy;
}

void ConstantApprox::computeGrad() {
    // clear gradients from last iteration
    for (Point &p : *points) {
        gradX[&p] = 0;
        gradY[&p] = 0;
    }
    for (Triangle &t : triangles) {
        for(int i = 0; i < 3; i++) { 
            double changeX, changeY;
            gradient(t, i, &changeX, &changeY);
            // constrain points on boundary of image
            if(t.vertices.at(i)->isBorderX()) {
                changeX = 0;
            }
            if(t.vertices.at(i)->isBorderY()) {
                changeY = 0;
            }
            gradX[t.vertices.at(i)] += changeX;
            gradY[t.vertices.at(i)] += changeY;
        }
    }
}

void ConstantApprox::gradient(Triangle &triangle, int movingPt, double *gradX, double *gradY) {
    // for integrating v dot n in the line integral when pt is moving
    // (x, y) will be on some side of the triangle
    auto vn = [&triangle, &movingPt](double x, double y, bool dirX) {
        Point current(x, y);
        // to more robustly determine which segment contains (x, y), compute
        // areas of triangles formed by current with triangle sides
        double areas[3];
        for(int i = 0; i < 3; i++) {
            // hold the area to the edge opposite vertices.at(i)
            areas[i] = abs(Triangle::getSignedArea((triangle.vertices.at((i+1)%3)), (triangle.vertices.at((i+2)%3)), &current));
        }
        double minArea = min(areas[0], min(areas[1], areas[2]));
        if (minArea == areas[movingPt]) { // point is closest to the stationary side
            return 0.0;
        }
        Point segmentEnd = (minArea == areas[(movingPt+1) % 3]) ? *triangle.vertices.at((movingPt+2)%3) : *triangle.vertices.at((movingPt+1)%3);
        // preserve orietation for outward normal
        Segment gamma = (segmentEnd == *triangle.vertices.at((movingPt+1)%3)) ? 
            Segment(triangle.vertices.at(movingPt), triangle.vertices.at((movingPt+1)%3)) : 
            Segment(triangle.vertices.at((movingPt+2)%3), triangle.vertices.at((movingPt)));
        // compute velocity at this point by scaling
        double distanceToVertex = current.distance(segmentEnd);
        double scale = distanceToVertex / gamma.length(); // 1 if at a, 0 if at opposite edge
        // determine whether x or y gradient is being computed
        double velX = (dirX) ? scale : 0;
        Matrix v(velX, scale - velX);
        Matrix n = gamma.unitNormal();
        return v.transpose().multiply(n).get(0,0);
    };
    // break down for x, y specifically
    auto vnx = [&vn](double x, double y) {
        return vn(x, y, true);
    };
    auto vny = [&vn](double x, double y) {
        return vn(x, y, false);
    };
    // to save time, only compute integrals if triangle is non-degenerate;
    // degenerate triangle has 0 energy and is locally optimal, set gradient to 0
    double area = triangle.getArea();
    double gradient[2] = {0, 0};
    if (area > tolerance) {
        // integral of fdA
        double imageIntegral = DoubleIntegral::evaluate(identity, image, &triangle);
        double dA[2] = {triangle.gradX(triangle.vertices.at(movingPt)), triangle.gradY(triangle.vertices.at(movingPt))};
        double boundaryChange[2];
        // compute gradient in x direction
        boundaryChange[0] = LineIntegral::evaluate(vnx, image, &triangle);
        // compute gradient in y direction
        boundaryChange[1] = LineIntegral::evaluate(vny, image, &triangle);
        for(int j = 0; j < 2; j++) {
            gradient[j] = (2 * area * imageIntegral * boundaryChange[j]
                - imageIntegral * imageIntegral * dA[j]) / (-area * area);
        }
    }
    // check for null pointers
    if (gradX && gradY) {
        *gradX = gradient[0];
        *gradY = gradient[1];
    }
}

bool ConstantApprox::gradUpdate() {
    // gradient descent update for each point
    for (Point &p : *points) {
        // assert(gradX.find(&p) != gradX.end());
        p.move(-stepSize * gradX.at(&p), -stepSize * gradY.at(&p));
    }
    // now check validity of result
    for (Triangle &t : triangles) {
        if (t.getSignedArea() < 0) {
            return false;
        }
    }
    return true;
}

void ConstantApprox::undo() {
    for (Point &p : *points) {
        p.move(stepSize * gradX.at(&p), stepSize * gradY.at(&p));
    }
    stepSize /= 2;
}

void ConstantApprox::updateApprox() {
    for(Triangle &t : triangles) {
        // cout << t;
        double val = 0;
        // compute total value of image by iterating over pixels
        for(int x = 0; x < maxX; x++) {
            for(int y = 0; y < maxY; y++) {
                Pixel *p = &(image->at(x).at(y));
                double area = p->intersectionArea(t);
                val += area * (p->getColor());
            }
        }
        // take average value
        double approxVal = val / t.getArea();
        /*
        cout << approxVal << endl;
        cout << "area " << t.getArea() << endl;
        cout << "computed area " << A << endl;
        */
        // handle degeneracy
        if (isnan(approxVal)) {
            assert(t.getArea() < tolerance);
            approxVal = 255; // TODO: do something better than this            
        }
        approx[&t] = approxVal;
    }
}

void ConstantApprox::run(int maxIter, double eps) {
    // track change in energy for stopping point
    double newEnergy = computeEnergy();
    // initialize to something higher than newEnergy
    double prevEnergy = newEnergy + 100 * eps;
    int iterCount = 0;
    while(iterCount < maxIter && abs(prevEnergy-newEnergy) > eps) {
        cout << "iteration " << iterCount << endl;
        computeGrad();
        while(!gradUpdate()) {
            undo(); // keep halving stepSize until it works
        }
        updateApprox();
        prevEnergy = newEnergy;
        newEnergy = computeEnergy();
        /*
        while(newEnergy > prevEnergy) { // overshot optimum?
            do {
                undo();
            } while (!gradUpdate());
            updateApprox();
            newEnergy = computeEnergy();
        }
        */
        cout << "new energy: " << newEnergy << endl;
        cout << "step size: " << stepSize << endl;
        iterCount++;
    }
}

double ConstantApprox::getStep() {
    return stepSize;
}

vector<Point> ConstantApprox::getVertices() {
    return *points;
}

vector<array<int,3>> ConstantApprox::getFaces() {
    return triangleInd;
}

vector<array<double,3>> ConstantApprox::getColors() {
    vector<array<double,3>> colors;
    for(int i = 0; i < triangles.size(); i++) {
        // scale to fit polyscope colors TODO: check that this is correct
        double approxColor = approx.at(&triangles.at(i))/255;
        // double grayScale[3] = {approxColor, approxColor, approxColor};
        colors.push_back({approxColor, approxColor, approxColor});
        // cout << approxColor * 255 << " on " << triangles.at(i);
    }
    return colors;
}