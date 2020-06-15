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

ConstantApprox::ConstantApprox(vector<vector<Pixel>> *img, int n, double step) : image(img), stepSize(step){
    // recall img.at(x).at(y) is pixel (x, y)
    maxX = img->size();
    maxY = img->at(0).size();
    // for now, initialize triangulation using a horizontal line 
    // of points across the middle of the image
    double centerY = maxY / 2 - 0.5;
    double distanceX = (double) maxX / (n+1);
    // corners of image
    for(int i = 0; i < 2; i++) {
        for(int j = 0; j < 2; j++) {
            Point pixelCorner(-0.5 + maxX * ((i+j)%2), -0.5 + maxY * i, true, true);
            points.push_back(pixelCorner);
        }
    }
    for(int i = 0; i < n; i++) {
        Point p(distanceX * (i+1) - 0.5, centerY);
        points.push_back(p);
    }
    // leave the one or two middle points to be set separately
    for(int i = 0; i < (n-2)/2.0; i++) {
        Triangle upperLeft(&points.at(0), &points.at(i+4), &points.at(i+5));
        Triangle lowerLeft(&points.at(3), &points.at(i+4), &points.at(i+5));
        Triangle upperRight(&points.at(1), &points.at(n-i+3), &points.at(n-i+2));
        Triangle lowerRight(&points.at(2), &points.at(n-i+3), &points.at(n-i+2));
        triangles.push_back(upperLeft);
        triangles.push_back(lowerLeft);
        triangles.push_back(upperRight);
        triangles.push_back(lowerRight);
        /*
        int uL[3] = {0, i+4, i+5};
        int lL[3] = {3, i+4, i+5};
        int uR[3] = {1, n-i+3, n-i+2};
        int lR[3] = {2, n-i+3, n-i+2};
        */
        triangleInd.push_back({0, i+4, i+5});
        triangleInd.push_back({3, i+4, i+5});
        triangleInd.push_back({1, n-i+3, n-i+2});
        triangleInd.push_back({2, n-i+3, n-i+2});
    }
    if (n % 2 == 0) {
        Triangle upperLeft(&points.at(0), &points.at(n/2+3), &points.at(n/2+4));
        Triangle upperTop(&points.at(0), &points.at(1), &points.at(n/2+4));
        Triangle lowerBottom(&points.at(2), &points.at(3), &points.at(n/2+3));
        Triangle lowerRight(&points.at(2), &points.at(n/2+3), &points.at(n/2+4));
        triangles.push_back(upperLeft);
        triangles.push_back(upperTop);
        triangles.push_back(lowerBottom);
        triangles.push_back(lowerRight);
        /*
        int uL[3] = {0, (n-2)/2+4, n/2+4};
        int uT[3] = {0, 1, n/2+4};
        int lB[3] = {2, 3, n/2+3};
        int lR[3] = {2, n/2+3, n/2+4};
        */
        triangleInd.push_back({0, n/2+3, n/2+4});
        triangleInd.push_back({0, 1, n/2+4});
        triangleInd.push_back({2, 3, n/2+3});
        triangleInd.push_back({2, n/2+3, n/2+4});
    } else {
        Triangle upper(&points.at(0), &points.at(1), &points.at((n-1)/2+4));
        Triangle lower(&points.at(2), &points.at(3), &points.at((n-1)/2+4));
        triangles.push_back(upper);
        triangles.push_back(lower);
        // int up[3] = {0, 1, (n-1)/2+4};
        // int lo[3] = {2, 3, (n-1)/2+4};
        triangleInd.push_back({0, 1, (n-1)/2+4});
        triangleInd.push_back({2, 3, (n-1)/2+4});
    }
    Triangle left(&points.at(0), &points.at(3), &points.at(4));
    Triangle right(&points.at(1), &points.at(2), &points.at(n+3));
    triangles.push_back(left);
    triangles.push_back(right);
    // int le[3] = {0, 3, 4};
    // int ri[3] = {1, 2, n+3};
    triangleInd.push_back({0, 3, 4});
    triangleInd.push_back({1, 2, n+3});
    // for testing: determine where these triangles actually are
    // compute initial approximation
    updateApprox();
}

double ConstantApprox::computeEnergy() {
    double energy = 0;
    for(Triangle &t: triangles) {
        double approxVal = approx[&t];
        /*
        if (isnan(approxVal)) {
            cout << t;
            throw domain_error("broken");
        }
        */
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
    for (Point &p : points) {
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
    for (Point &p : points) {
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
    for (Point &p : points) {
        p.move(stepSize * gradX.at(&p), stepSize * gradY.at(&p));
    }
    stepSize /= 2;
}

void ConstantApprox::updateApprox() {
    for(Triangle &t : triangles) {
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

vector<Point> ConstantApprox::getVertices() {
    return points;
}

vector<array<int,3>> ConstantApprox::getFaces() {
    return triangleInd;
}

vector<array<double,3>> ConstantApprox::getColors() {
    vector<array<double,3>> colors;
    for(int i = 0; i < triangles.size(); i++) {
        double approxColor = approx.at(&triangles.at(i));
        // double grayScale[3] = {approxColor, approxColor, approxColor};
        cout << approxColor << " on " << triangles.at(i);
        colors.push_back({approxColor, approxColor, approxColor});
        for(int ind : triangleInd[i]) {
            cout << ind << " ";
        }
        cout << endl;
    }
    return colors;
}

/*
CImg<unsigned char> ConstantApprox::show() {
    // if unassigned, fill with 0
    CImg<unsigned char> result(maxX, maxY, 1, 1, 0);
    for(Triangle &t : triangles) {
        cout << t;
        int coords[6];
        for(int i = 0; i < 3; i++) {
            coords[2 * i] = customRound(t.vertices.at(i)->getX());
            coords[2 * i + 1] = customRound(t.vertices.at(i)->getY());
        }
        int approxColor = (int) approx.at(&t);
        cout << "color: " << approxColor << endl;
        unsigned char color[] = {approxColor, approxColor, approxColor};
        result.draw_triangle(coords[0], coords[1], coords[2], coords[3],
            coords[4], coords[5], color, 1);
    }
    cout << "stepsize: " << stepSize << endl;
    return result;
}
*/