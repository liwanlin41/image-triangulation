#include "constant.hpp"

using namespace std;

// identity lambda function
auto identity = [](double x, double y) {return 1.0;};

void ConstantApprox::crossCheck() {
    for(Triangle &t : triangles) {
        cout << t;
        if (trianglesToPoints[&t].size() != 3) cout << "ERROR\n";
        for(Point *ptr : trianglesToPoints[&t]) {
            bool detected = false;
            for(Point &p : points) {
                if (ptr == &p) detected = true;
            }
            bool isCorner = false;
            for(Point &c : corners) {
                if (ptr == &c) isCorner = true;
            }
            if(!detected) {
                cout << *ptr << " of " << t;
            }
            assert(detected || isCorner);
        }
    }
}

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
            Point pixelCorner(-0.5 + maxX * ((i+j)%2), -0.5 + maxY * i);
            corners.push_back(pixelCorner);
        }
    }
    for(int i = 0; i < n; i++) {
        Point p(distanceX * (i+1) - 0.5, centerY);
        points.push_back(p);
    }
    int numTriangles = 0; // to more conveniently get pointers to triangles
    // leave the one or two middle points to be set separately
    for(int i = 0; i < (n-2)/2.0; i++) {
        vector<Point*> uL = {&corners.at(0), &points.at(i), &points.at(i+1)};
        vector<Point*> lL = {&corners.at(3), &points.at(i), &points.at(i+1)};
        vector<Point*> uR = {&corners.at(1), &points.at(n-i-1), &points.at(n-i-2)};
        vector<Point*> lR = {&corners.at(2), &points.at(n-i-1), &points.at(n-i-2)};
        Triangle upperLeft(uL);
        Triangle lowerLeft(lL);
        Triangle upperRight(uR);
        Triangle lowerRight(lR);
        triangles.push_back(upperLeft);
        triangles.push_back(lowerLeft);
        triangles.push_back(upperRight);
        triangles.push_back(lowerRight);
        trianglesToPoints.insert({&triangles[numTriangles], uL});
        trianglesToPoints.insert({&triangles.at(numTriangles+1), lL});
        trianglesToPoints.insert({&triangles.at(numTriangles+2), uR});
        trianglesToPoints.insert({&triangles.at(numTriangles+3), lR});
        numTriangles += 4;
    }
    if (n % 2 == 0) {
        vector<Point*> uL = {&corners.at(0), &points.at((n-2)/2), &points.at(n/2)};
        vector<Point*> uT = {&corners.at(0), &corners.at(1), &points.at(n/2)};
        vector<Point*> lB = {&corners.at(2), &corners.at(3), &points.at((n-2)/2)};
        vector<Point*> lR = {&corners.at(2), &points.at(n/2-1), &points.at(n/2)};
        Triangle upperLeft(uL);
        Triangle upperTop(uT);
        Triangle lowerBottom(lB);
        Triangle lowerRight(lR);
        triangles.push_back(upperLeft);
        triangles.push_back(upperTop);
        triangles.push_back(lowerBottom);
        triangles.push_back(lowerRight);
        trianglesToPoints[&triangles.at(numTriangles)] = uL;
        trianglesToPoints[&triangles.at(numTriangles+1)] = uT;
        trianglesToPoints[&triangles.at(numTriangles+2)] = lB;
        trianglesToPoints[&triangles.at(numTriangles+3)] = lR;
        numTriangles += 4;
    } else {
        vector<Point*> uu = {&corners.at(0), &corners.at(1), &points.at((n-1)/2)};
        vector<Point*> ll = {&corners.at(2), &corners.at(3), &points.at((n-1)/2)};
        Triangle upper(uu);
        Triangle lower(ll);
        triangles.push_back(upper);
        triangles.push_back(lower);
        trianglesToPoints[&triangles.at(numTriangles)] = uu;
        trianglesToPoints[&triangles.at(numTriangles+1)] = ll;
        numTriangles += 2;
    }
    vector<Point*> ll = {&corners.at(0), &corners.at(3), &points.at(0)};
    vector<Point*> rr = {&corners.at(1), &corners.at(2), &points.at(n-1)};
    Triangle left(ll);
    Triangle right(rr);
    triangles.push_back(left);
    triangles.push_back(right);
    trianglesToPoints[&triangles.at(numTriangles)] = ll;
    trianglesToPoints[&triangles.at(numTriangles+1)] = rr;
    numTriangles += 2;
    // for testing: determine where these triangles actually are
    for (map<Triangle*,vector<Point*>>::iterator it = trianglesToPoints.begin(); it != trianglesToPoints.end(); ++it)
    {
        Triangle* t =  it->first;
        cout << t << endl;
    }
    cout << "END TRIANGLES IN MAP\n";
    for(int i = 0; i < numTriangles; i++) {
        cout << &triangles[i] << endl;
    }
    cout << "END TRIANGLES IN VECTOR\n";
    for (auto& t : triangles)
    {
        cout << t;
        cout << &t << endl;
        assert(trianglesToPoints.find(&t) != trianglesToPoints.end());
    }
    // compute initial approximation
    crossCheck();
    updateApprox();
    for(Triangle &t : triangles) {
        cout << t;
        cout << approx[&t] << endl;
    }
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
                double diff = approxVal - (p->getColor());
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
        cout << t;
        vector<Point*> vertices = trianglesToPoints[&t];
        cout << "vertices: ";
        for(Point* p: vertices) {
            cout << *p << " ";
        }
        for(int i = 0; i < 3; i++) {
            double changeX, changeY;
            gradient(t, vertices, i, &changeX, &changeY);
            cout << "succeeded\n";
            gradX[vertices.at(i)] += changeX;
            gradY[vertices.at(i)] += changeY;
        }
    }
}

void ConstantApprox::gradient(Triangle &triangle, vector<Point*> &vertices, int movingPt, double *gradX, double *gradY) {
    // for integrating v dot n in the line integral when pt is moving
    // (x, y) will be on some side of the triangle
    cout << "POINT A\n";
    auto vn = [&vertices, &movingPt](double x, double y, bool dirX) {
        Point current(x, y);
        // to more robustly determine which segment contains (x, y), compute
        // areas of triangles formed by current with triangle sides
        double areas[3];
        for(int i = 0; i < 3; i++) {
            // hold the area to the edge opposite vertices.at(i)
            areas[i] = abs(Triangle::getSignedArea((vertices.at((i+1)%3)), (vertices.at((i+2)%3)), &current));
        }
        double minArea = min(areas[0], min(areas[1], areas[2]));
        if (minArea == areas[movingPt]) { // point is closest to the stationary side
            return 0.0;
        }
        Point segmentEnd = (minArea == areas[(movingPt+1) % 3]) ? *vertices.at((movingPt+2)%3) : *vertices.at((movingPt+1)%3);
        // preserve orietation for outward normal
        Segment gamma = (segmentEnd == *vertices.at((movingPt+1)%3)) ? 
            Segment(vertices.at(movingPt), vertices.at((movingPt+1)%3)) : 
            Segment(vertices.at((movingPt+2)%3), vertices.at((movingPt)));
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
    // integral of fdA
    double imageIntegral = DoubleIntegral::evaluate(identity, image, &triangle);
    double area = triangle.getArea();
    cout << vertices.size() << endl;
    double dA[2] = {triangle.gradX(vertices.at(movingPt)), triangle.gradY(vertices.at(movingPt))};
    cout << "POINT C\n";
    double boundaryChange[2];
    // compute gradient in x direction
    boundaryChange[0] = LineIntegral::evaluate(vnx, image, &triangle);
    // compute gradient in y direction
    boundaryChange[1] = LineIntegral::evaluate(vny, image, &triangle);
    cout << "POINT D\n";
    double gradient[2];
    for(int j = 0; j < 2; j++) {
        gradient[j] = (2 * area * imageIntegral * boundaryChange[j]
            - imageIntegral * imageIntegral * dA[j]) / (area * area);
    }
    // check for null pointers
    if (gradX && gradY) {
        *gradX = gradient[0];
        *gradY = gradient[1];
    }
}

int ConstantApprox::gradUpdate() {
    int badIndex = -1;
    return badIndex;
}

void ConstantApprox::undo(int ind) {

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
        approx[&t] = (val / t.getArea());
    }
}

void ConstantApprox::run(int maxIter, double eps) {

}