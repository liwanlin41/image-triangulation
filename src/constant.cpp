#include "constant.hpp"

using namespace std;

ConstantApprox::ConstantApprox(vector<vector<Pixel>> *img, int n, double step) : image(img), stepSize(step){
    // recall img.at(x).at(y) is pixel (x, y)
    maxX = img->size();
    maxY = img->at(0).size();
    // for now, initialize triangulation using a horizontal line 
    // of points across the middle of the image
    double centerY = maxY / 2 - 0.5;
    double distanceX = (double) maxX / (n+1);
    // corners of image
    Point topLeft(-0.5,-0.5);
    Point bottomLeft(-0.5, -0.5+maxY);
    Point topRight(-0.5+maxX, -0.5);
    Point bottomRight(-0.5+maxX, -0.5+maxY);
    corners.push_back(topLeft);
    corners.push_back(topRight);
    corners.push_back(bottomRight);
    corners.push_back(bottomLeft);
    for(int i = 0; i < n; i++) {
        Point p(distanceX * (i+1) - 0.5, centerY);
        points.push_back(p);
    }
    // leave the one or two middle points to be set separately
    for(int i = 0; i < (n-2)/2.0; i++) {
        Triangle upperLeft(&corners.at(0), &points.at(i), &points.at(i+1));
        Triangle lowerLeft(&corners.at(3), &points.at(i), &points.at(i+1));
        Triangle upperRight(&corners.at(1), &points.at(n-i-1), &points.at(n-i-2));
        Triangle lowerRight(&corners.at(2), &points.at(n-i-1), &points.at(n-i-2));
        triangles.push_back(upperLeft);
        triangles.push_back(lowerLeft);
        triangles.push_back(upperRight);
        triangles.push_back(lowerRight);
    }
    if (n % 2 == 0) {
        Triangle upperLeft(&corners.at(0), &points.at((n-2)/2), &points.at(n/2));
        Triangle upperTop(&corners.at(0), &corners.at(1), &points.at(n/2));
        Triangle lowerBottom(&corners.at(2), &corners.at(3), &points.at((n-2)/2));
        Triangle lowerRight(&corners.at(2), &points.at(n/2-1), &points.at(n/2));
        triangles.push_back(upperLeft);
        triangles.push_back(upperTop);
        triangles.push_back(lowerBottom);
        triangles.push_back(lowerRight);
    } else {
        Triangle upper(&corners.at(0), &corners.at(1), &points.at((n-1)/2));
        Triangle lower(&corners.at(2), &corners.at(3), &points.at((n-1)/2));
        triangles.push_back(upper);
        triangles.push_back(lower);
    }
    Triangle left(&corners.at(0), &corners.at(3), &points.at(0));
    Triangle right(&corners.at(1), &corners.at(2), &points.at(n-1));
    triangles.push_back(left);
    triangles.push_back(right);
    /*
    // for testing: determine where these triangles actually are
    for(Triangle &t : triangles) {
        cout << t;
    }
    */
    // compute initial approximation
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
                double diff = approxVal - (p->getColor());
                energy += (diff * diff) * area;
            }
        }
    }
    return energy;
}

void ConstantApprox::computeGrad() {

}

void ConstantApprox::gradient(Triangle &triangle, Point &pt, double *gradX, double *gradY) {
    auto identity = [](double x, double y) {return 1.0;};
    vector<Point> vertices = triangle.copyVertices();
    int movingInd;
    // determine index of pt
    for(int i = 0; i < 3; i++) {
        if (vertices.at(i) == pt) {
            movingInd = i;
        }
    }
    // for integrating v dot n in the line integral when pt is moving
    // (x, y) will be on some side of the triangle
    auto vn = [&vertices, &movingInd, &pt](double x, double y, bool dirX) {
        Point current(x, y);
        // to more robustly determine which segment contains (x, y)
        double areas[3];
        for(int i = 0; i < 3; i++) {
            // hold the area to the edge opposite vertices.at(i)
            areas[i] = abs(Triangle::getSignedArea(&(vertices.at((i+1)%3)), &(vertices.at((i+2)%3)), &current));
        }
        double minArea = min(areas[0], min(areas[1], areas[2]));
        if (minArea == areas[movingInd]) { // point is closest to the stationary side
            return 0.0;
        }
        Point segmentEnd = (minArea == areas[(movingInd+1) % 3]) ? vertices.at((movingInd+2)%3) : vertices.at((movingInd+1)%3);
        // preserve orietation for outward normal
        Segment gamma = (segmentEnd == vertices.at((movingInd+1)%3)) ? 
            Segment(&pt, &vertices.at((movingInd+1)%3)) : Segment(&vertices.at((movingInd+2)%3), &pt);
        // compute velocity at this point by scaling
        double distanceToVertex = current.distance(segmentEnd);
        double scale = distanceToVertex / gamma.length(); // 1 if at a, 0 if at opposite edge
        // determine whether x or y gradient is being computed
        double velX = (dirX) ? scale : 0;
        Matrix v(velX, scale - velX);
        Matrix n = gamma.unitNormal();
        return v.transpose().multiply(n).get(0,0);
    };
    auto vnx = [&vn](double x, double y) {
        return vn(x, y, true);
    };
    auto vny = [&vn](double x, double y) {
        return vn(x, y, false);
    };
    // integral of fdA
    double imageIntegral = DoubleIntegral::evaluate(identity, image, &triangle);
    double area = triangle.getArea();
    double dA[2] = {triangle.gradX(&pt), triangle.gradY(&pt)};
    double boundaryChange[2];
    // compute gradient in x direction
    boundaryChange[0] = LineIntegral::evaluate(vnx, image, &triangle);
    // compute gradient in y direction
    boundaryChange[1] = LineIntegral::evaluate(vny, image, &triangle);
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
    /*
    double val = 0;
    Pixel *p = &(image->at(15).at(49));
    for(Triangle &t : triangles) {
        cout << t;
        val += p->intersectionArea(t);
    }
    */
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