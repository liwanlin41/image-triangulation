#include "constant.hpp"

// gradient when moving point pt at velocity (vx, vy)
double linearGradient(Triangle &triangle, Point &pt, double vx, double vy, vector<vector<Pixel>> &image) {
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
    double imageIntegral = DoubleIntegral::evaluate(identity, &image, &triangle);
    double area = triangle.getArea();
    double dA[2] = {triangle.gradX(&pt), triangle.gradY(&pt)};
    double boundaryChange[2];
    // compute gradient in x direction
    boundaryChange[0] = LineIntegral::evaluate(vnx, &image, &triangle);
    // compute gradient in y direction
    boundaryChange[1] = LineIntegral::evaluate(vny, &image, &triangle);
    double gradient[2];
    for(int j = 0; j < 2; j++) {
        gradient[j] = (2 * area * imageIntegral * boundaryChange[j]
            - imageIntegral * imageIntegral * dA[j]) / (area * area);
    }
    double gradApprox = gradient[0] * vx + gradient[1] * vy;
    return gradApprox;
}