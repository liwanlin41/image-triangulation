#include <gtest/gtest.h>
#include <math.h>
#include "../src/triangle.hpp"
#include "../src/lineIntegral.hpp"
#include "../src/doubleIntegral.hpp"
#include "../src/matrix.hpp"
using namespace std;

const double eps = 0.001;
const int num_iter = 1000;
const double tolerance = 1e-6;

// for testing purposes, create an 100 x 100 image all of whose pixels have color 1
vector<vector<Pixel>> generateFakeImage() {
    vector<vector<Pixel>> image;
    for(int i = 0; i < 100; i++) {
        vector<Pixel> holder;
        for(int j = 0; j < 100; j++) {
            holder.push_back(Pixel(i, j, 1));
        }
        image.push_back(holder);
    }
    return image;
}

vector<vector<Pixel>> image = generateFakeImage();

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
    cout << "area: " << area << endl;
    cout << "imageIntegral: " << imageIntegral << endl; // this value is incorrect
    cout << "dA: " << dA[0] << ", " << dA[1] << endl;
    cout << "boundaryChange values " << boundaryChange[0] << ", " << boundaryChange[1] << endl;
    cout << "gradient values: " << gradient[0] << ", " << gradient[1] << endl;
    return gradApprox;
}

TEST(DerivativeTest, FixedTriangle) {
    double vx = 2 * ((double) rand() / RAND_MAX) - 1;
    // pick unit vector
    double vy = pow((1 - vx * vx), 0.5);
    // randomly pick sign
    if ((double) rand() / RAND_MAX < 0.5) {
        vy *= -1;
    }
    Point a(0,0);
    Point b(5,0);
    Point c(3,4);
    Triangle triangle(&a, &b, &c);
    double gradX = triangle.gradX(&c);
    double gradY = triangle.gradY(&c);
    // simulate eps timestep
    double curArea = triangle.getSignedArea();
    c.move(eps * vx, eps * vy);
    double futureArea = triangle.getSignedArea();
    c.move(-2 * eps * vx, -2 * eps * vy);
    double pastArea = triangle.getSignedArea();
    double change1 = (futureArea - curArea) / eps;
    double change2 = (curArea - pastArea) / eps;
    ASSERT_TRUE(abs(change1 - change2) < tolerance);
    double finiteDiff = (futureArea - pastArea) / (2 * eps);
    double gradApprox = (gradX * vx + gradY * vy);
    ASSERT_TRUE(abs(gradApprox - finiteDiff) < tolerance);
}

TEST(DerivativeTest, TriangleArea) {
    for(int i = 0; i < num_iter; i++) {
        // pick a random number from -1 to 1 for velocity
        double vx = 2 * ((double) rand() / RAND_MAX) - 1;
        double vy = 2 * ((double) rand() / RAND_MAX) - 1;
        // generate some random triangle with coordinates < 100
        double random_coords[6];
        for(int j = 0; j < 6; j++) {
            random_coords[j] = (double) rand() / RAND_MAX * 100;
        }
        Point a(random_coords[0], random_coords[1]);
        Point b(random_coords[2], random_coords[3]);
        Point c(random_coords[4], random_coords[5]);
        Segment e(&a, &b);
        Triangle triangle(&a, &b, &c);
        double gradX = triangle.gradX(&a);
        double gradY = triangle.gradY(&a);
        // simulate area(x + eps * v), will have to undo later
        a.move(eps * vx, eps * vy);
        double area1 = triangle.getSignedArea();
        a.move(- 2 * eps * vx, - 2 * eps * vy);
        double area2 = triangle.getSignedArea();
        double areaDiff = (area1 - area2) / (2 * eps);
        double gradApprox = gradX * vx + gradY * vy;
        ASSERT_TRUE(abs(areaDiff - gradApprox) < tolerance);
    }
}

TEST(ConstantTest, FixedTriangle) {
    auto identity = [](double x, double y) {return 1.0;};
    Point a(0,0);
    Point b(3,0);
    Point c(0,4);
    Triangle triangle(&a, &b, &c);
    // move vertex c
    double vx = 0.6;
    double vy = 0.8;
    auto vn = [&a, &b, &c](double x, double y, bool dirX) {
        Point current(x, y);
        // to more robustly determine which segment contains (x, y)
        double areaAB = Triangle::getSignedArea(&a, &b, &current);
        double areaBC = Triangle::getSignedArea(&b, &c, &current);
        double areaCA = Triangle::getSignedArea(&c, &a, &current);
        double minArea = min(areaAB, min(areaBC, areaCA));
        if (minArea == areaAB) { // point is closest to side AB, the stationary side
            return 0.0;
        }
        Point segmentEnd = (minArea == areaCA) ? a : b;
        double distanceToVertex = current.distance(segmentEnd);
        Segment gamma = (segmentEnd == a) ? Segment(&c, &a) : Segment(&b, &c);
        double scale = distanceToVertex / gamma.length();
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
    double dA[2] = {triangle.gradX(&c), triangle.gradY(&c)};
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
    
    // now compute central finite difference
    c.move(eps * vx, eps * vy);
    double futureImgInt = DoubleIntegral::evaluate(identity, &image, &triangle);
    double futureArea = triangle.getArea();
    double futureEnergy = futureImgInt * futureImgInt / futureArea;
    
    c.move(-2 * eps * vx, -2 * eps * vy);
    double pastImgInt = DoubleIntegral::evaluate(identity, &image, &triangle);
    double pastArea = triangle.getArea();
    double pastEnergy = pastImgInt * pastImgInt / pastArea;
    
    double finiteApprox = (futureEnergy - pastEnergy) / (2 * eps);
    ASSERT_TRUE(abs(gradApprox - finiteApprox) < tolerance);
}

// test broken case
TEST(ConstantTest, Bugged) {
    auto identity = [](double x, double y) {return 1.0;};
    Point a(46.8693, 78.12); 
    Point b(14.6368, 45.3289); 
    Point c(42.7346, 73.5944);
    Triangle triangle(&a, &b, &c);
    double vx = 0.353958; 
    double vy = -0.875836;
    double gradApprox = linearGradient(triangle, a, vx, vy, image);
    // now compute central finite difference
    a.move(eps * vx, eps * vy);
    double futureImgInt = DoubleIntegral::evaluate(identity, &image, &triangle);
    double futureArea = triangle.getArea();
    double futureEnergy = futureImgInt * futureImgInt / futureArea;
    
    a.move(-2 * eps * vx, -2 * eps * vy);
    double pastImgInt = DoubleIntegral::evaluate(identity, &image, &triangle);
    double pastArea = triangle.getArea();
    double pastEnergy = pastImgInt * pastImgInt / pastArea;
    
    double finiteApprox = (futureEnergy - pastEnergy) / (2 * eps);
    ASSERT_FLOAT_EQ(finiteApprox, gradApprox);
}

TEST(ConstantTest, Function) {
    auto identity = [](double x, double y) {return 1.0;};
    Point a(0,0);
    Point b(3,0);
    Point c(0,4);
    Triangle triangle(&a, &b, &c);
    // move vertex c
    double vx = 0.6;
    double vy = 0.8;
    double gradApprox = linearGradient(triangle, c, vx, vy, image);
    cout << gradApprox << endl;

    // now compute central finite difference
    c.move(eps * vx, eps * vy);
    double futureImgInt = DoubleIntegral::evaluate(identity, &image, &triangle);
    double futureArea = triangle.getArea();
    double futureEnergy = futureImgInt * futureImgInt / futureArea;
    
    c.move(-2 * eps * vx, -2 * eps * vy);
    double pastImgInt = DoubleIntegral::evaluate(identity, &image, &triangle);
    double pastArea = triangle.getArea();
    double pastEnergy = pastImgInt * pastImgInt / pastArea;
    
    double finiteApprox = (futureEnergy - pastEnergy) / (2 * eps);
    ASSERT_TRUE(abs(gradApprox - finiteApprox) < tolerance);
}

// test derivative for energy function with constant approximation
// over a single triangle
TEST(ConstantTest, ConstantApprox) {
    // for just integrating the given function
    auto identity = [](double x, double y) {return 1.0;};
    for(int i = 0; i < num_iter; i++) {
        cout << "iteration " << i << endl;
        // pick a random number from -1 to 1 for velocity directions
        // for testing, point a will be moved
        double vx = 2 * ((double) rand() / RAND_MAX) - 1;
        double vy = 2 * ((double) rand() / RAND_MAX) - 1;
        // generate some random triangle with coordinates < 100
        double random_coords[6];
        for(int j = 0; j < 6; j++) {
            random_coords[j] = (double) rand() / RAND_MAX * 99;
        }
        Point a(random_coords[0], random_coords[1]);
        Point b(random_coords[2], random_coords[3]);
        Point c(random_coords[4], random_coords[5]);
        // orientation
        if (Triangle::getSignedArea(&a, &b, &c) < 0) {
            Point temp = a;
            a = c;
            c = temp;
        }
        Triangle triangle(&a, &b, &c);
        // for integrating v dot n in the line integral when a is moving
        // at a velocity of (1, 0);
        // (x, y) will be on some side of the triangle
        auto vn = [&a, &b, &c](double x, double y, bool dirX) {
            Point current(x, y);
            // to more robustly determine which segment contains (x, y),
            // check which segment (x, y) is closest to
            double areaAB = Triangle::getSignedArea(&a, &b, &current);
            double areaBC = Triangle::getSignedArea(&b, &c, &current);
            double areaCA = Triangle::getSignedArea(&c, &a, &current);
            double minArea = min(areaAB, min(areaBC, areaCA));
            if (minArea == areaBC) { // point is closest to side BC, the stationary side
                return 0.0;
            }
            Point segmentEnd = (minArea == areaAB) ? b : c;
            // compute velocity at this point by scaling
            double distanceToVertex = current.distance(segmentEnd);
            // preserve orientation for outward normal
            Segment gamma = (segmentEnd == b) ? Segment(&a, &b) : Segment(&c, &a);
            double scale = distanceToVertex / gamma.length(); // 1 if at a, 0 if at opposite edge
            // determine whether x or y gradient is being computed
            double velX = (dirX) ? scale : 0;
            Matrix v(velX, scale - velX);
            Matrix n = gamma.unitNormal();
            return v.transpose().multiply(n).get(0,0);
        };
        // break down into x and y components
        auto vnx = [&vn](double x, double y) {
            return vn(x, y, true);
        };
        auto vny = [&vn](double x, double y) {
            return vn(x, y, false);
        };
        // integral of fdA
        double imageIntegral = DoubleIntegral::evaluate(identity, &image, &triangle);
        double area = triangle.getArea();
        double dA[2] = {triangle.gradX(&a), triangle.gradY(&a)};
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

        // now compute central finite difference
        a.move(eps * vx, eps * vy);
        double futureImgInt = DoubleIntegral::evaluate(identity, &image, &triangle);
        double futureArea = triangle.getArea();
        double futureEnergy = futureImgInt * futureImgInt / futureArea;

        a.move(-2 * eps * vx, -2 * eps * vy);
        double pastImgInt = DoubleIntegral::evaluate(identity, &image, &triangle);
        double pastArea = triangle.getArea();
        double pastEnergy = pastImgInt * pastImgInt / pastArea;

        double finiteApprox = (futureEnergy - pastEnergy) / (2 * eps);
        if (abs(finiteApprox - gradApprox) > 0.1) {
            cout << "vertices: (" << a.getX() << ", " << a.getY() << "), (" << b.getX() << ", " << b.getY() << "), (" << c.getX() << ", " << c.getY() << ")\n";
            cout << "velocity: " << vx << ", " << vy << endl;
            cout << "imageIntegral " << imageIntegral << endl;
            cout << "area: " << area << endl;
            cout << "dA: " << dA[0] << ", " << dA[1] << endl;
            cout << "boundaryChange values " << boundaryChange[0] << ", " << boundaryChange[1] << endl;
            cout << "gradient values: " << gradient[0] << ", " << gradient[1] << endl;
        }
        ASSERT_FLOAT_EQ(finiteApprox, gradApprox);
    }
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
