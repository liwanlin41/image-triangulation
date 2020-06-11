#include "lineIntegral.hpp"


// helper for extracting grayscale

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

LineIntegral::LineIntegral() {}

// double LineIntegral::evaluate(double (*func)(double, double), CImg<unsigned char> *img, Triangle *triangle) {
double LineIntegral::evaluate(double (*func)(double, double), Triangle *triangle) {
    int xMin, xMax, yMin, yMax;
    // compute bounding box of triangle, in pixels
    triangle->boundingBox(xMin, xMax, yMin, yMax);
    // iterate over bounding box
    vector<Point> vertices = triangle->copyVertices();
    // preserve ccw direction
    Segment seg01(&vertices.at(0), &vertices.at(1));
    Segment seg12(&vertices.at(1), &vertices.at(2));
    Segment seg20(&vertices.at(2), &vertices.at(0));
    double integral = 0;
    for(int i = xMin; i <= xMax; i++) {
        for(int j = yMin; j <= yMax; j++) {
            Pixel p(i, j, 1);
            integral += p.lineIntegral(func, seg01) + p.lineIntegral(func, seg12) + p.lineIntegral(func, seg20);
        }
    }
    return integral;
}

// double LineIntegral::evaluate(double (*func)(double, double), CImg<unsigned char> *img, Point *a, Point *b, Point *c) {
double LineIntegral::evaluate(double (*func)(double, double), Point *a, Point *b, Point *c) {
    // first compute bounding box
    double ax = a->getX();
    double ay = a->getY();
    double bx = b->getX();
    double by = b->getY();
    double cx = c->getX();
    double cy = c->getY();
    int xMin = round(min(ax, min(bx, cx)));
    int xMax = round(max(ax, max(bx, cx)));
    int yMin = round(min(ay, min(by, cy)));
    int yMax = round(max(ay, max(by, cy)));
    double integral = 0;
    // ensure order of points within segment is correct
    Segment ab(a, b);
    Segment bc(b, c);
    Segment ca(c, a);
    for(int i = xMin; i <= xMax; i++) {
        for(int j = yMin; j <= yMax; j++) {
            // set color to 1 just to test if this method works
            Pixel p(i, j, 1);
            integral += p.lineIntegral(func, ab) + p.lineIntegral(func, bc) + p.lineIntegral(func, ca);
        }
    }
    return integral;
}