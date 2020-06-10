#include "pixel.hpp"

Pixel::Pixel(int x_, int y_, int c) : x(x_), y(y_), color(c) {
    // confirmed: this is cast correctly 
    corners.push_back(Point(x-0.5,y-0.5));
    corners.push_back(Point(x+0.5,y-0.5));
    corners.push_back(Point(x+0.5,y+0.5));
    corners.push_back(Point(x-0.5,y+0.5));
}

int Pixel::getColor() {
    return color;
}

double Pixel::intersectionLength(Segment &e) {
    vector<Point> intersections(2); // hold intersections
    for(int i = 0; i < 4; i++) {
        // retrieve a side of the pixel; at most two will have an intersection
        Segment side(&corners.at(i), &corners.at((i+1) % 4));
    }
}