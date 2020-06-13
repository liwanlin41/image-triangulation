#include <assert.h>
#include "pixel.hpp"

// helper functions

// determine whether two points are "essentially" equal (floating point error)
bool approxEqual(Point &a, Point &b, double tolerance = 1e-12) {
    return (a.distance(b) < tolerance);
}

// check to ensure points is a valid ccw polygon
bool isCCW(vector<Point> &points) {
    double centroidX = 0;
    double centroidY = 0;
    int n = points.size();
    for(int i = 0; i < n; i++) {
        centroidX += points.at(i).getX();
        centroidY += points.at(i).getY();
    }
    centroidX /= n;
    centroidY /= n;
    Point centroid(centroidX, centroidY);
    for(int i = 0; i < n; i++) {
        if (Triangle::getSignedArea(&centroid, &points.at(i), &points.at((i+1)%n)) < 0) {
            cout << points.at(i) << "--" << points.at((i+1)%n) << " is bad\n";
            return false;
        }
    }
    return true;
}

// compute (unsigned) area of the polygon enclosed by points,
// where edges of the polygon are given by points.at(i) -- points.at(i+1)
double shoelace(vector<Point> &points) {
    if (points.size() < 3) {
        return 0;
    }
    double area = 0;
    int n = points.size();
    for(int i = 0; i < n; i++) {
        double x0 = points.at(i).getX();
        double y0 = points.at(i).getY();
        double x1 = points.at((i+1) % n).getX();
        double y1 = points.at((i+1) % n).getY();
        area += (x0 * y1 - x1 * y0);
    }
    // in practice points is supposed to be ccw up to floating point errors
    // that do not affect area
    /*
    if (!isCCW(points)) {
        cout << "POINTS START" << endl;
        for(Point pt : points) {
            cout << pt << endl;
        }
        cout << "POINTS END" << endl;
    }
    */
    assert(area >= 0);
    /*
    // account for sign
    if (area < 0) area *= -1;
    */
    // don't forget to divide by 2
    return area/2;
}

// compute integral of x over polygon points and store in totalX, similar for y
// center is a reference point inside the pixel; even if it lies outside the polygon,
// using signed areas means the result will still be correct
void integrateXY(double *totalX, double *totalY, vector<Point> &points, Point &center) {
    double sumX = 0;
    double sumY = 0;
    int n = points.size();
    for(int i = 0; i < n; i++) {
        // average over triangle is just the centroid
        double centroidX = (points.at(i).getX() + points.at((i+1) % n).getX() + center.getX())/3;
        double centroidY = (points.at(i).getY() + points.at((i+1) % n).getY() + center.getY())/3;
        double triangleArea = Triangle::getSignedArea(&center, &points.at(i), &points.at((i+1) % n));
        // weight the average
        sumX += centroidX * triangleArea;
        sumY += centroidY * triangleArea;
    }    
    *totalX = sumX;
    *totalY = sumY;
}

// compute average values of x, y over the polygon enclosed by points
// and put them in the given variables
// center is a reference point
void averageXY(double *avgX, double *avgY, vector<Point> &points, Point &center) {
    double totalX;
    double totalY;
    integrateXY(&totalX, &totalY, points, center);
    double totalArea = shoelace(points);
    *avgX = totalX / totalArea;
    *avgY = totalY / totalArea;
}

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

bool Pixel::containsPoint(Point &p) {
    double px = p.getX();
    double py = p.getY();
    return (-0.5 + x <= px && px <= 0.5 + x) && (-0.5 + y <= py && py <= 0.5 + y);
}

double Pixel::intersectionLength(Segment &e, double *xVal, double *yVal) {
    vector<Point> intersections; // hold intersections
    for(int i = 0; i < 4; i++) {
        // retrieve a side of the pixel; at most two will have an intersection
        // unless intersection is at corners
        Segment side(&corners.at(i), &corners.at((i+1) % 4));
        if (side.intersects(e)) {
            Point intersectionPoint = side.getIntersection(e);
            bool isNewPoint = true; // whether this intersection is a new distinct point
            for(Point pt : intersections) {
                if (pt == intersectionPoint) {
                    isNewPoint = false;
                }
            }
            if (isNewPoint) {
                intersections.push_back(intersectionPoint);
            }
        }
    }
    if (intersections.size() < 2) {
        Point start = e.getStart(); // creates a copy
        Point end = e.getEnd();
        if (containsPoint(start)) {
            intersections.push_back(start);
        }
        if (containsPoint(end)) {
            intersections.push_back(end);
        }
    }
    if (intersections.size() < 2) {
        return 0;
    }
    if (intersections.size() > 2) {
        cout << "BAD POINTS\n";
        for (Point pt : intersections) {
            cout << pt << endl;
        }
        cout << "END BAD POINTS\n";
    }
    assert(intersections.size() == 2);
    Segment contained(&intersections.at(0), &intersections.at(1));
    // check for null pointers, assign midpoint coords
    if (xVal && yVal) {
        *xVal = (intersections.at(0).getX() + intersections.at(1).getX())/2;
        *yVal = (intersections.at(0).getY() + intersections.at(1).getY())/2;
    }
    return contained.length();
}

double Pixel::lineIntegral(function<double(double, double)> func, Segment &e) {
    double midX, midY;
    double length = intersectionLength(e, &midX, &midY);
    if (length == 0) {
        return 0;
    }
    return func(midX, midY) * length * color;
}

double Pixel::intersectionArea(Triangle &t, vector<Point> *polygon) {
    vector<Point> boundary; // hold vertices of polygon formed by intersection
    vector<Point> triangleVertices = t.copyVertices();
    vector<Segment> triangleSides; // hold sides of triangle
    int inInd; // index of some triangle vertex that lies inside pixel (may not exist)
    for(int i = 0; i < 3; i++) {
        triangleSides.push_back(Segment(&triangleVertices.at(i),&triangleVertices.at((i+1) % 3)));
        // add triangle vertices which may be inside the pixel, but don't add corners
        bool isCorner = false;
        for(int j = 0; j < 4; j++) {
            if (triangleVertices.at(i) == corners.at(j)) isCorner = true;
        }
        if (!isCorner && containsPoint(triangleVertices.at(i))) {
            // cout << triangleVertices.at(i) << " is not a corner\n";
            inInd = i;
            boundary.push_back(triangleVertices.at(i));            
        }
    }
    // determine corner to start so as to preserve ccw property
    int start = 0;
    // do this by starting from a corner outside the triangle (if it exists);
    // if it doesn't exist start will stay at 0
    for(int i = 0; i < 4; i++) {
        // additionally, if there is exactly one point inside the triangle, make sure to start
        // at a corner on the same side of the interior point so that the first edge
        // interior point -- intersection point is correct (avoid issues of pixel corners inside
        // the triangle being non-adjacent)
        bool safelyOriented = (boundary.size() != 1) || 
            (Triangle::getSignedArea(&corners.at(i), &triangleVertices.at((inInd+1)%3), &triangleVertices.at((inInd+2)%3)) >= 0);
        if (!t.contains(corners.at(i)) && safelyOriented) {
            start = i;
        }
    }
    for(int i = 0; i < 4; i++) {
        // first determine if corner of pixel is inside
        Point corner = corners.at((i+start) % 4);
        Segment side(&corner, &corners.at((i+start+1)%4));
        if (t.contains(corner)) {
            // cout << "corner " << corner << " added in " << x << ", " << y << endl;
            boundary.push_back(corner);
        }
        // determine intersections with side (i, i+1)
        vector<Point> sideIntersections;
        for(Segment e : triangleSides) {
            if (side.intersects(e)) {
                Point intersectionPoint = side.getIntersection(e);
                // check to see if this point is already accounted for by corners
                // or by triangle vertices; if it isn't exactly equal it won't contribute to area
                // (and the lack of exact equality is likely due to floating point error)
                if (!approxEqual(intersectionPoint, corner) && !approxEqual(intersectionPoint, corners.at((i+start+1)%4))) {
                    bool isVertex = false;
                    for(Point tVertex : triangleVertices) {
                        if (approxEqual(intersectionPoint, tVertex)) {
                            isVertex = true;
                        }
                    }
                    if (!isVertex) {
                        sideIntersections.push_back(intersectionPoint);
                    }
                }
            }
        }
        // note a triangle can intersect a given side at most twice
        assert(sideIntersections.size() <= 2);
        /*
        if (sideIntersections.size() > 0) {
            cout << "points of intersection\n";
            for(Point inter : sideIntersections) {
                cout << inter << endl;
            }
            cout << "end\n";
        }
        */
        // handle normal case where there is only one intersection with this side
        if (sideIntersections.size() == 1) {
            boundary.push_back(sideIntersections.at(0));
        } else if (sideIntersections.size() == 2) {
            Point center(this->x, this->y); // center of this pixel
            double signedArea = Triangle::getSignedArea(&center, &sideIntersections.at(0), &sideIntersections.at(1));
            // if signedArea == 0, sideIntersections must contain two of the same point
            // which means one vertex of the triangle is on the side; this has
            // already been accounted for and shouldn't happen because of vertex check
            if (signedArea < 0) { // relative order of these two points is incorrect
            // (not ccw along pixel boundary)
                boundary.push_back(sideIntersections.at(1));
                boundary.push_back(sideIntersections.at(0));
            } else if (signedArea > 0) {
                boundary.push_back(sideIntersections.at(0));
                boundary.push_back(sideIntersections.at(1));
            }
        }
    }
    // check for null pointer
    if (polygon) {
        *polygon = boundary;
    }
    return shoelace(boundary);
}

double Pixel::doubleIntegral(function<double(double, double)> func, Triangle &t) {
    vector<Point> boundary;
    double area = intersectionArea(t, &boundary);
    if (area == 0) {
        return 0;
    }
    // approximate func by value at average x, average y
    double avgX, avgY;
    Point center(x,y);
    averageXY(&avgX, &avgY, boundary, center);
    return func(avgX, avgY) * area * color;
}