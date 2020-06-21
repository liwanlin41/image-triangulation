#include "pixel.cuh"

// helper functions

// determine whether two points are "essentially" equal (floating point error)
__device__ bool approxEqual(Point &a, Point &b, double tolerance = 1e-12) {
	return (a.distance(b) < tolerance);
}

// compute (unsigned) area of the polygon enclosed by points,
// where edegs of the polygon are given by points[i] -- points[i+1]
__device__ double shoelace(Point *points, int &size) {
	if (size < 3) {
		return 0;
	}
	double area = 0;
	for(int i = 0; i < size; i++) {
		double x0 = points[i].getX();
		double y0 = points[i].getY();
		double x1 = points[(i+1)%size].getX();
		double y1 = points[(i+1)%size].getY();
		area += (x0 * y1 - x1 * y0);
	}
	// in practice points is supposed to be ccw
	// up to floating point errors that don't affect area
	assert(area >= 0);
	return area/2;
}

// compute integral of x over polygon points and store it in totalX, sim for y
// center is a reference point inside the pixel; even if it lies outside the polygon,
// using signed areas means the result will still be correct
__device__ void integrateXY(double *totalX, double *totalY, Point *points, int &size, Point &center) {
	double sumX = 0;
	double sumY = 0;
	for(int i = 0; i < size; i++) {
		// average value over a triangle is just the centroid
		double centroidX = (points[i].getX() + points[(i+1)%size].getX() + center.getX())/3;
		double centroidY = (points[i].getY() + points[(i+1)%size].getY() + center.getY())/3;
		double triangleArea = Triangle::getSignedArea(&center, &points[i], &points[(i+1)%size]);
		// weight the average
		sumX += centroidX * triangleArea;
		sumY += centroidY * triangleArea;
	}
	*totalX = sumX;
	*totalY = sumY;
}

// compute average values of x, y over the polygon enclosed by points
// and put them in the given variables
// center is again a reference point
__device__ void averageXY(double *avgX, double *avgY, Point *points, int &size, Point &center) {
	double totalX;
	double totalY;
	integrateXY(&totalX, &totalY, points, size, center);
	double totalArea = shoelace(points, size);
	*avgX = totalX / totalArea;
	*avgY = totalY / totalArea;
}

Pixel::Pixel(int x_, int y_, int c) : x(x_), y(y_), color(c) {
	corners[0] = Point(x-0.5, y-0.5);
	corners[1] = Point(x+0.5, y-0.5);
	corners[2] = Point(x+0.5, y+0.5);
	corners[3] = Point(x-0.5, y+0.5);
}

int Pixel::getColor() {
	return color;
}

__device__ bool Pixel::containsPoint(Point &p) {
	double px = p.getX();
	double py = p.getY();
	return (-0.5+x <= px && px <= 0.5+x) && (-0.5+y <= py && py <= 0.5+y);
}

__device__ double Pixel::intersectionLength(Segment &e, double *xVal, double *yVal) {
	Point intersections[2]; // hold intersections
	int numPts; // track number of intersection points detected thus far
	for(int i = 0; i < 4; i++) {
		// retrieve a side of the pixel; at most two will have an 
		// intersection unless intersection is at corners
		Segment side(&corners[i], &corners[(i+1)%4]);
		if (side.intersects(e)) {
			Point intersectionPoint = side.getIntersection(e);
			bool isNewPoint = true; // whether this intersection is a new distinct point
			for(Point &pt : intersections) {
				if (approxEqual(pt, intersectionPoint)) {
					isNewPoint = false;
				}
			}
			if (isNewPoint) {
				intersections[numPts] = intersectionPoint;
				numPts++;
			}
		}
	}
	// handle segment endpoints potentially inside the pixel
	if (numPts < 2) {
		Point start = *(e.endpoint1);
		Point end = *(e.endpoint2);
		if (containsPoint(start)) {
			intersections[numPts] = start;
			numPts++;
		}
		if (containsPoint(end)) {
			intersections[numPts] = end;
			numPts++;
		}
	}
	if (numPts < 2) {
		return 0;
	}
	Segment contained(&intersections[0], &intersections[1]);
	// check for null pointers, assign midpoint coords
	if (xVal && yVal) {
		*xVal = (intersections[0].getX() + intersections[1].getX())/2;
		*yVal = (intersections[0].getY() + intersections[1].getY())/2;
	}
	return contained.length();
}

/*
__device__ double Pixel::lineIntegral(nvstd::function<double(double, double)> &func, Segment &e) {
	double midX, midY;
	double length = intersectionLength(e, &midX, &midY);
	if (length == 0) {
		return 0;
	}
	return func(midX, midY) * length * color;	
}

__device__ double Pixel::lineIntegral(nvstd::function<double(double, double)> &func, Triangle &t) {
	double total = 0;
	for(int i = 0; i < 3; i++) {
		Segment seg(t.vertices[i], t.vertices[(i+1)%3]); // preserve orientation
		total += lineIntegral(func, seg);
	}
	return total;
}
*/

__device__ double Pixel::intersectionArea(Triangle t, Point* polygon, int *size) {
	int numPoints = 0; // track number of points in polygon
	Point boundary[16]; // there should only be max 10 points on the boundary,
	// but allow some room for floating point error
	int inInd; // index of some triangle vertex that lies inside pixel (may not exist)
	Segment triangleSides[3]; // hold sides of triangle

	// goal: compute boundary of the intersection

	for(int i = 0; i < 3; i++) {
		triangleSides[i] = Segment(t.vertices[0], t.vertices[(i+1)%3]);
        // add triangle vertices which may be inside the pixel, but don't add corners
        bool isCorner = false;
        for(int j = 0; j < 4; j++) {
            if (*(t.vertices[i]) == corners[j]) isCorner = true;
        }
        if (!isCorner && containsPoint(*(t.vertices[i]))) {
            inInd = i;
			boundary[numPoints] = *(t.vertices[i]);
			numPoints++;
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
        bool safelyOriented = (numPoints != 1) || 
            (Triangle::getSignedArea(&corners[i], t.vertices[(inInd+1)%3], t.vertices[(inInd+2)%3]) >= 0);
        if (!t.contains(corners[i]) && safelyOriented) {
            start = i;
        }
    }
    for(int i = 0; i < 4; i++) {
        // first determine if corner of pixel is inside
        Point corner = corners[(i+start) % 4];
        Segment side(&corner, &corners[(i+start+1)%4]);
        if (t.contains(corner)) {
			boundary[numPoints] = corner;
			numPoints++;
        }
        // determine intersections with side (i, i+1)
		Point sideIntersections[2];
		int intersectNum; // track index in sideIntersections
        for(Segment e : triangleSides) {
            if (side.intersects(e)) {
                Point intersectionPoint = side.getIntersection(e);
                // check to see if this point is already accounted for by corners
                // or by triangle vertices; if it isn't exactly equal it won't contribute to area
                // (and the lack of exact equality is likely due to floating point error)
                if (!approxEqual(intersectionPoint, corner) && !approxEqual(intersectionPoint, corners[(i+start+1)%4])) {
                    bool isVertex = false;
                    for(Point *tVertex : t.vertices) {
                        if (approxEqual(intersectionPoint, *tVertex)) {
                            isVertex = true;
                        }
                    }
                    if (!isVertex) {
						sideIntersections[intersectNum] = intersectionPoint;
						intersectNum++;
                    }
                }
            }
        }
        // note a triangle can intersect a given side at most twice
        assert(intersectNum <= 2);
        // handle normal case where there is only one intersection with this side
        if (intersectNum == 1) {
			boundary[numPoints] = sideIntersections[0];
			numPoints++;
        } else if (intersectNum == 2) {
            Point center(this->x, this->y); // center of this pixel
            double signedArea = Triangle::getSignedArea(&center, &sideIntersections[0], &sideIntersections[1]);
            // if signedArea == 0, sideIntersections must contain two of the same point
            // which means one vertex of the triangle is on the side; this has
            // already been accounted for and shouldn't happen because of vertex check
            if (signedArea < 0) { // relative order of these two points is incorrect
            // (not ccw along pixel boundary)
				boundary[numPoints] = sideIntersections[1];
				boundary[numPoints + 1] = sideIntersections[0];
				numPoints += 2;
            } else if (signedArea > 0) {
				boundary[numPoints] = sideIntersections[0];
				boundary[numPoints + 1] = sideIntersections[1];
				numPoints += 2;
            }
        }
    }
    // check for null pointer
    if (polygon && size) {
        polygon = boundary;
		*size = numPoints;
    }
    return shoelace(boundary, numPoints);
}

/*
__device__ double Pixel::doubleIntegral(nvstd::function<double(double, double)> &func, Triangle &t) {
	Point boundary[16]; // again, leave ample space
	int size;
	double area = intersectionArea(t, boundary, &size);
	if (area == 0) {
		return 0;
	}
	// approximate func by value at average x, average y
	double avgX, avgY;
	Point center(x, y);
	averageXY(&avgX, &avgY, boundary, size, center);
	return func(avgX, avgY) * area * color;
}
*/