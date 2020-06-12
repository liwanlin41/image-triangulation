#ifndef segment_h
#define segment_h

#include "point.hpp"
#include "matrix.hpp"

/**
 * represent a segment (point1, point2) in the plane
 */

class Segment {
    private:
        Point *endpoint1, *endpoint2;
    public:
        Segment(Point *a, Point *b);
        double length();
        // return endpoint1
        Point getStart();
        // return endpoint2
        Point getEnd();

        // return true if this segment intersects other in exactly one point
        bool intersects(Segment other);

        // return intersection point if this intersects other
        // undefined behavior if intersection does not exist or if segments overlap with positive length
        Point getIntersection(Segment other);

        // return the 2x1 unit normal to this segment which is on the right when going from endpoint1 to endpoint2
        Matrix unitNormal();
        // return the 2x1 normal to this segment that has length |segment|/2
        Matrix scaledNormal();

        friend ostream& operator<<(ostream& os, const Segment &seg);
};

#endif