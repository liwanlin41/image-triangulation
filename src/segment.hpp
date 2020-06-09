#ifndef segment_h
#define segment_h

#include "point.hpp"

/**
 * represent a segment (point1, point2) in the plane
 */

class Segment {
    private:
        Point *endpoint1, *endpoint2;
    public:
        Segment(Point *a, Point *b);
        double length();

        // return true if this segment intersects other
        bool intersects(Segment other);

        // return intersection point if this intersects other
        // undefined behavior if intersection does not exist or if segments overlap with positive length
        Point getIntersection(Segment other);
};

#endif