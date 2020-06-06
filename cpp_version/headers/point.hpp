#ifndef point_h
#define point_h
/**
 * represent a point (x,y) on the plane
 * mutable; allows perturbing vertices
 */

class Point {
    private:
        double x, y;
    public:
        Point(double x, double y);
        double getX();
        double getY();
        void move(double deltaX, double deltaY);
        bool operator==(const Point &other) const;
};
#endif