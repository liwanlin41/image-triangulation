/**
 * represent a point (x,y) on the plane
 * mutable; allows perturbing vertices
 */

class Point {
    private:
        double x, y;
    public:
        Point(double _x, double _y);
        double getX();
        double getY();
        void move(double delta_x, double delta_y);
};