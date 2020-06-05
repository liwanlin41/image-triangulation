/**
 * generic class representing a point (x,y) on the plane
 * mutable; allows perturbing vertices
 */

template<typename T>
class Point {
    private:
        T x, y;
    public:
        Point(T _x, T _y) {
            x = _x;
            y = _y;
        }

        T getX() {
            return x;
        }

        T getY() {
            return y;
        }
        
        void move(T delta_x, T delta_y) {
            x += delta_x;
            y += delta_y;
        }
};