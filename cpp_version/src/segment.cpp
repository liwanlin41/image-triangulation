class Point {
    float x;
    float y;
    Point(float _x, float _y) {
        x = _x;
        y = _y;
    }
};

class Segment {
    Point endpoint1;
    Point endpoint2;
    Segment(Point &_endpoint1, Point &_endpoint2) {
        endpoint1 = _endpoint1;
        endpoint2 = _endpoint2;
    }
};