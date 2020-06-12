#include "point.hpp"

Point::Point(double _x, double _y) {
    x = _x;
    y = _y;
}

double Point::getX() {
    return x;
}

double Point::getY() {
    return y;
}

double Point::distance(Point &other) {
    double deltaX = x - other.getX();
    double deltaY = y - other.getY();
    return pow(deltaX * deltaX + deltaY * deltaY, 0.5);
}

void Point::move(double deltaX, double deltaY) {
    x += deltaX;
    y += deltaY;
}

bool Point::operator==(const Point& other) const {
    return x == other.x && y == other.y;
}

bool Point::operator!=(const Point &other) const {
    return !(*this == other);
}

ostream& operator<<(ostream& os, const Point &p) {
    os << "(" << p.x << ", " << p.y << ")";
    return os;
}