#include "point.hpp"

Point::Point(double _x, double _y, bool border_x, bool border_y) {
    x = _x;
    y = _y;
    borderX = border_x;
    borderY = border_y;
}

double Point::getX() const {
    return x;
}

double Point::getY() const {
    return y;
}

bool Point::isBorderX() const {
    return borderX;
}

bool Point::isBorderY() const {
    return borderY;
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
    os << setprecision(10) << "(" << p.x << ", " << p.y << ")";
    return os;
}