#include "../headers/point.hpp"

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