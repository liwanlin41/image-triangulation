#include "point.cuh"

Point::Point() {}

Point::Point(double x_, double y_, bool border_x, bool border_y):
	x(x_), y(y_), borderX(border_x), borderY(border_y) {}

double Point::getX() const {
	return x;
}

double Point::getY() const {
	return y;
}

__device__ bool Point::isBorderX() const {
	return borderX;
}

__device__ bool Point::isBorderY() const {
	return borderY;
}

double Point::distance(Point &other) {
	double deltaX = x - other.x;
	double deltaY = y - other.y;
	return pow(deltaX * deltaX + deltaY * deltaY, 0.5);
}

__device__ void Point::move(double deltaX, double deltaY) {
	x += deltaX;
	y += deltaY;
}

__device__ bool Point::operator==(const Point& other) const {
	return x == other.x && y == other.y;
}

__device__ bool Point::operator!=(const Point &other) const {
	return !(*this == other);
}

ostream& operator<<(ostream& os, const Point &p) {
	os << setprecision(10) << "(" << p.x << ", " << p.y << ")";
	return os;
}