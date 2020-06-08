#include <iostream>

#include "point.hpp"
#include "segment.hpp"
#include "matrix.hpp"
#include "triangle.hpp"

// eventually input will be the path to an image file?
int main(int argc, char* argv[]) {
	Point a(50, 100);
	std::cout << a.getX() << ", " << a.getY() << std::endl;
	return 0;
}
