#include <iostream>
#include <assert.h>
#include "../src/point.cpp"

using namespace std;

int main() {
    // simple test to make sure move behaves as expected
    Point<double> a(0,0);
    a.move(3, 4.5);
    assert (a.getX() == 3);
    assert (a.getY() == 4.5);

    cout << "all tests passed" << endl;
    return 0;
}