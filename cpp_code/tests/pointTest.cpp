#include <gtest/gtest.h>
#include "../src/point.hpp"

using namespace std;

TEST(PointTest, Move) {
    Point a(0,0);
    a.move(3, 4.5);
    ASSERT_EQ(a, a);
    ASSERT_EQ(3, a.getX());
    ASSERT_EQ(4.5, a.getY());
}

TEST(PointTest, Equality) {
    Point a(0,0);
    Point b(0,0);
    ASSERT_EQ(a, b);
    a.move(1,1);
    ASSERT_NE(a, b);
    b.move(1,1);
    ASSERT_EQ(a, b);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
