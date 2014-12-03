#include<iostream>
#include "Point.hpp"

typedef struct Collision {
	/** Constructor */
	Collision () {
	}

	/** Prints out a point to the console */
	void print_point(Point p) {
		std::cout << "(" << p.x << ", " << p.y << ", " << p.z << ")";
	}

	/** Runs a sequence of tests on the collision detector */
	void test() {
		std::cout << "This is a test" << std::endl;
		Point p = Point(1, 2, 3);
		std::cout<< "Here is a point: "; print_point(p);
		std::cout << std::endl;
	}
} Collision;
