#include<iostream>
#include<string>

#include "Point.hpp"

typedef struct Collision {
	/** Constructor */
	Collision () {
	}

	/** Runs a sequence of tests on the collision detector */
	void test() {
		print_line("=====> Testing print statements <======");
		test_prints();
	}

	/** Test the printing functionality */
	void test_prints() {
		Point p = Point(1, 2, 3);

		print_line("This is a line");
		print_string("Here is a point: ");
		print_point(p);
		end_line();
	}
	
	/** Prints out a point to the console */
	void print_point(Point p) {
		std::cout << "(" << p.x << ", " << p.y << ", " << p.z << ")";
	}

	/** Print a string to the console without starting a new line */
	void print_string(std::string s) {
		std::cout << s;
	}

	/** Print a string to the console and end the line */
	void print_line(std::string s) {
		std::cout << s << std::endl;
	}

	/** Prints a new line to the output stream */
	void end_line() {
		std::cout << std::endl;
	}

} Collision;
