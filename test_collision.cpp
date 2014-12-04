#include "CollisionDetector.hpp"
#include "Point.hpp"


// function prototypes
void test_is_colliding();
void test_plane_line_intersect();
void test_is_inside_triangle();
void test_on_same_side();
void test_plane_normal();
void test_print_statements();
void print_point(Point);
template <typename T>
void print_text(T);
void print_line(std::string);
void end_line();

void test_is_colliding() {
	print_line("Here will be some triangle-line segment tests");
}

void test_plane_line_intersect() {
	// Seed the random number generator
	srand(time(NULL));
	// Declare variables;
	Point v1, v2, v3, a, b, intersect;
	int sz = 10;
	int num_points = 100;
	int i = num_points;
	// Assert that each intersect point is on the line and plane
	while( i ) {
		// Randomize all the points
		v1 = Point(rand() % sz, rand() % sz, rand() % sz);
		v2 = Point(rand() % sz, rand() % sz, rand() % sz);
		v3 = Point(rand() % sz, rand() % sz, rand() % sz);
		a = Point(rand() % sz, rand() % sz, rand() % sz);
		b = Point(rand() % sz, rand() % sz, rand() % sz);
		// Make sure that a and b form a valid line
		if( equal(a, b) || equal(dot(a, b), 0) )
			continue;
		// Make sure that v1, v2, v3 form a valid plane
		if( equal(triangle_area(v1, v2, v3), 0) )
			continue;
		// Make sure that line is not parallel to the plane
		if( equal(dot(plane_normal(v1, v2, v3), a-b), 0) )
			continue;
		// Compute the intersection
		intersect = plane_line_intersect(v1, v2, v3, a, b);
		assert( is_on_plane(v1, v2, v3, intersect) );
		assert( is_on_line(a, b, intersect) );
		// Decrement the counter
		--i;
	}
	print_line("All assertions pass");
}

void test_is_inside_triangle() {
	// Seed the random number generator
	srand(time(NULL));
	// Declare variables
	Point t1, t2, t3, random_point;
	double hit_rate, hit_prob;
	int sz = 100;
	int hits = 0;
	int num_points = 10000;
	// Define a random triangle on the 10x10 grid on the x-y plane
	t1 = Point(rand() % sz, rand() % sz, 0);
	t2 = Point(rand() % sz, rand() % sz, 0);
	t3 = Point(rand() % sz, rand() % sz, 0);
	// Create an array of points in a szxsz grid on x-y plane
	int i = num_points;
	while( i ) {
		random_point = Point(rand() % sz, rand() % sz, 0);
		if( is_inside_triangle(t1, t2, t3, random_point) )
			++hits;
		--i;
	}
	// Print out the hit rate and the hit probability
	print_line("Hit rate and hit probabilities");
	hit_rate = (double) hits / num_points;
	hit_prob = triangle_area(t1, t2, t3) / (sz * sz);
	std::cout << "Hit rate: " << hit_rate << std::endl;
	std::cout << "Hit probability: " << hit_prob << std::endl;
}

void test_on_same_side() {
	// Seed the random number generator
	srand(time(NULL));
	// Declare variables
	Point a, b, p0, p1;
	double hit_rate, hit_prob;
	int sz = 100;
	int hits = 0;
	int num_points = 10000;
	// Define a line that cuts szXsz grid on the x-y plane in half
	a = Point(0, 0, 0);
	b = Point(sz, sz, 0);
	// Create an array of points in a szxsz grid on x-y plane
	int i = num_points;
	while( i ) {
		p0 = Point(rand() % sz, rand() % sz, 0);
		p1 = Point(rand() % sz, rand() % sz, 0);
		if( on_same_side(a, b, p0, p1) )
			++hits;
		--i;
	}
	// Print out the hit rate and the hit probability
	print_line("Hit rate and hit probabilities");
	hit_rate = (double) hits / num_points;
	hit_prob = 0.5;
	std::cout << "Hit rate: " << hit_rate << std::endl;
	std::cout << "Hit probability: " << hit_prob << std::endl;
}

void test_plane_normal() {
	Point p1 = Point(rand() % 100, rand() % 100, rand() % 100);
	Point p2 = Point(rand() % 100, rand() % 100, rand() % 100);
	Point p3 = Point(rand() % 100, rand() % 100, rand() % 100);

	// Print out the three starting points
	print_text("Point 1: ");
	print_point(p1);
	end_line();
	print_text("Point 2: ");
	print_point(p2);
	end_line();
	print_text("Point 3: ");
	print_point(p3);
	end_line();

	// Compute the normal
	Point the_normal = plane_normal(p1, p2, p3);
	print_text("Plane normal: ");
	print_point(p1);
	end_line();

	// Check that all the dot products are 0
	assert( equal(dot(the_normal, p2 - p1), 0) );
	assert( equal(dot(the_normal, p3 - p1), 0) );
	assert( equal(dot(the_normal, p3 - p2), 0) );
	print_line("Passing all asserts");
}

/** Test the printing functionality */
void test_print_statements() {
	Point p = Point(1, 2, 3);

	print_line("This is a line");
	print_text("Here is a point: ");
	print_point(p);
	end_line();
}

/** Prints out a point to the console */
void print_point(Point p) {
	std::cout << "(" << p.x << ", " << p.y << ", " << p.z << ")";
}

/** Print a string to the console without starting a new line */
template<typename T>
void print_text(T s) {
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

/** Runs a sequence of tests on the collision detector */
void test() {
	print_line("=====> Testing print statements <======");
	test_print_statements();
	print_line("=====> Testing plane normals <======");
	test_plane_normal();
	print_line("=====> Testing on same side functionality<=====");
	test_on_same_side();
	print_line("=====> Testing triangle interior checks <=====");
	test_is_inside_triangle();
	print_line("=====> Testing plane-line intersection <=====");
	test_plane_line_intersect();
	print_line("=====> Testing triangle-segment collision <=====");
	test_is_colliding();
}

int main () {
	typedef Mesh<int, int, int> MeshType;
	typedef CollisionDetector<MeshType> collider;
	collider c = collider();
	test();
	return 0;
}

