#include<iostream>
#include<string>
#include<cmath>
#include<cassert>
#include<vector>

#include "Point.hpp"
#include "Mesh.hpp"

template <typename MeshType>
typedef struct CollisionDetector {

	// used types ----------------------
	struct Collision;
	typedef MeshType::Triangle Triangle;
	typedef MeshType::Node Node;
	typedef MeshType::Edge Edge;
	typedef std::vector<MeshType>::iterator MeshIter;
	typedef std::vector<Collision>::iterator CollIter;

	// public interfaces ------------------------------------

	/** Constructor */
	CollisionDetector () {
	}

	/** main function
	 * takes in iterator over vector of meshes
	 * stores results in internal vector of colisions
	 * overwrites old collisions
	 * */
	void check_collisions(MeshIter begin, MeshIter end);

	/** accessors to found collisions
	 */
	CollIter begin();
	CollIter end();

	/** struct to store collision information
	 * @a n_ the node that is inside another mesh
	 * @a t_ the triangle in the mesh that is closest to the node
	 * 	(most likely guess for the collision)
	 */
	struct Collision {
		Node n_;
		Triangle t_;

		Collision(Node n, Triangle t) : n_(n), t_(t) {}
	} Collision;




	/** Determines whether a line segment passes through triangle
	 * @param[in] Triangle is defined by t1, t2, t3
	 * @param[in] Line segment is defined by a and b
	 * @returns true if the line segement intersects plane of triangle
	 */
	bool is_colliding(Point t1, Point t2, Point t3, Point p0, 
						 Point p1) {
		// Declare variables
		Point n;
		double t;

		// r(t) = p0 + t*(p1 - p0) is the equation of the line; find
		// t at the point of intersection; t between 0 and 1 indicates
		// that the intersection point lies on the line segment
		n = plane_normal(t1, t2, t3);
		t = (dot(n, t1) - dot(p0, n)) / (dot(n, p1 - p0));
		if( 0 <= t <= 1 ) {
			Point intersect = plane_line_intersect(t1, t2, t3, p0, p1);
			return is_inside_triangle(t1, t2, t3, intersect);
		}
		else {
			return false;
		}
	}

	/** Returns true if a given point is inside the triangle
	 * @param[in] Triangle defined by t1, t2, t3
	 * @param[in] Point to check is p
	 * @returns true is point is inside triangle, false otherwise
	 */
	bool is_inside_triangle(Point t1, Point t2, Point t3, Point p) {
		return on_same_side(t1, t2, t3, p) && 
			   on_same_side(t1, t3, t2, p) && 
			   on_same_side(t2, t3, t1, p);
	}

	/** Returns true if two points are on the same side of a line
	 * @param[in] Line is defined by points a and b
	 * @param[in] The two points are p0, and p1
	 * @returns true if p0 and p1 are on the same side of line defined
	 * 	by b - a
	 */
	bool on_same_side(Point a, Point b, Point p0, Point p1) {
		Point cross1 = cross(b - a, p0 - a);
		Point cross2 = cross(b - a, p1 - a);
		return dot(cross1, cross2) >= 0;
	}

	/** Compute the intersection of a plane and a line
	 * @pre p0 and p1 must form a valid line, p0 != p1 && 
	 * 	dot(p0, p1) != 0
	 * @pre line formed by p0 and p1 is not parallel to the plane
	 * @pre p0 and p1 must form a line not on the plane
	 * @param[in] Line is defined by the two points p0 and p1
	 * @param[in] Plane is defined by three points t1, t2, t3
	 * @return a Point of intersection between the point and plane
	 */
	Point plane_line_intersect(Point t1, Point t2, Point t3, 
	                           Point p0, Point p1) {
		// Declare variables
		Point n;
		double t;

		n = plane_normal(t1, t2, t3);
		// Check preconditions
		assert( !equal(p0, p1) );
		assert( !equal(dot(p0, p1), 0) );
		assert( !equal(dot(n, p1 - p0), 0) );

		// r(t) = p0 + t*(p1 - p0) is the equation of the line; find
		// t at the point of intersection; then compute the intersect
		t = (dot(n, t1) - dot(p0, n)) / (dot(n, p1 - p0));
		return p0 + t * (p1 - p0);
	}

	/** Compute the normal to the plane defined by the three
	 *	points p1, p2, and p3.
	 */
	Point plane_normal(Point p1, Point p2, Point p3) {
		Point v1 = p3 - p1;
		Point v2 = p2 - p1;
		return cross(v1, v2);
	}

	/* Return true if the point is on the plane defined by v1, v2, v3*/
	bool is_on_plane(Point v1, Point v2, Point v3, Point p) {
		// Check if the vectors joining the point p to every vertex
		// on the plane dotted with the normal vector is 0
		Point n = plane_normal(v1, v2, v3);
		return equal(dot(n, v1 - p), 0) && 
			   equal(dot(n, v2 - p), 0) &&
			   equal(dot(n, v3 - p), 0);
	}

	/* Return true if the point is on the line defined by a and b */
	bool is_on_line(Point a, Point b, Point p) {
		// Declare variables
		Point v1, v2, crpr;
		// If b and p jut out from a in the same or opposite 
		// directions then p lies on the line formed by a and b
		v1 = b - a;
		v2 = p - a;
		crpr = cross(v1, v2);
		return equal(crpr, Point(0, 0, 0));
	}

	/* Return the area of the triangle formed by t1, t2, t3 */
	double triangle_area(Point t1, Point t2, Point t3) {
		return 0.5 * norm(cross(t3 - t1, t3 - t2));
	}

	/* Returns true if the two doubles are within epsilon */
	bool equal(double a, double b, double epsilon=0.001) {
		return std::abs(a - b) < epsilon;
	}

	/* Returns true if the two points x, y, z coords satisfy equal */
	bool equal(Point a, Point b) {
		return equal(a.x, b.x) && equal(a.y, b.y) && equal(a.z, b.z);
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

} Collision;
