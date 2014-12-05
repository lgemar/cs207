#include "CollisionDetector.hpp"

// function prototypes
void test_is_colliding();
void test_plane_line_intersect();
void test_is_inside_triangle();
void test_on_same_side();
void test_plane_normal();

void test_is_colliding() {
	db("Here will be some triangle-line segment tests");
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
	dbg("All assertions pass");
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
	db("Hit rate and hit probabilities");
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
	db("Hit rate and hit probabilities");
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
	db("Point 1: ",p1);
	db("Point 2: ",p2);
	db("Point 3: ",p3);


	// Compute the normal
	Point the_normal = plane_normal(p1, p2, p3);
	db("Plane normal: ",p1);


	// Check that all the dot products are 0
	assert( equal(dot(the_normal, p2 - p1), 0) );
	assert( equal(dot(the_normal, p3 - p1), 0) );
	assert( equal(dot(the_normal, p3 - p2), 0) );
	dbg("Passing all asserts");
}

/** Runs a sequence of tests on the collision detector */
void test_geometry() {
	db("=====> Testing plane normals <======");
	test_plane_normal();
	db("=====> Testing on same side functionality<=====");
	test_on_same_side();
	db("=====> Testing triangle interior checks <=====");
	test_is_inside_triangle();
	db("=====> Testing plane-line intersection <=====");
	test_plane_line_intersect();
	db("=====> Testing triangle-segment collision <=====");
	test_is_colliding();
}

void test_add_remove() {
	db("======> Testing Add Remove <=============");
	db("creating collider");
	typedef Mesh<char, char, char> MeshType;
	typedef CollisionDetector<MeshType> collider;
	typedef collider::Tag Tag;
	collider c = collider();

	db("creating a few meshes");
	MeshType m1;
	MeshType m2;
	MeshType m3;
	MeshType m4;
	MeshType m5;

	db("adding a node to each");
	m1.add_node(Point());
	m2.add_node(Point());
	m3.add_node(Point());
	m4.add_node(Point());
	m5.add_node(Point());

	db("adding objects");
	c.add_object(m1);
	c.add_object(m2);
	c.add_object(m3);
	c.add_object(m4);
	c.add_object(m5);

	db("drawing graph");
	db("should see 5 fully connected nodes........");
	c.print_graph();

	db("removing objects");
	c.remove_object(m4);
	c.remove_object(m1);
	c.remove_object(m2);
	c.remove_object(m3);
	c.remove_object(m5);

	dbg("No errors!");
}

void test_tags() {
	db("======> Testing tags <=============");
		db("creating collider");
		typedef Mesh<char, char, char> MeshType;
		typedef CollisionDetector<MeshType> collider;
		typedef collider::Tag Tag;
		collider c = collider();

		db("creating a few meshes");
		MeshType m0;
		MeshType m1;
		MeshType m2;
		MeshType m3;
		MeshType m4;
		MeshType m5;
		MeshType m6;
		MeshType m7;
		MeshType m8;

		db("adding a node to each");
		m0.add_node(Point());
		m1.add_node(Point());
		m2.add_node(Point());
		m3.add_node(Point());
		m4.add_node(Point());
		m5.add_node(Point());
		m6.add_node(Point());
		m7.add_node(Point());
		m8.add_node(Point());

		db("creating some tags");
		Tag t0 = Tag(); // default tag
		Tag t1 = c.getNoneTag(); // checks against nothing
		Tag t2 = c.getOtherTag(); // checks against not self
		Tag t3 = c.getSelfTag(); // checks against self
		Tag t4 = c.get_tag(true); // only on t5
		Tag t5 = c.get_tag(true); // only on t4

		t4.add(t5);
		t4.add(t2);
		t4.add(t3);
		t5.add(t4);
		t5.add(t0);
		t5.add(t1);

		db("printing tags");
		db("t0 id:", t0.id_);
		db("t1 id:", t1.id_);
		db("t2 id:", t2.id_);
		db("t3 id:", t3.id_);
		db("t4 id:", t4.id_);
		db("t5 id:", t5.id_);
		db("");

		db("printing tags");
		db("t0 w:", t0.white_);
		db("t1 w:", t1.white_);
		db("t2 w:", t2.white_);
		db("t3 w:", t3.white_);
		db("t4 w:", t4.white_);
		db("t5 w:", t5.white_);

		db("printing lists");
		db_vec("t0 l:", t0.list_);
		db_vec("t1 l:", t1.list_);
		db_vec("t2 l:", t2.list_);
		db_vec("t3 l:", t3.list_);
		db_vec("t4 l:", t4.list_);
		db_vec("t5 l:", t5.list_);

		db("adding objects");
		c.add_object(m0); // connect to all
		c.add_object(m1,t0); // connect to all
		c.add_object(m2,t1); // connects to nothing
		c.add_object(m3,t4); // checks against not t2
		c.add_object(m4,t5); // checks against not t2
		c.add_object(m5,t4); // checks against only t3
		c.add_object(m6,t5); // checks against only t3
		c.add_object(m7,t4); // checks against only t5
		c.add_object(m8,t5); // checks aginst only t4

		db("printing lists");
		db_vec("t0 l:", t0.list_);
		db_vec("t1 l:", t1.list_);
		db_vec("t2 l:", t2.list_);
		db_vec("t3 l:", t3.list_);
		db_vec("t4 l:", t4.list_);
		db_vec("t5 l:", t5.list_);

		db("printing graph");
		c.print_graph();
		db("look at test_collision:228 for reference");
		dbg("No errors!");
}

int main () {

	//test_geometry();
	//db("");
	//test_add_remove();
	//db("");
	test_tags();
	db("");


	return 0;
}

