#include "CollisionDetector.hpp"
#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"

typedef Mesh<char, char, char> MeshType;
typedef CollisionDetector<MeshType> collider;
typedef typename collider::Collision collision;
typedef typename MeshType::node_type Node;
typedef collider::Tag Tag;

// Functor for find function
struct InCollision {
	MeshType* m_;
	InCollision(MeshType* m) : m_(m) {
	}
	bool operator()(collision c) {
		return m_ == c.mesh1 || m_ == c.mesh2;
	}
};

// generate a mesh ball around the point
void add_ball(MeshType& m, Point p, double r = .3) {

	size_t slices = 6;
	double dth = 2*3.14 / slices;
	std::vector<Node> tops;
	std::vector<Node> bots;

	tops.reserve(100);
	bots.reserve(100);

	// creating points for geodesic dome
	Node tip1 = m.add_node(p + Point(0,0,r));
	Node tip2 = m.add_node(p + Point(0,0,-r));
	double tz = .45*r;
	double tr = .8*r;
	for (size_t i = 0; i < slices; ++i) {
		double th = i*dth;
		tops.push_back(m.add_node(p + Point(tr*cos(th), tr*sin(th), tz)));
	}
	for (size_t i = 0; i < slices; ++i) {
		double th = dth/2 + i*dth;
		bots.push_back(m.add_node(p + Point(tr*cos(th), tr*sin(th), -tz)));
	}

	assert(tops.size() == slices);
	assert(bots.size() == slices);

	// creating triangles
	for (size_t i = 0; i < slices; ++i) {
		int j = (i + 1) % slices;
		// top and bottom
		m.add_triangle(tip1, tops[i], tops[j]);
		m.add_triangle(tip2, bots[i], bots[j]);

		// middle
		m.add_triangle(tops[i], bots[i], tops[j]);
		m.add_triangle(bots[i], tops[j], bots[j]);
	}
}

// moves all points in a mesh
void move_mesh(MeshType& m, Point v, double dt) {
	v = v * dt;
	for (auto it = m.vertex_begin(); it != m.vertex_end(); ++it) {
		(*it).position() = (*it).position() + v;
	}
}

double drand(double lim) {
	return (rand() % 10000) * lim / 10000;
}

int main () {

	// initialization
	srand(time(NULL));
	int N = 5;

	// create meshes
	std::vector<MeshType> meshes;

	for (int i = 0; i < N; ++i)
		meshes.push_back(MeshType());

	// create checker
	collider c;
	double pos = 0;
	for (auto it = meshes.begin(); it != meshes.end(); ++it) {
		add_ball(*it, Point(pos++, 0, 0));
		c.add_object(*it);
	}

	// debug initial state
	c.print_graph();

	// Launch the SDLViewer
	CS207::SDLViewer viewer;
	viewer.launch();

	// Add the vertices of the mesh to the viewer
	auto vertex_map = viewer.empty_vertex_map(meshes[0]);
	for (auto it = meshes.begin(); it != meshes.end(); ++it) {
		viewer.add_nodes((*it).vertex_begin(), (*it).vertex_end(),
				   CS207::GreenColor(), NodeToPoint(), vertex_map);
		viewer.add_edges((*it).edge_begin(), (*it).edge_end(), vertex_map);
	}
	viewer.center_view();
	// Display viewer
	double dt = .003;
	//move_mesh(meshes[0], Point(drand(2)-1,drand(2)-1,drand(2)-1), dt);

	db("nodes:", meshes[0].num_nodes());
	db("edges:", meshes[0].num_edges());
	db("tris:", meshes[0].num_triangles());

	for (double t = 0; t < 10; t += dt) {

		// moving
		for (auto it = ++meshes.begin(); it != meshes.end(); ++it)
			move_mesh(*it, Point(-1,drand(2)-1,drand(2)-1), dt);

		// moving other
		move_mesh(*(++meshes.begin()), Point(2,0,0), dt);

	    // check for collision
	    c.check_collisions();
	    size_t count = std::distance(c.begin(), c.end());

	    // redraw

	    for (int i = 0; i < N; ++i) {
			//If there is a collision display those meshes in red
			if( std::find_if(c.begin(),c.end(),InCollision(&meshes[i]))!= c.end()) {
				viewer.add_nodes(meshes[i].vertex_begin(),
				meshes[i].vertex_end(),
				CS207::RedColor(), NodeToPoint(),
				vertex_map);
			}
			else {
				viewer.add_nodes(meshes[i].vertex_begin(), 
					meshes[i].vertex_end(),
					CS207::GreenColor(), NodeToPoint(),
					vertex_map);
			}
		}
	    viewer.set_label(count);
	}
	return 0;
}

