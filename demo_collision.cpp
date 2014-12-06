#include "CollisionDetector.hpp"
#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"

typedef Mesh<char, char, char> MeshType;
typedef CollisionDetector<MeshType> collider;
typedef typename MeshType::node_type Node;
typedef collider::Tag Tag;

// generate a mesh ball around the point
void add_ball(MeshType& m, Point p, double r = .3) {

	size_t slices = 10;
	std::vector<Node> tops;
	std::vector<Node> bots;

	// creating points for geodesic dome
	Node tip1 = m.add_node(p + Point(0,0,r));
	Node tip2 = m.add_node(p + Point(0,0,-r));
	double tz = .45*r;
	double tr = .8*r;
	for (double th = 0; th < 2 * 3.14; th += (2*3.14 / slices)) {
		tops.push_back(m.add_node(p + Point(tr*cos(th), tr*sin(th), tz)));
	}
	for (double th = 3.14 / slices; th < 2 * 3.14; th += (2*3.14 / slices)) {
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

int main () {

	// initialization
	srand(time(NULL));
	int N = 3;

	// create meshes
	std::vector<MeshType> meshes;
	for (int i = 0; i < N; ++i)
		meshes.push_back(MeshType());




	// create checker
	collider c;
	int pos = 0;
	for (auto it = meshes.begin(); it != meshes.end(); ++it) {
		add_ball(*it, Point(pos,0,0));
		c.add_object(*it);
		++pos;
	}

	// debug initial state
	c.print_graph();

	// Launch the SDLViewer
	CS207::SDLViewer viewer;
	viewer.launch();

	// Add the vertices of the mesh to the viewer
	auto vertex_map = viewer.empty_vertex_map(meshes[0]);
	for (int i = 0; i < N; ++i) {
		viewer.add_nodes(meshes[i].vertex_begin(), meshes[i].vertex_end(),
				   CS207::GreenColor(), NodeToPoint(), vertex_map);
		viewer.add_edges(meshes[i].edge_begin(), meshes[i].edge_end(), vertex_map);
	}
	// Display viewer
	viewer.center_view();

	double dt = .002;
	for (double t = 0; t < 10; t += dt) {

		// update
		for (int i = 1; i < N; ++i)
			move_mesh(meshes[i], Point(-1,0,0), dt);

	    // check for collision
	    c.check_collisions();
	    size_t count = std::distance(c.begin(), c.end());

	    // redraw
	    for (int i = 0; i < N; ++i)
	    	viewer.add_nodes(meshes[i].vertex_begin(), meshes[i].vertex_end(),
					   CS207::GreenColor(), NodeToPoint(), vertex_map);

	    viewer.set_label(count);

	  }
	return 0;
}

