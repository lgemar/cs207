/**
 * @file shallow_water.cpp
 * Implementation of a shallow water system using Mesh
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D point list (one per line) defined by three doubles
 * Second file: Triangles (one per line) defined by 3 indices into the point list
 */

#include <fstream>
#include <cmath>
#include <math.h>
#include <string.h>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"

#include "Point.hpp"
#include "Mesh.hpp"

/* Prints out the variable s */
template <typename in, typename in2>
void db(in s, in2 s2) {
	std::cout << s << s2 << std::endl;
}

/* Prints out the variable s */
template <typename in>
void db(in s) {
	std::cout << s << std::endl;
}

/* Prints out the point p */
void db(Point p) {
	(void) p;
	std::cout << "(" << p.x << ", " << p.y << ", " << p.z << ")" << std::endl;
}

/* Prints out the point p */
template <typename in>
void db(in s, Point p) {
	std::cout << s << "(" << p.x << ", " << p.y << ", " << p.z << ")" << std::endl;
}

/** Water column characteristics */
struct QVar {
  double h;   // Height of column
  double hx;  // Height times average x velocity of column
  double hy;  // Height times average y velocity of column

  /** Default constructor.
   * A default water column is 1 unit high with no velocity. 
   */
  QVar() : h(1), hx(0), hy(0) {}
  /** Construct the given water column. */
  QVar(double h_, double hx_, double hy_) : h(h_), hx(hx_), hy(hy_) { }
  /* Implement QVar operators */
  QVar operator+(QVar o) const { return QVar(h + o.h, hx + o.hx, hy + o.hy); }
  QVar operator-(QVar o) const { return QVar(h - o.h, hx - o.hx, hy - o.hy); }
  void operator+=(QVar o) { *this = (*this + o); }
  void operator-=(QVar o) { *this = (*this - o); }
  QVar operator*(double s) const { return QVar(h*s, hx*s, hy*s); }
  QVar operator/(double s) const { return *this * 1/s; }
};

// Standard gravity (average gravity at Earth's surface) in meters/sec^2
static constexpr double grav = 9.80665;

// Define NodeData, EdgeData, TriData, etc
typedef struct my_triangle_data { double area_; QVar qvar_; } my_triangle_data;
typedef struct my_edge_data { 
	Point normal_; 
	// my_edge_data() : normal_(Point(0, 0, 0)) {}
} my_edge_data;
typedef struct my_vertex_data { double h; } my_vertex_data;

// Mesh types
typedef Mesh<my_vertex_data, my_edge_data, my_triangle_data> MeshType;
typedef MeshType::Link Link;
typedef MeshType::Triangle Triangle;
typedef MeshType::Vertex Vertex;
typedef MeshType::Edge Edge;
typedef MeshType::triangle_iterator triangle_iterator;
typedef std::set<Triangle> triangle_set;


/** Function object for calculating shallow-water flux.
 *          |n
 *   T_k    |---> n = (nx,ny)   T_m
 *   QBar_k |                   QBar_m
 *          |
 * @param[in] nx, ny Defines the 2D outward normal vector n = (@a nx, @a ny)
 *            from triangle T_k to triangle T_m. The length of n is equal to the
 *            the length of the edge, |n| = |e|.
 * @param[in] dt The time step taken by the simulation. Used to compute the
 *               Lax-Wendroff dissipation term.
 * @param[in] qk The values of the conserved variables on the left of the edge.
 * @param[in] qm The values of the conserved variables on the right of the edge.
 * @return The flux of the conserved values across the edge e
 */
struct EdgeFluxCalculator {
  QVar operator()(double nx, double ny, double dt,
                  const QVar& qk, const QVar& qm) {
    // Normalize the (nx,ny) vector
    double n_length = std::sqrt(nx*nx + ny*ny);
    nx /= n_length;
    ny /= n_length;

    // The velocities normal to the edge
    double wm = (qm.hx*nx + qm.hy*ny) / qm.h;
    double wk = (qk.hx*nx + qk.hy*ny) / qk.h;

    // Lax-Wendroff local dissipation coefficient
    double vm = sqrt(grav*qm.h) + sqrt(qm.hx*qm.hx + qm.hy*qm.hy) / qm.h;
    double vk = sqrt(grav*qk.h) + sqrt(qk.hx*qk.hx + qk.hy*qk.hy) / qk.h;
    double a  = dt * std::max(vm*vm, vk*vk);

    // Helper values
    double scale = 0.5 * n_length;
    double gh2   = 0.5 * grav * (qm.h*qm.h + qk.h*qk.h);

    // Simple flux with dissipation for stability
    return QVar(scale * (wm*qm.h  + wk*qk.h)           - a * (qm.h  - qk.h),
                scale * (wm*qm.hx + wk*qk.hx + gh2*nx) - a * (qm.hx - qk.hx),
                scale * (wm*qm.hy + wk*qk.hy + gh2*ny) - a * (qm.hy - qk.hy));
  }
};

/** Node position function object for use in the SDLViewer. */
struct NodePosition {
  template <typename NODE>
  Point operator()(const NODE& n) const {
    return n.position();
  }
};

/** Calculates the position of the vertex */
struct VertexPosition {
	Point operator()(Vertex& v) const {
		return Point(v.position().x, v.position().y, v.value().h - 1);
	}
};

/** Integrate a hyperbolic conservation law defined over the mesh m
 * with flux functor f by dt in time.
 */
template <typename MESH, typename FLUX>
double hyperbolic_step(MESH& m, FLUX& f, double t, size_t n_triangles, double dt) {

	// Step the finite volume model in time by dt.
	std::vector<QVar> flux_list;
	flux_list.reserve(n_triangles);

	for (auto tri = m.triangles_begin(); tri != m.triangles_end(); ++tri ) {

		// initial values for this triangle
		Triangle this_tri = *tri;
		QVar flux = QVar(0,0,0);
		auto these_edges = m.edges(this_tri);

		// Integrating over edges of triangle
		for (auto eit = these_edges.begin(); eit != these_edges.end(); ++eit) {

			// Grab the edge
			Edge this_edge = *eit;

			// getting the normal
			Point normal;
			normal = this_edge.value().normal_;

			// getting the flux
			QVar this_qvar = this_tri.value().qvar_;
			QVar other_q; 

			// getting triangle on edge that's not this one
			triangle_set others = this_edge.other_triangle(this_tri);
			if (others.size() == 0)
				other_q = QVar(this_qvar.h, 0, 0);
			else {
				other_q = (*others.begin()).value().qvar_;
				if (this_tri > *others.begin()) normal = -normal;
			}
			
			// Increment the overall flux for this triangle by the flux 
			// contribution across this edge
			flux = flux + f( normal.x, normal.y, dt,
					   this_tri.value().qvar_, other_q);
		}
			
		// saving the flux
		flux_list.push_back(flux);
	}

	// now updating all fluxes
	auto fli = flux_list.begin();
	for (auto tri = m.triangles_begin(); tri != m.triangles_end(); ++tri, ++fli ) {
		(*tri).value().qvar_ -= (*fli) * (dt / (*tri).value().area_);
	}

	return t + dt;
}

/** Convert the triangle-averaged values to node-averaged values for viewing. */
template <typename MESH>
void post_process(MESH& m) {
	for (auto vit = m.vertex_begin(); vit != m.vertex_end(); ++vit) {
		double total_h = 0;
		int count = 0;
		for (auto adj = (*vit).triangles_begin(); 
		    	adj != (*vit).triangles_end(); ++adj) {
			total_h += (*adj).value().qvar_.h;
			++count;
		}
		(*vit).value().h = total_h / count;
	}
}

/** Creates initial conditions for dam simulation */
double step(double x) { if( x < 0 ) return 1; else return 0; }

/** Creates initial conditions for dam simulation */
struct DamInitializer {
	QVar operator()(Point p) { return QVar(1 + 0.75*step(p.x), 0, 0); }
};

/** Creates initial conditions for the giant wave simulation */
struct GiantWaveInitializer {
	QVar operator()(Point p) {
		double temp = (p.x - 0.75)*(p.x - 0.75) + p.y*p.y - 0.15*0.15;
		return QVar(1 + 0.75*step(temp), 0, 0);
	}
};

/** Creates initial conditions for the pebble simulation */
struct PebbleInitializer {
	QVar operator()(Point p) {
		double temp = 1-0.75 * std::exp(-80*((p.x-0.75)*(p.x - 0.75) + p.y*p.y));
		return QVar(temp, 0, 0);
	}
};

/* Initialzes the qvar_ values for all the triangles in the @a mesh using the
 * 	given initializer functor
 * @pre INIT is a type that implements the operator()(Point p) and returns
 * 	QVar with the corresponding initial condition
 */
template <typename INIT>
void initialize_mesh(MeshType& mesh, INIT initializer) {
  // Precompute the Qvars in the triangles
  for(auto it = mesh.triangles_begin(); it != mesh.triangles_end(); ++it) {
  	Triangle this_triangle = (*it);
	this_triangle.value().qvar_ = initializer(this_triangle.position());
  }
}

int main(int argc, char* argv[])
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: shallow_water NODES_FILE TRIS_FILE INIT\n";
    exit(1);
  }

  MeshType mesh;
  std::vector<typename MeshType::node_type> mesh_node;
  // Read all Points and add them to the Mesh
  std::ifstream nodes_file(argv[1]);
  Point p;
  while (CS207::getline_parsed(nodes_file, p)) {
    mesh_node.push_back(mesh.add_node(p));
  }

  // Read all mesh triangles and add them to the Mesh
  std::ifstream tris_file(argv[2]);
  std::array<int,3> t;
  while (CS207::getline_parsed(tris_file, t)) {
    mesh.add_triangle(mesh_node[t[0]], mesh_node[t[1]], mesh_node[t[2]]);
  }

  // Print out the stats
  std::cout << mesh.num_nodes() << " "
            << mesh.num_edges() << " "
            << mesh.num_triangles() << std::endl;

  // Set the initial conditions
  // Perform any needed precomputation

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  viewer.launch();

  auto vertex_map = viewer.empty_vertex_map(mesh);
  viewer.add_nodes(mesh.vertex_begin(), mesh.vertex_end(),
                   CS207::DefaultColor(), NodePosition(), vertex_map);
  viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), vertex_map);


  /* Uncomment this code to draw a graph of links between triangles of the graph
  // Add triangles to the the graph to test the triangle iterator
  auto triangle_map = viewer.empty_triangle_map(mesh);
  viewer.add_nodes(mesh.triangles_begin(), mesh.triangles_end(),
                   CS207::RedColor(), NodePosition(), triangle_map);
  viewer.add_edges(mesh.link_begin(), mesh.link_end(), triangle_map);
  viewer.center_view();
  */


  // Initialize an initial condition structures
  DamInitializer di;
  GiantWaveInitializer pi;
  PebbleInitializer pebi;

  // select initial conditions from 3rd command line argument
  std::string choice;
  if (argc < 4)
	  choice = ""; // go to default
  else
	  choice = argv[3];
  if (choice.compare("dam") == 0)
	  initialize_mesh(mesh, di);
  else if (choice.compare("pebble") == 0)
	  initialize_mesh(mesh, pebi);
  else
	  initialize_mesh(mesh, pi); // wave is default case


  // Compute Timestep
  // CFL stability condition requires dt <= dx / max|velocity|
  // For the shallow water equations with u = v = 0 initial conditions
  //   we can compute the minimum edge length and maximum original water height
  //   to set the time-step
  // Compute the minimum edge length and maximum water height for computing dt
  double max_height = 0; // hardcoded from knowledge about initial conditions
  for (auto it = mesh.triangles_begin(); it != mesh.triangles_end(); ++it) {
	  double this_height = (*it).value().qvar_.h;
	  if (this_height > max_height)
		  max_height = this_height;
  }

  double min_edge_length = 10000;
  for(auto it = mesh.edge_begin(); it != mesh.edge_end(); ++it) {
	double this_edge_length = (*it).length();
  	if(this_edge_length < min_edge_length)
		min_edge_length = this_edge_length;
  }
  double dt = 0.25 * min_edge_length / (sqrt(grav * max_height));
  double t_start = 0;
  double t_end = 10;

  db("time step:", dt);
  db("min edge length: ", min_edge_length);
  db("max height: ", max_height);

  // Preconstruct a Flux functor
  EdgeFluxCalculator f;

  // cache the areas of triangles  	
  for(auto it = mesh.triangles_begin(); it != mesh.triangles_end(); ++it) {
	(*it).value().area_ = (*it).area();
  }
  // Compute the normals across each of the edges
  for(auto it = mesh.edge_begin(); it != mesh.edge_end(); ++it) {
  		Edge this_edge = *it;
		this_edge.value().normal_ = mesh.normal(this_edge);
  }

  size_t num_triangles = mesh.num_triangles();

  // Begin the time stepping
  for (double t = t_start; t < t_end; t += dt) {

    // Step forward in time with forward Euler
    hyperbolic_step(mesh, f, t, num_triangles, dt);

    // Update node values with triangle-averaged values
    post_process(mesh);

    // Update the viewer with new node positions
    viewer.add_nodes(mesh.vertex_begin(), mesh.vertex_end(),
                     CS207::DefaultColor(), VertexPosition(), vertex_map);
    viewer.add_nodes(mesh.vertex_begin(), mesh.vertex_end(),
                     CS207::DefaultColor(), VertexPosition(), vertex_map);
    viewer.set_label(t);

  }


  return 0;
}