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

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"

#include "Point.hpp"
#include "Mesh.hpp"

template <typename in>
void db(in s) {
	// std::cout << s << std::endl;
}

void db(Point p) {
	// std::cout << "(" << p.x << ", " << p.y << ", " << p.z << ")" << std::endl;
}

/** Water column characteristics */
// This goes in as the UserTriangleData
struct QVar {
  double h;   // Height of column
  double hx;  // Height times average x velocity of column
  double hy;  // Height times average y velocity of column

  /** Default constructor.
   *
   * A default water column is 1 unit high with no velocity. */
  QVar()
      : h(1), hx(0), hy(0) {
  }
  /** Construct the given water column. */
  QVar(double h_, double hx_, double hy_)
      : h(h_), hx(hx_), hy(hy_) {
  }

  /** sum up two Qs*/
  QVar operator+(QVar o) const {
	  return QVar(h + o.h, hx + o.hx, hy + o.hy);
  }

  QVar operator-(QVar o) const {
	  return QVar(h - o.h, hx - o.hx, hy - o.hy);
  }

  void operator+=(QVar o) {
	  *this = (*this + o);
  }

  void operator-=(QVar o) {
	  *this = (*this - o);
  }

  QVar operator*(double s) const {
	  return QVar(h*s, hx*s, hy*s);
  }

  QVar operator/(double s) const {
	  return *this * 1/s;
  }
};

// Standard gravity (average gravity at Earth's surface) in meters/sec^2
static constexpr double grav = 9.80665;

typedef struct my_triangle_data {
	// TODO: Qvar stuff
	double area_;
	QVar qvar_;
} my_triangle_data;

typedef struct my_link_data {
	Point normal_;
} my_link_data;

typedef struct my_vertex_data {
	double h;
} my_vertex_data;


// HW4B: Placeholder for Mesh Type!
// Define NodeData, EdgeData, TriData, etc
// or redefine for your particular Mesh
typedef Mesh<my_vertex_data, my_link_data, my_triangle_data> MeshType;
typedef MeshType::Link Link;
typedef MeshType::Triangle Triangle;
typedef MeshType::Vertex Vertex;
typedef MeshType::triangle_iterator triangle_iterator;


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

/** For debugging purposes */
typedef struct Tilt {
	QVar operator()(double nx, double ny, double dt, const QVar& qk, const QVar& qm) {
		return qk;
	}
} Tilt;

/** Node position function object for use in the SDLViewer. */
struct NodePosition {
  template <typename NODE>
  Point operator()(const NODE& n) const {
    return n.position();
  }
};

struct VertexPosition {
	Point operator()(Vertex& v) const {
		return Point(v.position().x, v.position().y, v.value().h);
	}
};

/** Integrate a hyperbolic conservation law defined over the mesh m
 * with flux functor f by dt in time.
 */
template <typename MESH, typename FLUX>
double hyperbolic_step(MESH& m, FLUX& f, double t, double dt) {
	// HW4B: YOUR CODE HERE
	// Step the finite volume model in time by dt.

	std::vector<QVar> flux_list;
	for (auto tri = m.triangles_begin(); tri != m.triangles_end(); ++tri ) {
		QVar flux = QVar(0,0,0);
		for (auto link = (*tri).link_begin(); link != (*tri).link_end(); ++link) {
			// getting the normal
			Point normal;
			if(*tri < (*link).triangle2())
				normal = (*link).value().normal_;
			else
				normal = -(*link).value().normal_;

			normal = m.normal(*tri, (*link).triangle2());

			db("getting normal");
			db(normal.x);
			db(normal.y);
			db("this triangle");
			db((*tri).value().qvar_.h);
			db((*tri).value().qvar_.hx);
			db((*tri).value().qvar_.hy);
			db("other triangle");
			db((*link).triangle2().value().qvar_.h);
			db((*link).triangle2().value().qvar_.hx);
			db((*link).triangle2().value().qvar_.hy);

			// getting the flux
			QVar temp = f(normal.x, normal.y, dt, (*tri).value().qvar_, (*link).triangle2().value().qvar_);
			db("temp");
			db(temp.h);
			db(temp.hx);
			db(temp.hy);
			flux = flux + temp;
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
		for (auto adj = (*vit).triangles_begin(); adj != (*vit).triangles_end(); ++adj) {
			total_h += (*adj).value().qvar_.h;
			++count;
		}
		// (*vit).value().h = total_h / count;
		(*vit).value().h += (*vit).position().y;
	}
}



int main(int argc, char* argv[])
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: shallow_water NODES_FILE TRIS_FILE\n";
    exit(1);
  }

  MeshType mesh;
  // HW4B: Need node_type before this can be used!
#if 1
  std::vector<typename MeshType::node_type> mesh_node;
#endif
  std::cerr << "started 1" << std::endl;
  // Read all Points and add them to the Mesh
  std::ifstream nodes_file(argv[1]);
  Point p;
  while (CS207::getline_parsed(nodes_file, p)) {
    // HW4B: Need to implement add_node before this can be used!
#if 1
    mesh_node.push_back(mesh.add_node(p));
#endif
  }
  std::cerr << "started 2" << std::endl;

  // Read all mesh triangles and add them to the Mesh
  std::ifstream tris_file(argv[2]);
  std::array<int,3> t;
  while (CS207::getline_parsed(tris_file, t)) {
    // HW4B: Need to implement add_triangle before this can be used!
#if 1
    mesh.add_triangle(mesh_node[t[0]], mesh_node[t[1]], mesh_node[t[2]]);
#endif
  }

  std::cerr << "started 3" << std::endl;

  // Print out the stats
  std::cout << mesh.num_nodes() << " "
            << mesh.num_edges() << " "
            << mesh.num_triangles() << std::endl;

  // HW4B Initialization
  // Set the initial conditions
  // Perform any needed precomputation

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  viewer.launch();

  // HW4B: Need to define Mesh::node_type and node/edge iterator
  // before these can be used!
#if 1
  auto vertex_map = viewer.empty_vertex_map(mesh);
  viewer.add_nodes(mesh.vertex_begin(), mesh.vertex_end(),
                   CS207::DefaultColor(), NodePosition(), vertex_map);
  viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), vertex_map);
#endif
// Add triangles to the the graph to test the triangle iterator
#if 1
  auto triangle_map = viewer.empty_triangle_map(mesh);
  viewer.add_nodes(mesh.triangles_begin(), mesh.triangles_end(),
                   CS207::RedColor(), NodePosition(), triangle_map);
  viewer.add_edges(mesh.link_begin(), mesh.link_end(), triangle_map);
#endif
  viewer.center_view();


  // HW4B: Timestep
  // CFL stability condition requires dt <= dx / max|velocity|
  // For the shallow water equations with u = v = 0 initial conditions
  //   we can compute the minimum edge length and maximum original water height
  //   to set the time-step
  // Compute the minimum edge length and maximum water height for computing dt

#if 0
  double dt = 0.25 * min_edge_length / (sqrt(grav * max_height));
#else
  // Placeholder!! Delete me when min_edge_length and max_height can be computed!
  double dt = 0.1;
#endif
  double t_start = 0;
  double t_end = 10;

  // Preconstruct a Flux functor
  EdgeFluxCalculator f;

  // cache the areas of triangles  	
  for(auto it = mesh.triangles_begin(); it != mesh.triangles_end(); ++it) {
	(*it).value().area_ = (*it).area();
  }
  // Compute the normals for all of the links
  for(auto it = mesh.link_begin(); it != mesh.link_end(); ++it) {
  	Link this_link = (*it);

  	Triangle t1 = this_link.triangle1();
  	Triangle t2 = this_link.triangle2();

	if(this_link.triangle1() < this_link.triangle2()) {
		this_link.value().normal_ = mesh.normal(this_link.triangle1(), this_link.triangle2());
	}
	else {
		this_link.value().normal_ = mesh.normal(this_link.triangle2(), this_link.triangle1());
	}

	db("precompute done");
	db(mesh.normal(this_link.triangle1(), this_link.triangle2()));
	db(this_link.value().normal_);


  }
  // Precompute the Qvars in the triangles
  for(auto it = mesh.triangles_begin(); it != mesh.triangles_end(); ++it) {
  	Triangle this_triangle = (*it);
	if( this_triangle.position().x < 0)
		this_triangle.value().qvar_ = QVar(1.75, 0, 0);
	else
		this_triangle.value().qvar_ = QVar(1.0, 0, 0);
  }

  for (auto it = mesh.triangles_begin(); it != mesh.triangles_end(); ++it) {
	  db("precomputed:");
	  db((*it).value().qvar_.h);
  }

  // Begin the time stepping
  for (double t = t_start; t < t_end; t += dt) {
    // Step forward in time with forward Euler
    hyperbolic_step(mesh, f, t, dt);

    for (auto it = mesh.triangles_begin(); it != mesh.triangles_end(); ++it) {
    	  db("after hyperbolic:");
    	  db((*it).value().qvar_.h);
      }

    // Update node values with triangle-averaged values
    post_process(mesh);

    for (auto it = mesh.triangles_begin(); it != mesh.triangles_end(); ++it) {
    	  db("after post process:");
    	  db((*it).value().qvar_.h);
      }

    // Update the viewer with new node positions
    // HW4B: Need to define node_iterators before these can be used!
#if 1
    viewer.add_nodes(mesh.vertex_begin(), mesh.vertex_end(),
                     CS207::DefaultColor(), VertexPosition(), vertex_map);
#endif
    viewer.set_label(t);

    // These lines slow down the animation for small meshes.
    // Feel free to remove them or tweak the constants.
    if (mesh.num_nodes() < 100)
      CS207::sleep(0.05);
  }


  return 0;
}
