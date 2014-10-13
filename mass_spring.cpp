/**
 * @file mass_spring.cpp
 * Implementation of mass-spring system using Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <fstream>
#include <numeric>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"

#include "Graph.hpp"
#include "Point.hpp"


// Define "value_type" of Point to be a scalar
typedef Point::value_type scalar;

// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point velocity;  //< Node velocity
  double mass;     //< Node mass
};

struct EdgeData {
	int test;
};

// HW2 #1 YOUR CODE HERE
// Define your Graph type
typedef Graph<NodeData, EdgeData> GraphType;
typedef typename GraphType::node_type Node;
typedef typename GraphType::edge_type Edge;

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on Node
 *           at time @a t.
 */
template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
  // Compute the {n+1} node positions
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().velocity * dt;
  }

  // Compute the {n+1} node velocities
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().velocity += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}


/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force being applied to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. 
   * The force on a given node is computed by adding the forces on it from 
   * all its connected nodes. Represent the node adjacency list by A. Then, 
   * the force on a given node n is given by the sum of forces from nodes
   * n0, n1, ..., ni, ..., nm in the node adjacency list. To calculate the 
   * force from any given node, use the spring equation. This is given simply
   * as F = -kx, where k is the spring constant and x the displacement. The
   * displacement is a vector. Its direction is computed by finding the vector
   * difference between @a n and node xi in the adjacency list. The magnitude
   * is calculated by finding the Euclidean distance between the nodes in the 
   * adjacency list and the spring length, L. 
   * @returns a Point that represents the force vector
   */
  Point operator()(Node n, double t) {
	// Initialize variables
	(void) t;//suppress compiler warning
	Node adjacent_node;
  	scalar K = 100.0; // Spring constant
	scalar displacement; // displacement from spring rest-length
	Point direction; // direction of the force
	Point total_force;
	Point xi, xj; // xi: position of node n; xj position of adjacent node

	// Calculate force on node
  	if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
		return Point(0, 0, 0);
	}
	total_force = Point(0, 0, -grav * n.value().mass);
	xi = n.position();
	for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
		adjacent_node = (*it).node2();
		xj = adjacent_node.position();
		displacement = distance(xi, xj) - (*it).length();
		direction = (xi - xj) / distance(xi, xj);
		total_force += -K * displacement * direction;
	}
	return total_force;
  }
};

class Stimulus {
	public: 
		virtual Point apply(Node, double)=0;
};

class Force {
	public: 
		typedef typename std::list<Stimulus*> f_composition;
		f_composition forces_;
		Force(Stimulus* s) {
			forces_.push_front(s);
		}

		Force(f_composition composite) : forces_(composite) {
		}
		
		Point operator()(Node n, double t) const {
			Point total_force = Point(0,0,0);
			for(auto it = forces_.begin(); it != forces_.end(); ++it) {
				total_force = total_force + (*it)->apply(n, t);
			}
			return total_force;
		}

		Force operator+(Force f) const {
			f_composition force_list = f.forces_;
			for(auto it = forces_.begin(); it != forces_.end(); ++it) {
				force_list.push_front(*it);	
			}
			return Force(force_list);
		}
};

class GravityForce : public Stimulus {
	public:
		virtual Point apply(Node n, double t) {
			(void) t;
			Point total_force;
			// Calculate force on node
			if (n.position() == Point(0, 0, 0) || 
							n.position() == Point(1, 0, 0)) {
				return Point(0, 0, 0);
			}
			total_force = Point(0, 0, -grav * n.value().mass);
			return total_force;
		}
};

class MassSpringForce : public Stimulus {
  public: 
	  virtual Point apply(Node n, double t) {
		// Initialize variables
		(void) t;//suppress compiler warning
		Node adjacent_node;
		scalar K = 100.0; // Spring constant
		scalar displacement; // displacement from spring rest-length
		Point direction; // direction of the force
		Point total_force = Point(0, 0, 0);
		Point xi, xj; // xi: position of node n; xj position of adjacent node

		// Calculate force on node
		if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
			return Point(0, 0, 0);
		}
		xi = n.position();
		for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
			adjacent_node = (*it).node2();
			xj = adjacent_node.position();
			displacement = distance(xi, xj) - (*it).length();
			direction = (xi - xj) / distance(xi, xj);
			total_force += -K * displacement * direction;
		}
		return total_force;
	  }
};

int main(int argc, char** argv) {
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a graph
  GraphType graph;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  std::vector<Node> nodes;
  Point p;
  while (CS207::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CS207::getline_parsed(tets_file, t)) {
    for (unsigned i = 1; i < t.size(); ++i) {
      graph.add_edge(nodes[t[0]], nodes[t[1]]);
      graph.add_edge(nodes[t[0]], nodes[t[2]]);
      // Diagonal edges: include as of HW2 #2
      graph.add_edge(nodes[t[0]], nodes[t[3]]);
      graph.add_edge(nodes[t[1]], nodes[t[2]]);
      graph.add_edge(nodes[t[1]], nodes[t[3]]);
      graph.add_edge(nodes[t[2]], nodes[t[3]]);
    }
  }

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
  for( auto it = graph.node_begin(); it != graph.node_end(); ++it ) {
  	// Initialize velocities
	(*it).value().velocity = Point(0, 0, 0);
	// Initialize mass
	(*it).value().mass = (scalar) 1 / graph.size();
	// Initialize edge lengths
  }
  // Construct Forces/Constraints

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  auto node_map = viewer.empty_node_map(graph);
  viewer.launch();

  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  viewer.center_view();

  // Begin the mass-spring simulation
  double dt = 0.001;
  double t_start = 0.0;
  double t_end   = 5.0;

  // Construct the problem 1 force using new Force structure

  GravityForce g_force;
  MassSpringForce spring_force;
  Force gravity_f (&g_force);
  Force spring_f (&spring_force);
  Force problem1_f = spring_f + gravity_f;

  for (double t = t_start; t < t_end; t += dt) {
    //std::cout << "t = " << t << std::endl;
    symp_euler_step(graph, t, dt, problem1_f);

    // Update viewer with nodes' new positions
    viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
    viewer.set_label(t);

    // These lines slow down the animation for small graphs, like grid0_*.
    // Feel free to remove them or tweak the constants.
    if (graph.size() < 100)
      CS207::sleep(0.001);
  }

  return 0;
}
