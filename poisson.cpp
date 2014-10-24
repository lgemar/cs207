/**
 * @file poisson.cpp
 * Test script for treating the Graph as a MTL Matrix
 * and solving a Poisson equation.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles.
 * Second file: Eges (one per line) defined by 2 indices into the point list
 *              of the first file.
 *
 * Launches an SDLViewer to visualize the solution.
 */

// HW3: Need to install/include Boost and MTL in Makefile
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"

#include "Graph.hpp"
#include "Point.hpp"
#include "BoundingBox.hpp"
#include <fstream>
#include <list>

/** Useful type information */
// Define Node data and Edge data types
typedef struct NodeData {
	double poisson;
} node_data;

typedef struct EdgeData {
} edge_data;

typedef Graph<node_data,edge_data> GraphType;
typedef typename GraphType::node_type Node;
typedef typename GraphType::edge_type Edge;

bool on_boundary(const Point& p) {
	bool trial0 = norm_inf(p - Point(0.6, 0.6, 0)) < 0.2;
	bool trial1 = norm_inf(p - Point(-0.6, 0.6, 0)) < 0.2;
	bool trial2 = norm_inf(p - Point(0.6, -0.6, 0)) < 0.2;
	bool trial3 = norm_inf(p - Point(-0.6, -0.6, 0)) < 0.2;
	if( trial0 || trial1 || trial2 || trial3 )
		return true;
	else if (BoundingBox(Point(-0.6, -0.2, -1), Point(0.6, 0.2, 1)).contains(p))
		return true;
	else if ( norm_inf(p) == 1.0 ) 
		return true;
	else 
		return false;
}

// Define a GraphSymmetricMatrix that maps
// your Graph concept to MTL's Matrix concept. This shouldn't need to copy or
// modify Graph at all!
class GraphSymmetricMatrix {

	public:
		const GraphType* g;
		GraphSymmetricMatrix(const GraphType* graph) : g(graph) {
		};
		/** TODO: 
		 * [ ] Fix node iterator problem that causes invalidation of the last
		 *  : node in graph upon node removal
		 *	: when a node gets removed it becomes the end node
		 *	: bug screws up viewer when there are no nodes in the graph
		 * [ ] test whether adjacency iteration is working properly by trying to
		 * 	   remove adjacent nodes. This will prove that the correct nodes
		 * 	   are having their boundary value set to true
		 *	: test did not work. Some adjacencies are removed, others not.
		 * [x] flag whether a node is on a boundary
		 *	: maybe I can do this within "remove_box" ? 
		 *	: iterate along the adjacency list of the removed node and 
		 * 	  tag the "boundary" flag for that node
		 * [x] determine whether two nodes share an edge
		 *	: graph.has_edge(Node i, Node j);  
		 * [x] determine the degree of a node 
		 *	: n.degree();
		 * [x] check viewer to see if remove_box is working with "remove_list" 
		 *		implementation as a workaround for the way "remove_node" 
		 *		invalidates iterators
		 * [x] Use the first page of the Problem 2 description to start 
		 * 		developing this GraphSymmetricMatrix as the linear operator A
		 * [x] Write the operators that determine A and L for some i and j
		 * 		that represent indices into the graph
		 * [x] Write the constructor for a graph symmetric matrix
		 * [x] copy over the Identity matrix multiply functions
		 * [x] Write size, num_rows, and num_cols functions for GraphSymmetric
		 * [x] Modify multiply function to reflect the A(i, j) matrix elements
		 * [x] Get code to compile with new matrix multiplication definition
		 * [x] Figure out why g this class things g is  undefined
		 * [x] Figure out why "mtl" is an undeclared identifier
		 * [x] Solve dummy A * 0 = 0 equation to make sure there are no 
		 	wierd segfaults with indices
		 * [x] construct the b vector in the main function
		 * [x] implement the f and g functions
		 * [x] provide conditional branch in g function base on whether node
		 		is a boundary node
		 * [x] Take advantage of the fact that A is symmetric to reduce
		 	multiplications from O(N^2) to O(N)
				--> This is actually necessary because my solver is too slow
		 * [x] write "on_boundary" function using the fact that nodes on 
		 		in the domain of g are boundary nodes
		 * [x] go through multiply logic to ensure that calculation proceeds 
		 	properly
		 * [ ] go through generation of b matrix and make sure logic is correct
		 */
		 /** Calculate the A(i, j) value for indices i and j */

		/** Helper function to perform delayed evalutation of multiplication. 
		 * Assign::apply(a, b) resolves to an assignment such as a += b, a-=b, 
		 * 	or a = b
		 * @pre size(v) == size(w)
		 */
		template<typename VectorIn, typename VectorOut, typename Assign>
		void mult(const VectorIn& v, VectorOut& w, Assign) const {
			assert( size(v) == size(w) );
			assert( g->size() == size(v) );

			// w must be a vector of all zeros
			size_t highest_index = g->size();
			for (size_t i = 0; i < highest_index; i++) {
				Node n = g->node(i);
				Point pi = n.position();
				double sum = 0.0;
				if( on_boundary( pi ) )
					sum += 1.0 * v[i];
				else
					sum += -1.0 * (double) n.degree() * v[i];
				for(auto it = n.edge_begin(); it != n.edge_end(); ++it) {
					Node adj_node = (*it).node2();
					Point pj = adj_node.position();
					if( !on_boundary(pj) && !on_boundary(pi))
						sum += 1.0 * v[adj_node.index()];
				}
				Assign::apply(w[i], sum);
			}
		}

		/** Matrix-vector multiplication forwards to MTL's lazy mat_cvec_mult
		 * operator */
		template <typename Vector>
		mtl::vec::mat_cvec_multiplier<GraphSymmetricMatrix, Vector> operator*(const Vector& v) const {
			return mtl::vec::mat_cvec_multiplier
						<GraphSymmetricMatrix, Vector>(*this, v);
		}
};

/** Size helper functions for the GraphSymmetricMatrix class */
inline std::size_t size(const GraphSymmetricMatrix& A) { 
	return A.g->size() * A.g->size(); 
}
inline std::size_t num_rows(const GraphSymmetricMatrix& A) { 
	return A.g->size(); 
}
inline std::size_t num_cols(const GraphSymmetricMatrix& A) { 
	return A.g->size(); 
}

/** Traits that mtl uses to determine properties of our Identity matrix */
namespace mtl {

/** Define IdentityMatrix to be a non-scalar type */
namespace ashape {
template<>
struct ashape_aux<GraphSymmetricMatrix> {
	typedef nonscal type;
};
} // end namespace ashape

/** IdentityMatrix implements Collection concept with value type 
 * and size type */
template<>
struct Collection<GraphSymmetricMatrix> {
	typedef double value_type;
	typedef unsigned size_type;
};
}

/** Remove all the nodes in graph @a g whose posiiton is contained within
 * BoundingBox @a bb
 * @post For all i, 0 <= i < @a g.num_nodes(),
 *        not bb.contains(g.node(i).position())
 * @post all nodes adjacent to the removed edges are flagged as boundary nodes
 */
void remove_box(GraphType& g, const BoundingBox& bb) {
  std::list<Node> remove_list;

  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
	Node n = (*it);
	Point p = n.position();
  	if( bb.contains(p) ) {
		// Remove the node contained within the bounding box
		remove_list.push_front(n);
	}
  }

  /** Remove nodes in the remove list */
  for (auto it = remove_list.begin(); it != remove_list.end(); ++it)
		g.remove_node((*it));

  return;
}

double function_f(const Point& p) {
	return 5 * std::cos(norm_1(p));
}

double function_g(const Point& p) {
	bool trial0 = norm_inf(p - Point(0.6, 0.6, 0)) < 0.2;
	bool trial1 = norm_inf(p - Point(-0.6, 0.6, 0)) < 0.2;
	bool trial2 = norm_inf(p - Point(0.6, -0.6, 0)) < 0.2;
	bool trial3 = norm_inf(p - Point(-0.6, -0.6, 0)) < 0.2;
	if( trial0 || trial1 || trial2 || trial3 )
		return -0.2;
	else if (BoundingBox(Point(-0.6, -0.2, -1), Point(0.6, 0.2, 1)).contains(p))
		return 1.0;
	else if ( norm_inf(p) == 1.0 ) 
		return 0;
	assert(false); // g(x) is undefined
}

int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE EDGES_FILE\n";
    exit(1);
  }

  // Define an empty Graph
  GraphType graph;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  std::vector<typename GraphType::node_type> node_vec;
  Point p;
  while (CS207::getline_parsed(nodes_file, p))
    node_vec.push_back(graph.add_node(2*p - Point(1,1,0)));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CS207::getline_parsed(tets_file, t)) {
    graph.add_edge(node_vec[t[0]], node_vec[t[1]]);
    graph.add_edge(node_vec[t[0]], node_vec[t[2]]);
    graph.add_edge(node_vec[t[1]], node_vec[t[3]]);
    graph.add_edge(node_vec[t[2]], node_vec[t[3]]);
  }
  
  // Get the edge length, should be the same for each edge
  double h = graph.edge(0).length();

  // Make holes in our Graph
  remove_box(graph, BoundingBox(Point(-0.8+h,-0.8+h,-1), 
  								Point(-0.4-h,-0.4-h,1)));
  remove_box(graph, BoundingBox(Point( 0.4+h,-0.8+h,-1), 
  								Point( 0.8-h,-0.4-h,1)));
  remove_box(graph, BoundingBox(Point(-0.8+h, 0.4+h,-1), 
  								Point(-0.4-h, 0.8-h,1)));
  remove_box(graph, BoundingBox(Point( 0.4+h, 0.4+h,-1), 
  								Point( 0.8-h, 0.8-h,1)));
  remove_box(graph, BoundingBox(Point(-0.6+h,-0.2+h,-1), 
  								Point( 0.6-h, 0.2-h,1)));

  // Define b using the graph, f, and g.
  // Construct the GraphSymmetricMatrix A using the graph
  // Solve Au = b using MTL.

  /** Set up equations */
  typedef GraphSymmetricMatrix matrix_type;

  // Set up an identity matrix, A
  matrix_type A(&graph);

  // Create a preconditioner
  itl::pc::identity<matrix_type> P(A);

  // Set up a matrix
  mtl::dense_vector<double> x(graph.size(), 0.0), b(graph.size(), 0.0);

  // Make the b vector
  for(auto it = graph.node_begin(); it != graph.node_end(); ++it) {
	Node n = (*it);
  	Point p = n.position();
	double bi = 0.0;
	if( on_boundary(p) ) {
		bi = function_g(p);
	}
	else {
		Edge first_edge = *(graph.edge_begin());
		double h = first_edge.length(); //all edges have same length
		double f_res = function_f(p);
		double adj_sum = 0.0;
		for(auto adj_it = n.edge_begin(); adj_it != n.edge_end(); ++adj_it) {
			Node adj_node = (*adj_it).node2();
			Point adj_p = adj_node.position();
			if( on_boundary(adj_p) )
				adj_sum += function_g(adj_p);
		}
		bi = h * h * f_res - adj_sum;
	}
	b[n.index()] = bi;
  }

  itl::cyclic_iteration<double> iter(b, 10000, 1e-10, 0, 50, std::cout);

  cg(A, x, b, P, iter);

  /** Display the results for sanities sake
  std::cout << "b matrix: " << b << std::endl;
  b = A * x;
  std::cout << "A * x: " << b << std::endl;
  */

  // Set the values of the nodes to the corresponding solutions
  for(auto it = graph.node_begin(); it != graph.node_end(); ++it) {
  	Node n = (*it);
	n.value().poisson = x[n.index()];
  }

  // Launch a viewer
  CS207::SDLViewer viewer;
  viewer.launch();

  // Create a graph
  auto node_map = viewer.empty_node_map(graph);
  viewer.add_nodes(graph.node_begin(), graph.node_end(), 
  					CS207::PoissonColor(), CS207::PoissonPosition(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
  viewer.center_view();

  return 0;
}
