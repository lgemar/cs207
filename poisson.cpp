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
	bool boundary;
} node_data;

typedef struct EdgeData {
} edge_data;

typedef Graph<node_data,edge_data> GraphType;
typedef typename GraphType::node_type Node;
typedef typename GraphType::edge_type Edge;

// Define a GraphSymmetricMatrix that maps
// your Graph concept to MTL's Matrix concept. This shouldn't need to copy or
// modify Graph at all!
class GraphSymmetricMatrix {

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
	 */
};

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
		// Iterate through all the adjacent nodes and flag adjacent nodes
		// as boundary nodes; need to figure out if node is one edge length
		// away from the boundary box
		for (auto adj_it = n.edge_begin(); adj_it != n.edge_end(); ++adj_it) {
			Edge e = *(adj_it);
			Node adj_node = e.node2();
			adj_node.value().boundary = true;
		}
		// Remove the node contained within the bounding box
		remove_list.push_front(n);
	}
  }

  /** Remove nodes in the remove list */
  for (auto it = remove_list.begin(); it != remove_list.end(); ++it)
		g.remove_node((*it));

  return;
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
  
  /** Initialize the boundary values of the nodes */
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
	Node n = (*it);
  	n.value().boundary = false;
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

  // HW3: YOUR CODE HERE
  // Define b using the graph, f, and g.
  // Construct the GraphSymmetricMatrix A using the graph
  // Solve Au = b using MTL.

  // Launch a viewer
  CS207::SDLViewer viewer;
  viewer.launch();

  // Create a graph
  auto node_map = viewer.empty_node_map(graph);
  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
  viewer.center_view();

  return 0;
}
