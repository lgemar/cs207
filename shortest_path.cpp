/**
 * @file shortest_path.cpp
 * Test script for using our templated Graph to determine shortest paths.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <vector>
#include <fstream>
#include <math.h>
#include <queue>
#include <set>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"

#include "Graph.hpp"


/** Define custom operator to return a color object for a node */
struct MyColorFunc {
	int lp_;

	/** Constructor */
	MyColorFunc(int longest_path) : lp_(longest_path) {
	}

	template <typename NODE>
	// Return a color object to color the graph
	CS207::Color operator()(const NODE& node) {
		// return CS207::Color::make_heat((float)node.value() / lp_);
		return CS207::Color::make_heat(0.5);
	}
};

/** Comparator that compares the distance from a given point p.
 */
struct MyComparator {
   Point p_;
   MyComparator(const Point& p) : p_(p) {
   };

   template <typename NODE>
   // Return true if node1 is closer to p than node2
   bool operator()(const NODE& node1, const NODE& node2) const {
		double dist1 = sqrt(pow(node1.position().x, 2) 
									+ pow(node1.position().y, 2));
		double dist2 = sqrt(pow(node2.position().x, 2) 
									+ pow(node2.position().y, 2));
		return dist1 < dist2;
   }
};

template <typename NODE>
int breadth_search(Graph<int>& g, 
				   Graph<int>::node_type& node, 
				   std::set<NODE> visited) {
	Graph<int>::node_type adjacent_node;
	std::queue<NODE> to_visit;
	int max = 0;
	for (auto it = node.edge_begin(); it != node.edge_end(); ++it) {
		adjacent_node = (*it).node2();
		if (visited.find(adjacent_node) == visited.end()) {
			double real_distance = sqrt(pow(node.position().x, 2) + 
									   pow(adjacent_node.position().y, 2));
			int integer_distance = (int) (real_distance * 10);
			adjacent_node.value() = integer_distance + node.value();
			if (adjacent_node.value() > max) {
				max = adjacent_node.value();
			}
			visited.insert(adjacent_node);
			to_visit.push(adjacent_node);
		}
	}
	while (!to_visit.empty()) {
		adjacent_node = to_visit.front();
		int temp_max = breadth_search(g, adjacent_node, visited);
		if (temp_max > max) {
			max = temp_max;
		}
		to_visit.pop();
	}
	std::cout << "max: " << max << std::endl;
	return max;
}

/** Calculate shortest path lengths in @a g from the nearest node to @a point.
 * @param[in,out] g Input graph
 * @param[in] point Point to find the nearest node to.
 * @post Graph has modified node values indicating the minimum path length
 *           to the nearest node to @a point
 * @post Graph nodes that are unreachable to the nearest node to @a point have
 *           the value() -1.
 * @return The maximum path length found.
 *
 * Finds the nearest node to @a point and treats that as the root node for a
 * breadth first search.
 * This sets node's value() to the length of the shortest path to
 * the root node. The root's value() is 0. Nodes unreachable from
 * the root have value() -1.
 */
int shortest_path_lengths(Graph<int>& g, const Point& point) {
	int longest_path;
	std::set<Graph<int>::node_type> visited;
	auto closest = std::min_element(g.node_begin(), 
									g.node_end(), MyComparator(point));
	Graph<int>::node_type closest_node = *closest;
	// Set the root value to 0
	closest_node.value() = 0;
	// Put the root in the "visited" list
	visited.insert(closest_node);
	// Recursively search all of the roots nodes
	longest_path = breadth_search(g, closest_node, visited);
	return longest_path;
}


int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a Graph
  typedef Graph<int> GraphType;
  GraphType graph;
  std::vector<GraphType::node_type> nodes;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  while (CS207::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CS207::getline_parsed(tets_file, t))
    for (unsigned i = 1; i < t.size(); ++i)
      for (unsigned j = 0; j < i; ++j)
        graph.add_edge(nodes[t[i]], nodes[t[j]]);

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  viewer.launch();

  // Create empty node map
  auto node_map = viewer.empty_node_map(graph);

  // Use shortest_path_lengths to set the node values to the path lengths
  int longest_path = shortest_path_lengths(graph, Point(-1, 0, 1));

  // Construct a Color functor and view with the SDLViewer
  viewer.add_nodes(graph.node_begin(), 
  				   graph.node_end(), 
				   MyColorFunc(longest_path), 
				   node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
  viewer.center_view();
  return 0;
}
