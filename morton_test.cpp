/**
 * @file morton_test.cpp
 * Example of a spatial search using Morton codes.
 *
 * In this example, SpaceSearcher is used to iterate over the Nodes contained
 * inside a small BoundingBox close to the center of a Graph. The performance
 * of SpaceSearcher (Morton codes) is compared to that of a linear search.
 * This example can be added to the Makefile and executed as usual:
 *
 * $ make morton_test
 * $ ./morton_test data/large.*
 *
 * 37813 217974
 * Morton search: Found 49 items in 1.3977e-05 seconds.
 * Linear search: Found 49 items in 0.000195776 seconds.
 * Speedup: 14.007
 *
 * The speedup depends to a large extent on the application, and generally
 * becomes worse as the size of the BoundingBox increases. In other words,
 * if the user needs to iterate over most of the items in a Graph or Mesh anyway,
 * then a linear search may actually be more efficient in some cases. But for
 * smaller BoundingBoxes, containing ~20% of the items or less, SpaceSearcher
 * can provide a significant boost in performance. This makes SpaceSearcher
 * an ideal tool for collision detection between small objects and for nearest-
 * neighbor searches, among other possible applications.
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
#include "Graph.hpp"
#include "SpaceSearcher.hpp"

// Define Graph type
typedef Graph<char, char> GraphType;
typedef typename GraphType::node_type Node;
typedef typename GraphType::edge_type Edge;

/** Simple mapping from a node to its corresponding point. */
struct NodeToPoint {
  template <typename NODE>
  Point operator()(const NODE& n) {
    return n.position();
  }
};

// Define spatial searcher type
typedef SpaceSearcher<Node, NodeToPoint> SpaceSearcherType;

int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a Graph
  GraphType graph;
  std::vector<typename GraphType::node_type> nodes;

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

  // Print number of nodes and edges
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Create spatial searcher for nodes
  SpaceSearcherType s(graph.node_begin(), graph.node_end(), NodeToPoint());

  // Find points inside this bounding box
  BoundingBox bb(Point(-0.1, -0.1, -0.05), Point(0.1, 0.1, 0.05));

  // For performance checks
  CS207::Clock clock;

  // Iterate over nodes in the bounding box using SpaceSearcher.
  std::size_t count = 0;
  clock.start();
  for (auto it = s.begin(bb); it != s.end(bb); ++it)
    ++count;
  double time_morton = clock.seconds();
  std::cout << "Morton search: Found " << count << " items in " <<
      time_morton << " seconds." << std::endl;

  // Compare with linear search.
  count = 0;
  clock.start();
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    if (bb.contains((*it).position()))
      ++count;
  }
  double time_linear = clock.seconds();
  std::cout << "Linear search: Found " << count << " items in " <<
      time_linear << " seconds." << std::endl;
  std::cout << "Speedup: " << time_linear/time_morton << std::endl;

  return 0;
}
