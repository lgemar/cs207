#ifndef CS207_GRAPH_HPP
#define CS207_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

#include "CS207/Util.hpp"
#include "Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
 private:


 public:

  /////////////////////////////
  // PUBLIC TYPE DEFINITIONS //
  /////////////////////////////

  /** Type of this graph. */
  typedef Graph graph_type;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  typedef Node node_type;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  typedef Edge edge_type;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  typedef unsigned size_type;

  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
  }
  /** Default destructor */
  ~Graph() = default;

  /////////////
  // General //
  /////////////

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
  	return num_nodes_;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
  }

  /////////////////
  // GRAPH NODES //
  /////////////////

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node {
   public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Graph::node_type x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */
    Node() {
    }

    /** Return this node's position. */
    const Point& position() const {
	  return graph_->nodes_[uid_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return uid_;
    }

    /** Test whether this node and @a x are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& x) const {
		return (uid_ == x.index());
    }

    /** Test whether this node is less than @a x in the global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& x) const {
		return (uid_ < x.index());
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
	const Graph* graph_;
	size_type uid_ = 0;
	// Construct a node as just a pointer to the graph and an id number
	Node(const Graph* graph, size_type uid) : graph_(graph), uid_(uid) {
	}
  };

  /** Synonym for size(). */
  size_type num_nodes() const {
		  return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new size() == old size() + 1
   * @post result_node.index() == old size()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position) {
	nodes_.push_back(position);
  	Node new_node = Node(this, num_nodes_);
	++num_nodes_;
	return new_node;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, i);
  }

  /////////////////
  // GRAPH EDGES //
  /////////////////

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
		return Node(graph_, graph_->edges_[index_].node1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
		return Node(graph_, graph_->edges_[index_].node2);
    }

    /** Test whether this edge and @a x are equal.
     *
     * Equal edges are from the same graph and have the same nodes.
     */
    bool operator==(const Edge& x) const {
		bool equal_check1 = (node1() == x.node1() && node2() == x.node2());
		bool equal_check2 = (node2() == x.node1() && node1() == x.node2());
		return (equal_check1 || equal_check2);
    }

    /** Test whether this edge is less than @a x in the global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The edge ordering relation must obey trichotomy: For any two edges x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Edge& x) const {
		return (index() < x.index());
    }

	size_type index() const {
		return index_;
	}

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
	
	// data members
	const Graph* graph_;
	size_type index_;

	// Constructor available to the Graph class for constructing edges
  	Edge(const Graph* graph, size_type index) : graph_(graph), index_(index) {
	}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
  	return num_edges_;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b) {
  	edge_data new_edge;
	size_type edge_index;

	edge_index = has_edge(a, b);
	if (edge_index == num_edges_) {
		new_edge.node1 = a.index();
		new_edge.node2 = b.index();
		edges_.push_back(new_edge);
		++num_edges_;
	}
	return Edge(this, edge_index);
  }

  /** Check to see if there is an edge between nodes a and b
   * @a a and @a b are distince valid nodes in the graph
   * @return index of the edge if it exists, else return the number of edges
   */
   size_type has_edge(const Node& a, const Node& b) {
	edge_data edge_of_interest;
	bool equal_check1, equal_check2;
	size_type i;
   	for (i = 0; i < num_edges_; i++) {
		edge_of_interest = edges_[i];
		equal_check1 = (edge_of_interest.node1 == a.index() && 
						edge_of_interest.node2 == b.index());
		equal_check2 = (edge_of_interest.node2 == a.index() &&
						edge_of_interest.node1 == b.index());
		if (equal_check1 || equal_check2) {
			return i;
		}
	}
	return num_edges_; // if no edge matches the two nodes, then return false
   }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this, i);
  }

 private:
	typedef struct edge_data {
		size_type node1;
		size_type node2;
	} edge_data;

	size_type num_nodes_ = 0;
	size_type num_edges_ = 0;
 	std::vector<Point> nodes_;
	std::vector<edge_data> edges_;
};
#endif
