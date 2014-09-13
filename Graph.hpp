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

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  struct node_element;

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
    // HW0: YOUR CODE HERE
    return size_nodes_;
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
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return fetch().position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      for (size_type i = 0; i < graph_->num_nodes(); ++i)
        if (graph_->graph_nodes_.at(i).uid == uid_)
          return size_type(i);
		assert(false);
    }

    /** Test whether this node and @a x are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& x) const {
      // HW0: YOUR CODE HERE
      // (void) x;          // Quiet compiler warning
	  bool result;
	  result = (index() == x.index());
	  return result;
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
      // HW0: YOUR CODE HERE
      //(void) x;           // Quiet compiler warning
	  bool result;
	  result = index() < x.index();
	  return result;
	}

   private:
	// Allow Graph to access Node's private member data and functions.
	friend class Graph;
	// HW0: YOUR CODE HERE
	// Use this space to declare private data members and methods for Node
	// that will not be visible to users, but may be useful within Graph.
	// i.e. Graph needs a way to construct valid Node objects
	// Pointer back to Graph
	Graph* graph_;
	// Identification Number for the Node
	size_type uid_;
	/** Private Constructor */
	Node(const Graph* graph, size_type uid)
			: graph_(const_cast<Graph*>(graph)), uid_(uid) {
	}

	/** Helper method that returns the corresponding element of Graph **/
	node_element& fetch() const {
		size_type idx = index();
  		return graph_->graph_nodes_.at(idx);
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
		  // HW0: YOUR CODE HERE
		  // (void) position;      // Quiet compiler warning
		  // return Node();        // Invalid node
		  node_element element;
		  Node node_element;
		  element.position = position;
		  element.uid = uid();
		  node_element = Node(this, element.uid);
		  element.node = node_element;
		  graph_nodes_.push_back(element);
		  size_nodes_++;
		  return node_element;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
		  // HW0: YOUR CODE HERE
		  // (void) i;             // Quiet compiler warning
		  // return Node();        // Invalid node
		  Node node_element = graph_nodes_.at(i).node;
		  return node_element;
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
				  // HW0: YOUR CODE HERE
				  return Node();      // Invalid Node
		  }

		  /** Return the other node of this Edge */
		  Node node2() const {
				  // HW0: YOUR CODE HERE
				  return Node();      // Invalid Node
		  }

    /** Test whether this edge and @a x are equal.
     *
     * Equal edges are from the same graph and have the same nodes.
     */
    bool operator==(const Edge& x) const {
      // HW0: YOUR CODE HERE
      // (void) x;          // Quiet compiler warning
      // return false;
	  bool check_same1 = (nodea == x.nodea && nodeb == x.nodeb);
	  bool check_same2 = (nodeb == x.nodea && nodea == x.nodeb);
	  return (check_same1 || check_same2);
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
      // HW0: YOUR CODE HERE
      // (void) x;           // Quiet compiler warning
      // return false;
	  Node pair1_min;
	  Node pair2_min;
	  if (nodea < nodeb) 
	  	pair1_min = nodea;
	  else
	  	pair1_min = nodeb;
	  if (x.nodea < x.nodeb)
	  	pair2_min = x.nodea;
	  else
	  	pair2_min = x.nodeb;
	  return pair1_min < pair2_min;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
	Graph* graph_;
	Node nodea;
	Node nodeb;
	/** Private Constructor */
	Edge(const Graph* graph, const Node& a, const Node& b)
			: graph_(const_cast<Graph*>(graph)), nodea(a), nodeb(b) {
	}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    // return 0;
	return size_edges_;
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
    // HW0: YOUR CODE HERE
    // (void) a, (void) b;   // Quiet compiler warning
	// return Edge();        // Invalid Edge
    for (size_type i = 0; i < num_nodes(); ++i) {
	  Edge x;
	  x = edge(i);
	  bool check_same1 = (a == x.nodea && b == x.nodeb);
	  bool check_same2 = (b == x.nodea && a == x.nodeb);
	  bool result = (check_same1 || check_same2);
      if (result)
        return graph_edges_.at(i);
	  Edge temp = Edge(this, a, b);
  	  graph_edges_.push_back(temp);
  	  size_edges_++;
  	  return temp;
	}
	return Edge();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    // (void) i;             // Quiet compiler warning
    // return Edge();        // Invalid Edge
	return graph_edges_.at(i);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  int uid() {
  	static size_type x = size_type(0);
	x++;
	return x;
  }

  struct node_element {
	size_type uid;
	Point position;
	Node node;
  };

  size_type size_nodes_ = 0;
  size_type size_edges_ = 0;

  std::vector<node_element> graph_nodes_; 
  std::vector<Edge> graph_edges_;
};

#endif
