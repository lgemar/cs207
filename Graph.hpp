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
template <typename V>
class Graph {
 private:


 public:

  /////////////////////////////
  // PUBLIC TYPE DEFINITIONS //
  /////////////////////////////

  /** Define the node value type in terms of template parameter */
  typedef V node_value_type;

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

  /** Define uid_type */
  typedef size_type uid_type;

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  typedef NodeIterator node_iterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  typedef EdgeIterator edge_iterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  typedef IncidentIterator incident_iterator;

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
	  return graph_->nodes_[uid_].p;
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

	/** Returns the node_value_type value associated with this Node */
	node_value_type& value() {
		return graph_->nodes_[uid_].v;
	}

	/** Returns the node_value_type value associated with this Node */
	const node_value_type& value() const {
		return graph_->nodes_[uid_].v;
	}

	/** Returns the number of edges associated with this Node */
	size_type degree() const {
		return graph_->nodes_[uid_].degree;
	}

   	/** Returns an iterator to beginning of incident iterator list */
	IncidentIterator& edge_begin() const {
		return IncidentIterator(graph_, uid_).begin();
	}

	/** Returns an iterator to the end of incident iterator list */
	IncidentIterator& edge_end() const {
		return IncidentIterator(graph_, uid_).end();
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
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
	node_data new_node_data;

	new_node_data.p = position;
	new_node_data.v = value;
	new_node_data.degree = 0;
	nodes_.push_back(new_node_data);

  	Node new_node = Node(this, num_nodes_);
	++num_nodes_;
	return new_node;
  }

  /** Determine if this Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
  	return n.index() < num_nodes_;
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
    }

    /** Return a node of this Edge */
    Node node1() const {
		return Node(graph_, uid1_);
	}

    /** Return the other node of this Edge */
    Node node2() const {
		return Node(graph_, uid2_);
    }

    /** Test whether this edge and @a x are equal.
     * @pre RI for Edges must hold: uid1_ < uid2_
     * Equal edges are from the same graph and have the same nodes.
     */
    bool operator==(const Edge& x) const {
		return (node1() == x.node1() && node2() == x.node2());
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
		if (node1() < x.node1())
			return true;
		else
			return node2() < x.node2();
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
	
	// data members
	const Graph* graph_;

	/** RI: uid1_ < uid2 */
	uid_type uid1_;
	uid_type uid2_;

	// Constructor available to the Graph class for constructing edges
  	Edge(const Graph* graph, uid_type node1, uid_type node2) {
		// Make sure that there are no self edges
		assert(node1 != node2);
		graph_ = graph;
		// Enforce RI that uid1_ < uid2_
		if (node1 < node2) {
			uid1_ = node1;
			uid2_ = node2;
		}
		else {
			uid2_ = node1;
			uid1_ = node2;
		}
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
	assert(a.index() == b.index());
  	if (!has_edge(a, b)) {
		nodes_[a.index()].adj.push_back(b.index());
		nodes_[b.index()].adj.push_back(a.index());
	}
	return Edge(this, a.index(), b.index());
  }

  /** Check to see if there is an edge between nodes a and b
   * @a a and @a b are distince valid nodes in the graph
   * @return true if it exists, else return the number of edges
   */
   size_type has_edge(const Node& a, const Node& b) {
	Edge test_edge = Edge(this, a.index(), b.index());
	for (auto it = a.edge_begin(); it != a.edge_end(); ++it) {
		if (test_edge == *it)
			return true;
	}
	// If no equal equivalent edge is found, return false
	return false;
   }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
	assert(i < size());
	EdgeIterator it = edge_begin();
	for (int counter = 0; counter < i; counter++) {
		it++;
	}
	return *it;
  }

  ///////////////
  // Iterators //
  ///////////////

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Node value_type;
    /** Type of pointers to elements. */
    typedef Node* pointer;
    /** Type of references to elements. */
    typedef Node& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

	/** Returns the node referenced by the iterator
	 */
	Node operator*() const{
		return Node(graph_, index_);
	}

	/** Returns the node iterator that points to the next node in the graph
	 */
	NodeIterator& operator++() {
		++index_;
		return *this;
	}

	/** Returns true if this iterator points to the same element as the 
	 * parameterized iterator
	 * 
	 * @param[in] @a it is another node iterator
	 * @returns true if the two iterators point to the same node
	 */
	bool operator==(const NodeIterator& it) const {
		return (Node(graph_, index_) == *it);
	}

   private:
    friend class Graph;
	const Graph* graph_;
	size_type index_;
  	NodeIterator(const Graph* graph, size_type index) : graph_(graph), index_(index) {
	}
  };

  /** Returns a node_iterator pointing to the beginning of the node list
   */
  NodeIterator node_begin() const {
  	return NodeIterator(this, 0);
  }

  /** Returns a node_iterator pointing to the end of the node list
   */
  NodeIterator node_end() const {
  	return NodeIterator(this, num_nodes_);
  }

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

	/** Returns the edge pointed to by the iterator
	 */
	Edge operator*() const {
		assert(uid_ < graph_->size());
		return *it_;
	}
	
	/** Returns an edge iterator that points to the next edge in the graph
	 */
	EdgeIterator& operator++() {
		++it_;
		fix();
	}

	/* Returns true if this iterator points to the same edge, false otherwise
	 */
	bool operator==(const EdgeIterator& it) const {
		return uid_ == it.uid_;
	}

   private:
    friend class Graph;
	const Graph* graph_;
	uid_type uid_;
	IncidentIterator it_;
	EdgeIterator(const Graph* graph, uid_type i) : graph_(graph), uid_(i) {
	}

	/** Private function to maintain representation invariants */
	void fix() {
		Node incident_node = Node(graph_, uid_);
		if (it_ == incident_node.edge_end()) {
			++uid_;
		}
		if (uid_ < graph_->size()) {
			it_ = IncidentIterator(graph_, uid_);
		}
	}
  };

  /** Returns an iterator to the first edge in the graph
   */
  EdgeIterator edge_begin() const {
  	return EdgeIterator(this, 0);
  }
  
  /** Return an iterator to one past the last edge in the graph
   */
  EdgeIterator edge_end() const {
  	return EdgeIterator(this, num_nodes_);
  }

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;
	/** Iterator type */
	typedef std::vector<uid_type>::iterator iterator;

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

	// Return iterator to front
	IncidentIterator& begin() {
		return *this;
	}

	// Return iterator to back
	IncidentIterator& end() {
		it_ = node_.degree();
		return *this;
	}

	/** Return the edge to which the iterator is pointing */
	Edge operator*() const {
		return Edge(graph_, 
					node_.index(), 
					graph_->nodes_[node_.index()].adj[it_]);
	}

	/** Return the iterator to the next element in the indicent list */
	IncidentIterator& operator++() {
		it_++;
		return *this;
	}

	/** Return true if two iterators point to the same element, else false */
	bool operator==(const IncidentIterator& other) const {
		return (node_ == other.node_ && it_ == other.it_);
	}

   private:
    friend class Graph;
	const Graph* graph_;

	// Represenation of node
	Node node_;

	// Index to box that contains uid of second node in edge
	size_type it_;

	IncidentIterator(const Graph* graph, uid_type uid) {
		graph_ = graph;
		node_ = Node(graph, uid);
		it_ = 0;
	}
  };


 private:
	typedef struct node_data {
		Point p;
		node_value_type v;
		size_type degree;
		std::vector<uid_type> adj;
	} node_data;

	size_type num_nodes_ = 0;
	size_type num_edges_ = 0;
 	std::vector<node_data> nodes_;
};
#endif
