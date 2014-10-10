#ifndef CS207_GRAPH_HPP
#define CS207_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <set>
#include <cassert>
#include <iterator>

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


  /** define a scalar type */
  typedef double scalar;

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
  typedef unsigned uid_type;

  /** Define idx type */
  typedef unsigned idx_type;

  /** Define imap index type */
  typedef unsigned imap_idx_type;

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

  /** custom type to hold node data */
  typedef struct node_data {
	size_type degree;
	imap_idx_type imap_idx;
	Point p_orig;
	mutable Point p;
	mutable node_value_type v;
  } node_data;

  typedef struct imap_data {
  	uid_type uid;
	size_type idx;
  } imap_data;

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

	/** Return this nodes position as a modifiable reference */
	Point& position() {
		return graph_->nodes_[uid_].p;
	}

    /** Return this node's index, a number in the range [0, graph_size).
	 * @pre The node represented by this uid must be a valid node
	 * 		If this precondition is not met, the behavior is undefined. 
	 * 		The assert may fail, or the function could return a valid, 
	 *		incorrect index.
	 */
    idx_type index() const {
		return graph_->u2i_(uid_);
    }

    /** Test whether this node and @a x are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& x) const {
		return (index() == x.index());
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
		return (index() < x.index());
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
	// TODO: fix this once IncidentIterator is updated with new representation
	IncidentIterator& edge_begin() const {
		return IncidentIterator(graph_, uid_).begin();
	}

	/** Returns an iterator to the end of incident iterator list */
	// TODO: fix this once IncidentIterator is updated with new representation
	IncidentIterator& edge_end() const {
		return IncidentIterator(graph_, uid_).end();
	}

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
	const Graph* graph_;
	uid_type uid_ = 0;
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
  Node add_node(const Point& position, 
  					const node_value_type& value = node_value_type()) {
	node_data temp_node_data;
	imap_data temp_imap_data;
  	Node new_node;
	if ( size() < imap_.size() ) {
		// First case: there have been deleted nodes and we can reuse uids
		// Set all the data fields of the node data structure in nodes
		uid_type reusable_uid = imap_[size()].uid;
		nodes_[reusable_uid].p = position;
		nodes_[reusable_uid].p_orig = position;
		nodes_[reusable_uid].v = value;
		nodes_[reusable_uid].degree = 0;
		// Update the index of the node
		imap_[size()].idx = size();
		// Update the index of the indices vector point back to this imap entry
		indices_[size()] = size();
		// Make new node using the reusable uid
		new_node = Node(this, reusable_uid);
	}
	else {
		/** Second case: all uids are valid and we must push back nodes_, 
		 * imap_, and indices_ */
		// Set all the data fields of the new node data structure
		temp_node_data.p = position;
		temp_node_data.p_orig = position;
		temp_node_data.v = value;
		temp_node_data.degree = 0;
		nodes_.push_back(temp_node_data);
		// Set all the data fields of the new imap data structure
		temp_imap_data.uid = size();
		temp_imap_data.idx = size();
		imap_.push_back(temp_imap_data);
		// Push the new imap index onto the back indices
		indices_.push_back(size());
		// Create new node off of the given uid
		new_node = Node(this, size());
	}

	++num_nodes_; //new size() = old size() + 1
	return new_node;
  }

  /** Determine if this Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
  	return n.index() < size();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(idx_type i) const {
	assert( i < size() );
    return Node(this, i2u_(i));
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

	scalar length() const {
		return distance(graph_->nodes_[ uid1_ ].p_orig, 
								graph_->nodes_[ uid2_ ].p_orig);
	}
    /** Test whether this edge and @a x are equal.
     * @pre RI for Edges must hold: uid1_ < uid2_
     * Equal edges are from the same graph and have the same nodes.
     */
    bool operator==(const Edge& x) const {
		return ((node1() == x.node1() && node2() == x.node2()) ||
				(node2() == x.node1() && node1() == x.node2()));
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
		bool result; 
		if (node1() < node2())
			if (x.node1() < x.node2())
				result = (node1() < x.node1() || node2() < x.node2());
			else
				result = (node1() < x.node2() || node2() < x.node1());
		else
			if (x.node1() < x.node2())
				result = (node2() < x.node1() || node1() < x.node2());
			else
				result = (node2() < x.node2() || node1() < x.node1());
		return result;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
	// data members
	const Graph* graph_;
	uid_type uid1_;
	uid_type uid2_;
	// Constructor available to the Graph class for constructing edges
  	Edge(const Graph* graph, idx_type idx1, idx_type idx2) {
		graph_ = graph;
		uid1_ = graph->i2u_(idx1);
		uid2_ = graph->i2u_(idx2);
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
	assert(a.index() != b.index()); // no self edges
	// Compute the corresponding uid's of a and b
   	uid_type uid_a = i2u_( a.index() );
	uid_type uid_b = i2u_( b.index() );
	// Insert a and b into each others adjacency lists
	edges_[uid_a].insert(uid_b);
	edges_[uid_b].insert(uid_a);
	// Add one to the degree of the edges
	nodes_[uid_a].degree++;
	nodes_[uid_b].degree++;
	// Increment the number of edges by one
	num_edges_++;
	return Edge(this, a.index(), b.index());
  }

  /** Check to see if there is an edge between nodes a and b
   * @a a and @a b are distinct valid nodes in the graph
   * @return true if it exists, else return the number of edges
   */
   bool has_edge(const Node& a, const Node& b) {
   	uid_type uid_a = i2u_( a.index() );
	uid_type uid_b = i2u_( b.index() );
	/** RI: a must be in the adjacency list of b and vice versa
	 * so just check one case: is a in the adjacency list of b */
	return !(edges_[uid_a].find( uid_b ) == edges_[uid_a].end());
   }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
	assert(i < size());
	EdgeIterator it = edge_begin();
	for (size_type counter = 0; counter < i; counter++) {
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
		idx_type this_index = *it_;
		return Node(graph_, graph_->i2u_(this_index));
	}

	/** Returns the node iterator that points to the next node in the graph
	 */
	NodeIterator& operator++() {
		++it_;
		return *this;
	}

	/** Returns true if this iterator points to the same element as the 
	 * parameterized iterator
	 * 
	 * @param[in] @a it is another node iterator
	 * @returns true if the two iterators point to the same node
	 */
	bool operator==(const NodeIterator& it) const {
		idx_type this_index = *it_;
		return (Node(graph_, graph_->i2u_(this_index)) == *it);
	}

   private:
    friend class Graph;
	const Graph* graph_;
	std::vector<idx_type>::iterator it_;
  	NodeIterator(Graph* graph, idx_type index) {
		graph_ = graph;
		it_ = graph->indices_.begin() + index;
	}
  };

  /** Returns a node_iterator pointing to the beginning of the node list
   */
  NodeIterator node_begin() {
	// Return an node iterator into the graph at the specified index
  	return NodeIterator(this, 0);
  }

  /** Returns a node_iterator pointing to the end of the node list
   */
  NodeIterator node_end() {
	// Return an node iterator into the graph at the specified index
  	return NodeIterator(this, size());
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
	 * @pre this->it_nodes_ < nodes_end_
	 */
	Edge operator*() const {
		return Edge(graph_, *outer_, graph_->u2i_( *inner_ ));
	}
	
	/** Returns an edge iterator that points to the next edge in the graph
	 */
	EdgeIterator& operator++() {
		next();
		fix();
		return *this;
	}

	/* Returns true if this iterator points to the same edge, false otherwise
	 */
	bool operator==(const EdgeIterator& it) const {
		return (outer_ == it.outer_ && inner_ == it.inner_);
	}

   private:
    friend class Graph;
	const Graph* graph_;
	/** Define node iterator to loop over all nodes and edge iterator to 
	 * iterate over adjacent nodes that form edges */
	std::vector<idx_type>::iterator outer_;
	std::set<uid_type>::iterator inner_;
	std::vector<idx_type>::iterator end_;
	/** Initialize an edge iterator to point to the first valid edge */
	EdgeIterator(Graph* graph, idx_type index) {
		graph_ = graph;
		outer_ = graph->indices_.begin() + index;
		end_ = graph->indices_.begin() + graph->size();
		if (outer_ != end_) 
			inner_ = graph->edges_[ graph_->i2u_(index) ].begin();
		fix();
	}

	/** Private function to maintain representation invariants 
	  * @pre All iterator_category iterators must point at valid node/edge
	  * or the end of the vector
	  * @post Iterator category iterator must point at a valid edge unless
	  * it is pointing at the very end of the edge list
	  */
	void fix() {
		while ((outer_ != end_) && (graph_->imap_[ *outer_ ].uid > *inner_)) {
			next();
		}
	}

	void next() {
		++inner_;
		uid_type this_uid = graph_->imap_[ *outer_ ].uid;
		if (inner_ == graph_->edges_[ this_uid ].end()) {
			++outer_;
			uid_type next_uid = graph_->imap_[ *outer_ ].uid;
			inner_ = graph_->edges_[ next_uid ].begin();
		}
	}
  };

  /** Returns an iterator to the first edge in the graph
   */
  EdgeIterator edge_begin() {
  	return EdgeIterator(this, 0);
  }
  
  /** Return an iterator to one past the last edge in the graph
   */
  EdgeIterator edge_end() {
  	return EdgeIterator(this, size());
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
		it_ = graph_->edges_[ uid_ ].end();
		return *this;
	}

	/** Return the edge to which the iterator is pointing */
	Edge operator*() const {
		return Edge(graph_, 
					graph_->u2i_( uid_ ), 
					graph_->u2i_( *it_ ));
	}

	/** Return the iterator to the next element in the indicent list */
	IncidentIterator& operator++() {
		++it_;
		return *this;
	}

	/** Return true if two iterators point to the same element, else false */
	bool operator==(const IncidentIterator& other) const {
		return (graph_ == other.graph_ && 
				uid_ == other.uid_ && it_ == other.it_);
	}

   private:
    friend class Graph;
	const Graph* graph_;

	// UID of this node
	uid_type uid_;

	// Index to box that contains uid of second node in edge
	std::set<uid_type>::iterator it_;

	IncidentIterator(const Graph* graph, uid_type uid) {
		// Set private variables
		graph_ = graph;
		uid_ = uid;
		// Set the iterator to the beginning of the adjacency list by default
		it_ = graph->edges_[ uid ].begin();
	}
  };


 private:
	// Utility functions that maps an indices to uids and vice versa
	uid_type i2u_(idx_type index) const {
		return imap_[ indices_[index] ].uid;
	}
	idx_type u2i_(uid_type u) const {
		return imap_[ nodes_[u].imap_idx ].idx;
	}
	// Keep track of the number of nodes and edges in the graph
	size_type num_nodes_ = 0;
	size_type num_edges_ = 0;
	// Stores the mapping between uids and indices
	std::vector<imap_data> imap_;
	// Maps between the uid's and imap nodes
 	std::vector<node_data> nodes_;
	// Maps from indices to imap nodes
	std::vector<imap_idx_type> indices_;
	// Maps from uids to sets of uids that represent adjacency lists
	std::vector<std::set<uid_type>> edges_;
};
#endif
