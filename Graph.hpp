#ifndef CS207_GRAPH_HPP
#define CS207_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <iostream>
#include <algorithm>
#include <vector>
#include <list>
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
template <typename V, typename E>
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
  /** Define the value type for edge values */
  typedef E edge_value_type;

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

  /** Define edge index type */
  typedef int edge_idx_type;

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
	uid_type uid;
  	idx_type idx;
	Point p_orig;
	mutable Point p;
	mutable node_value_type v;
  } node_data;

  typedef struct edge_data {
  	uid_type uid;
	edge_value_type v;
  } edge_data;

  typedef struct adjacency_data {
	uid_type uid;
	mutable std::list<edge_data> adj_list;
  } adjacency_data;

  typedef typename std::list<edge_data>::iterator adj_list_iterator;

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
  	num_nodes_ = 0;
	num_edges_ = 0;
  }

  /////////////////
  // GRAPH NODES //
  /////////////////

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node : private totally_ordered<Node> {
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
		return (uid_ == x.uid_);
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
		graph_->edges_[ uid_ ].adj_list.size();
	}

   	/** Returns an iterator to beginning of incident iterator list */
	IncidentIterator edge_begin() const {
		IncidentIterator it (graph_, uid_);
		return it;
	}

	/** Returns an iterator to the end of incident iterator list */
	IncidentIterator edge_end() const {
		IncidentIterator it (graph_, uid_);
		return it.to_end_();
	}

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
	const Graph* graph_;
	uid_type uid_ = 0;
	// Construct a node as just a pointer to the graph and an id number
	Node(const Graph* graph, uid_type uid) : graph_(graph), uid_(uid) {
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
  	Node new_node;
	if ( size() < i2u_vect_.size() ) {
		// First case: there have been deleted nodes and we can reuse uids
		// Set all the data fields of the node data structure in nodes
		uid_type reusable_uid = i2u_vect_[size()];
		nodes_[reusable_uid].idx = size();
		nodes_[reusable_uid].p = position;
		nodes_[reusable_uid].p_orig = position;
		nodes_[reusable_uid].v = value;
		new_node = Node(this, reusable_uid);
	}
	else {
		/** Second case: all uids are valid and we must push back nodes_, 
		 * imap_, and indices_ */
		// Set all the data fields of the new node data structure
		node_data temp_node_data;

		temp_node_data.uid = size();
		temp_node_data.idx = size();
		temp_node_data.p = position;
		temp_node_data.p_orig = position;
		temp_node_data.v = value;
		nodes_.push_back(temp_node_data);
		// Push back onto i2u mapping vector
		i2u_vect_.push_back(size());
		// Create an empty adjacency list for this node to put in edges
		adjacency_data adj;
		adj.uid = size();
		edges_.push_back(adj);
		// Create new node off of the given uid
		new_node = Node(this, size());
	}

	++num_nodes_; //new size() = old size() + 1
	return new_node;
  }

  node_iterator remove_node(node_iterator n_it) const {
	Node del_node = *n_it;
	idx_type del_idx = del_node.index(); // this will be the index of the "next"
  	remove_node(*n_it); // Remove the node at n_it
	return NodeIterator(this, i2u_(del_idx)); // return iterator to next el.
  }

  size_type remove_node(const Node& n) {
	idx_type n_idx, prev_idx;
	for( n_idx = n.index() + 1, prev_idx = n.index(); 
		 	n_idx < size(); ++n_idx, ++prev_idx) {
		swap_(prev_idx, n_idx); // After first swap node is invalid
	}
	--num_nodes_;
	return num_nodes_;
  }

  void swap_(idx_type a, idx_type b) {
  	uid_type uid_a = i2u_(a);
	uid_type uid_b = i2u_(b);
	i2u_vect_[a] = uid_b;
	i2u_vect_[b] = uid_a;
	nodes_[uid_a].idx = b;
	nodes_[uid_b].idx = a;
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
  class Edge : private totally_ordered<Edge> {
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

	edge_value_type& value() const {
		Node this_node (graph_, uid1_);
		Edge this_edge (graph_, uid1_, uid2_);
		for(auto it = this_node.edge_begin();it != this_node.edge_end();++it) { 
			if( this_edge == (*it) )
				return it.value_();
		}
		assert(false); // edge is invalid 
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
		if (node1() < node2()) {
			if (x.node1() < x.node2())
				result = (node1() < x.node1() || node2() < x.node2());
			else
				result = (node1() < x.node2() || node2() < x.node1());
		}
		else {
			if (x.node1() < x.node2())
				result = (node2() < x.node1() || node1() < x.node2());
			else
				result = (node2() < x.node2() || node1() < x.node1());
		}
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
  	Edge(const Graph* graph, Node node1, Node node2) {
		graph_ = graph;
		uid1_ = node1.uid_;
		uid2_ = node2.uid_;
	}

  	Edge(const Graph* graph, uid_type uid1, uid_type uid2) : 
								graph_(graph), uid1_(uid1), uid2_(uid2) {
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
	assert(a.index() != b.index());
	// Compute the corresponding uid's of a and b
	uid_type uid_a = i2u_( a.index() );
	uid_type uid_b = i2u_( b.index() );
	if ( !has_edge(a, b) ) {
		// This is the case where a and b have more than 0 edges
		// Insert a and b into each others adjacency lists
		edge_data edge_a_data;
		edge_data edge_b_data;
		edge_a_data.uid = uid_a;
		edge_b_data.uid = uid_b;
		edges_[uid_a].adj_list.push_back( edge_b_data );
		edges_[uid_b].adj_list.push_back( edge_a_data );
		num_edges_++;
	}
	return Edge(this, a, b);
  }

  edge_iterator remove_edge(edge_iterator e_it) const {
	Edge e = *e_it;
	uid_type a = u2i_(e.node1());
	uid_type b = u2i_(e.node1());
  	adj_list_iterator it = remove_edge_( a, b );
	return EdgeIterator(this, a, it);
  }

  size_type remove_edge(const Edge& e) {
  	return remove_edge(e.node1(), e.node2());
  }

  // Remove the nodes' uids from each other's adjacency lists
  size_type remove_edge(const Node& a, const Node& b) {
  	uid_type uid_a = u2i_( a.index() );
  	uid_type uid_b = u2i_( b.index() );
	remove_edge_(uid_a, uid_b);
	return 0;
  }

  adj_list_iterator remove_edge_(uid_type uid_a, uid_type uid_b) const {
	adj_list_iterator it_a, it_b;
	for(it_a = edges_[uid_a].adj_list.begin(); 
								it_a != edges_[uid_a].adj_list.end(); ++it_a) {
		if((*it_a).uid == uid_b) {
			edges_[uid_a].adj_list.erase(it_a);
		}
	}
	for(it_b = edges_[uid_b].adj_list.begin(); 
								it_b != edges_[uid_b].adj_list.end(); ++it_b) {
		if((*it_b).uid == uid_a) {
			edges_[uid_a].adj_list.erase(it_b);
		}
	}
	return ((uid_a < uid_b) ? it_a : it_b);
  }

  /** Check to see if there is an edge between nodes a and b
   * @a a and @a b are distinct valid nodes in the graph
   * @return true if it exists, else return the number of edges
   */
   //TODO: make sure that we don't try to access invalid index of edges
   bool has_edge(const Node& a, const Node& b) {
   	uid_type uid_a = i2u_( a.index() );
	uid_type uid_b = i2u_( b.index() );
	/** RI: a must be in the adjacency list of b and vice versa
	 * so just check one case: is a in the adjacency list of b */
	 for (auto i = edges_[uid_a].adj_list.begin(); 
	 	  i != edges_[uid_a].adj_list.end(); ++i) {
	 	if( (*i).uid == uid_b )
			return true;
	}
	return false; // Returns false if we get to end of list and don't find edge
   }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
	assert(i < size());
	edge_iterator it = edge_begin();
	for (size_type counter = 0; counter < i; counter++) {
		++it;
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
	Node operator* () const {
		return Node(graph_, uid_);
	}

	/** Returns the node iterator that points to the next node in the graph
	 */
	NodeIterator& operator++() {
		idx_type this_index = graph_->u2i_(uid_);
		++this_index;
		if(this_index < graph_->size())
			uid_ = graph_->i2u_(this_index);
		else 
			uid_ = -1; // The end of the graph is represented by a uid -1
		//std::cout << "The uid of this node is" << uid_ << std::endl;
		//std::cout << "The index of this node is" << this_index << std::endl;
		return *this;
	}

	/** Returns true if this iterator points to the same element as the 
	 * parameterized iterator
	 * 
	 * @param[in] @a it is another node iterator
	 * @returns true if the two iterators point to the same node
	 */
	bool operator==(const NodeIterator& it) const {
		return uid_ == it.uid_;
	}

   private:
    friend class Graph;
	const Graph* graph_;
	uid_type uid_;
  	NodeIterator(const Graph* graph, idx_type index) {
		graph_ = graph;
		if(index < graph_->size())
			uid_ = graph_->i2u_(index);
		else 
			uid_ = -1; // The end of graph is represented by uid of -1
	}
  };

  /** Returns a node_iterator pointing to the beginning of the node list
   */
  NodeIterator node_begin() const {
	// Return an node iterator into the graph at the specified index
  	return NodeIterator(this, 0);
  }

  /** Returns a node_iterator pointing to the end of the node list
   */
  NodeIterator node_end() const {
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
	/** Type of the inner iterator */
	typedef IncidentIterator inner_it_type;
	/** Type of outer iterator */
	typedef NodeIterator outer_it_type;

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

	/** Returns the edge pointed to by the iterator
	 */
	Edge operator*() const {
		Node node1 = *outer_pos_;
		Node node2 = (*inner_pos_).node2(); // node 2 is the adjacent node
		//std::cout << "This edge is: " << "(" << node1.index() << "," << node2.index() << ")" << std::endl;
		return Edge(graph_, node1, node2);
	}
	
	/** Returns an edge iterator that points to the next edge in the graph
	 */
	EdgeIterator& operator++() {
		++inner_pos_;
		fix();
		return (*this);
	}

	/* Returns true if this iterator points to the same edge, false otherwise
	 */
	bool operator==(const EdgeIterator& it) const {
		return ((outer_pos_ == graph_->node_end() 
			&& it.outer_pos_ == graph_->node_end()) || 
				(inner_pos_ == it.inner_pos_ && outer_pos_ == it.outer_pos_));
	}

   private:
    friend class Graph;
	const Graph* graph_;
	/** Define node iterator to loop over all nodes and edge iterator to 
	 * iterate over adjacent nodes that form edges */
	outer_it_type outer_pos_;
	inner_it_type inner_pos_;
	/** Initialize an edge iterator to point to the first valid edge */
	EdgeIterator(const Graph* graph) {
		graph_ = graph;
		outer_pos_ = graph_->node_begin();
		inner_pos_ = (*outer_pos_).edge_begin();
		fix();
	}

	EdgeIterator(const Graph* graph, uid_type uid, adj_list_iterator adj_it) {
		graph_ = graph;
		outer_pos_ = graph_->node_begin() + graph_->u2i_(uid);
		inner_pos_ = adj_it;
	}

	EdgeIterator& to_end_() {
		outer_pos_ = graph_->node_end();	
		return *this;
	}
	/** Private function to maintain representation invariants 
	  * @pre All iterator_category iterators must point at valid node/edge
	  * or the end of the vector
	  * @post Iterator category iterator must point at a valid edge unless
	  * it is pointing at the very end of the edge list
	  */
	void fix() {
		while(!outer_end_()) {
			if(!inner_end_() && !outer_end_()) {
				if((*outer_pos_) < (*inner_pos_).node2())
					break; // both inner and outer are valid
				else
					++inner_pos_;
			}
			else if(inner_end_() && !outer_end_()) {
				++outer_pos_;
				if(!outer_end_())
					inner_pos_ = (*outer_pos_).edge_begin();	
			}
		}
	}

	/** Helper functions for deciding when iterator is at end of adj list 
	 * or at the end of a file */
	bool inner_end_() {
		return inner_pos_ == (*outer_pos_).edge_end();
	}
	bool outer_end_() {
		return outer_pos_ == graph_->node_end();
	}
  };

  /** Returns an iterator to the first edge in the graph
   */
  EdgeIterator edge_begin() const {
  	return EdgeIterator(this);
  }
  
  /** Return an iterator to one past the last edge in the graph
   */
  EdgeIterator edge_end() const {
  	return EdgeIterator(this).to_end_();
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
	/** Adjacency list iterator type */
	typedef typename std::list<edge_data>::iterator list_it_type;

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

	/** Return the edge to which the iterator is pointing */
	Edge operator*() const {
		assert(pos_ != graph_->edges_[uid_].adj_list.end());
		return Edge(graph_, uid_, (*pos_).uid);
	}

	/** Return the iterator to the next element in the indicent list */
	IncidentIterator& operator++() {
		++pos_;
		return *this;
	}

	/** Return true if two iterators point to the same element, else false */
	bool operator==(const IncidentIterator& other) const {
		return (graph_ == other.graph_ && 
				uid_ == other.uid_ && pos_ == other.pos_);
	}

	edge_value_type& value_() const {
		return (*pos_).v;
	}

   private:
    friend class Graph;
	const Graph* graph_;

	// UID of this node
	uid_type uid_;
	list_it_type pos_;

	IncidentIterator(const Graph* graph, uid_type uid) {
		// Set private variables
		graph_ = graph;
		uid_ = uid;
		pos_ = graph_->edges_[ uid_ ].adj_list.begin();
	}

	IncidentIterator& to_end_() {
		pos_ = graph_->edges_[ uid_ ].adj_list.end();
		return *this;
	}

  };

 private:
	// Utility functions that maps an indices to uids and vice versa
	uid_type i2u_(idx_type index) const { return i2u_vect_[index];}
	idx_type u2i_(uid_type u) const {return nodes_[u].idx;}
	// Keep track of the number of nodes and edges in the graph
	size_type num_nodes_ = 0;
	size_type num_edges_ = 0;
	// Stores the mapping between uids and indices
	std::vector<uid_type> i2u_vect_;
	// Maps between the uid's and imap nodes
 	std::vector<node_data> nodes_;
	// Maps from indices to imap nodes
	std::vector<imap_idx_type> indices_;
	// Maps from uids to sets of uids that represent adjacency lists
	std::vector<adjacency_data> edges_;
};
#endif
