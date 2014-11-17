#ifndef CS207_GRAPH_HPP
#define CS207_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <tuple>

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
 public:

  /////////////////////////////
  // PUBLIC TYPE DEFINITIONS //
  /////////////////////////////

  //
  typedef V node_value_type;
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

  /** Package for Node's data*/
   struct NodeData : private equality_comparable<NodeData>{
 	  size_type uid;
 	  Point position;
 	  node_value_type value;

 	 /** provide comparison to other structs and ints */
	  bool operator==(const NodeData& x) const {
		return uid == x.uid;
	  }
	  bool operator==(size_type i) const {
		  return uid == i;
	  }
   };

  /** Container type for Adjacency list */
  struct EdgeData : private equality_comparable<EdgeData>{
	  size_type n2_id;
	  edge_value_type value;

	  /** provide comparison to other structs and ints */
	  bool operator==(const EdgeData& x) const {
		return n2_id == x.n2_id;
	  }
	  bool operator==(size_type i) const {
		  return n2_id == i;
	  }
  };
  typedef std::vector<EdgeData> adj_list__subtype;
  typedef std::vector<adj_list__subtype> adj_list__type;

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
    edge_count_ = 0;
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
    return nodes_.size();
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    adj_list_.clear();
    uid2idx_.clear();
    available_uids_.clear();
    edge_count_ = 0;

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

    /** Return this node's position.
     * O(1)
     * */
    const Point& position() const {
    	return graph_->node_data_by_uid(uid_).position;
    }


    /** Return non const node position.
     * O(1)
     */
	Point& position() {
		return graph_->node_data_by_uid(uid_).position;
	}

    /** Return this node's index, a number in the range [0, graph_size).
     * O(1)
     * */
    size_type index() const {
    	return graph_->uid2idx_[uid_];
    }

    /** Test whether this node and @a x are equal.
     *
     * Equal nodes have the same graph and the same index.
     * @return true iff the nodes have the same graph and same uid
     */
    bool operator==(const Node& x) const {
      return (graph_ == x.graph_) && (uid_ == x.uid_);
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
    	return (graph_ < x.graph_ || (graph_ == x.graph_ &&  uid_ < x.uid_ ));
    }

    /** Returns how many edges are connected to this node
     * @pre this is a valid node object
     * @return the number of edges incident to this node
     * O(1)
     */
    size_type degree() const {
    	return graph_->adj_list_[uid_].size();
    }

    /** returns an incident iterator to the beginning of this node's edges
     * O(1)
     */
    incident_iterator edge_begin() const {
    	return IncidentIterator(graph_, uid_, 0);
    }

    /** returns an iterator pointing right past the last edge for this node
     * O(1)
     */
    incident_iterator edge_end() const {
    	return IncidentIterator(graph_, uid_, degree());
    }

    /** Return the value of this node
     * O(1)
     */
    node_value_type& value() {
    	return graph_->node_data_by_uid(uid_).value;
    }

    /** Return the value of this node */
    const node_value_type& value() const {
    	return graph_->node_data_by_uid(uid_).value;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // data members
    Graph* graph_;
    size_type uid_; //no delete function, so uids will remain in order

    // private constructor for Graph
    Node(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
      assert(graph_ != nullptr);
      if (!graph_->valid_node_uid(uid_)) {
    	  std::cerr << "invalid node uid" << std::endl;
    	  graph_->debug_uid(uid_);
      }
      assert(graph_->valid_node_uid(uid_));
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
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
	  size_type next_uid = take_next_uid();
	  size_type next_idx = num_nodes();
	  NodeData nd;
	  nd.position = position;
	  nd.value = val;
	  nd.uid = next_uid;
	  // recording
	  nodes_.push_back(nd);
	  uid2idx_[next_uid] = next_idx;


	  // adding new row to adjecency list
	  adj_list_.push_back(adj_list__subtype());

	  return Node(this,next_uid);
  }

  /** Determine if this Node belongs to this Graph
   * @return True iff @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return n.graph_ == this;
  }

  /** Determine if this is a valid index for a node
   *
   * @param i possible index for a node
   * @return True if there exists a node in this graph with index i
   *
   * complexity: O(1)
   */
  bool valid_node_index(const size_type i) const {
	  return (i >= 0 && i < num_nodes());
  }



  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // just a reference
	//assert(valid_node_uid(i));
    return Node(this,nodes_[i].uid);
  }




  /** Removes a node from the graph
   *
   * @param[in] n Node to be removed
   * @return 1 if old has_node(@a n), else 0
   *
   * @post new size() == old_size() - result
   * @post new num_edges() = old num_edges() - n.degree()
   *
   * complexity: O(d)
   *
   * Can invalidate iterators
   * if old has_node(@a n), @a n becomes invalid.
   * All other nodes remain valid, because they are referenced
   * by static uids, rather than indices
   */
  size_type remove_node(const Node& n) {
	  size_type degree = n.degree();
	  size_type oldn = num_nodes();
	  size_type olde = num_edges();
	  size_type uid = n.uid_;


	  if (!valid_node_uid(uid)) {
		  // this is not a good node
		  std::cerr << "not a node anyway" << std::endl;
		  return 0;
	  }
	  size_type idx = n.index();

	  // deleting edges
	  // deleting from other adj_list spots first
	  for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
		  // finding it in this list
		  size_type uid2 = (*it).node2().uid_;
		  auto it2 = std::find(adj_list_[uid2].begin(), adj_list_[uid2].end(), uid);

		  adj_list_[uid2].erase(it2);
	  }

	  // updating edge count
	  edge_count_ -= degree;

	  // deleting this row
	  adj_list_[uid].clear();

	  // copy last node to this spot
	  nodes_[idx] = nodes_.back();
	  uid2idx_[nodes_.back().uid] = idx;
	  nodes_.pop_back();

	  // updating look up tables
	  uid2idx_[uid] = -1;

	  available_uids_.push_back(uid);

	  assert(num_nodes() == oldn - 1);
	  assert(num_edges() == olde - degree);

	  return 1;
  }

  /** See function above for description */
  node_iterator remove_node(node_iterator nit) {
	  assert(valid_node_index(nit.idx_));
	  remove_node(*nit);
	  return nit;
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
      return graph_->node(n1_uid_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node(n2_uid_);
    }

    /** return the value type of this edge */
	edge_value_type& value() {
		std::cout << "Finding a value in the edge" << std::endl;
		auto it = std::find(graph_->adj_list_[n1_uid_].begin(),
				graph_->adj_list_[n1_uid_].end(), n2_uid_);
		return (*it).value;
	}

	/** Return the value of this edge */
	const edge_value_type& value() const {
		auto it = std::find(graph_->adj_list_[n1_uid_].begin(),
				graph_->adj_list_[n1_uid_].end(), n2_uid_);
		return (*it).value;
	}

	/** Set the value of this edge */
	void set_value(edge_value_type v) {
		auto it = std::find(graph_->adj_list_[n1_uid_].begin(),
				graph_->adj_list_[n1_uid_].end(), n2_uid_);
		(*it).value = v;
		it = std::find(graph_->adj_list_[n2_uid_].begin(),
				graph_->adj_list_[n2_uid_].end(), n1_uid_);
		(*it).value = v;
	}

    /** returns the vector between node1 and node2; */
    Point difference() const {
    	return node1().position() - node2().position();
    }

    /** return the length of this edge */
    double length() const {
    	return norm(difference());
    }

    /** Test whether this edge and @a x are equal.
     *
     * Equal edges are from the same graph and have the same nodes.
     */
    bool operator==(const Edge& x) const {
      bool nodes_same = node1() == x.node1() && node2() == x.node2();
      bool nodes_flipped = node1() == x.node2() && node2() == x.node1();
      bool nodes_equivalent = nodes_same || nodes_flipped;
      bool graphs_same = graph_ == x.graph_;
      return graphs_same && nodes_equivalent;
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
      // first handle graphs
      if (graph_ < x.graph_)
        return true;
      else if (graph_ > x.graph_)
        return false;

      size_type big, sm, x_big, x_sm;
      big = std::max(n1_uid_,n2_uid_);
      sm = std::min(n1_uid_,n2_uid_);
      x_big = std::max(x.n1_uid_,x.n2_uid_);
      x_sm = std::min(x.n1_uid_,x.n2_uid_);

      return std::tie(sm, big) < std::tie(x_sm, x_big);

      // graphs are equal, fall back to node 1
      if (sm < x_sm)
        return true;
      else if (sm > x_sm)
        return false;

      // fall back to node 2
      if (big < x_big)
        return true;
      else if (big > x_big)
        return false;

      // they are equal
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // data members
    Graph* graph_;

    // node uids
    size_type n1_uid_;
    size_type n2_uid_;

    /** private constructor for Graph
     * @pre n1_id, n2_id are both >= 0 and < num_edges
     * @pre graph_ is a pointer to valid graph object
     * @pre graph_ already has this edge
     */
    Edge(const Graph* graph, size_type n1_id, size_type n2_id)
        : graph_(const_cast<Graph*>(graph)), n1_uid_(n1_id), n2_uid_(n2_id) {

      assert(graph_ != nullptr);
      assert(graph_->valid_node_uid(n1_uid_));
      assert(graph_->valid_node_uid(n2_uid_));
      assert(graph_->has_edge(n1_uid_, n2_uid_));
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_count_;
  }



  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post has_edge(@a b, @a a) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: O(deg(a))
   */
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& val = edge_value_type()) {
    // checking pre
    assert(a.graph_ == this);
    assert(b.graph_ == this);
    assert(has_node(a));
    assert(has_node(b));
    assert(a.uid_ != b.uid_);

    // check if node already exists
    if (has_edge(a.uid_,b.uid_)) {
      return edge(a.uid_, b.uid_);
    }

    size_type old = num_edges();

    // putting into adjacency list
    EdgeData ed_b, ed_a;
    ed_b.n2_id = b.uid_;
    ed_a.n2_id = a.uid_;
    ed_b.value = val;
    ed_a.value = val;

    adj_list_[a.uid_].push_back(ed_b);
    adj_list_[b.uid_].push_back(ed_a);
    ++edge_count_;

    // creating proxy to return
    Edge e = edge(a.uid_, b.uid_);

    // checking post
    assert(has_edge(a.uid_,b.uid_));
    assert(has_edge(b.uid_,a.uid_));
    assert(num_edges() == old + 1);

    return e;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid indices of nodes of this graph
   * @pre all indices between 0 and size() are valid nodes
   * @return true if, for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: O(deg(ai))
   */
  bool has_edge(const size_type uid1, const size_type uid2) const
  {
	assert(valid_node_uid(uid1));
	assert(valid_node_uid(uid2));

    // at least one must be in map
	auto it = std::find(adj_list_[uid1].begin(), adj_list_[uid1].end(), uid2);
	return it != adj_list_[uid1].end();
  }

  bool has_edge(const Node& n1, const Node& n2) const {
	  return has_edge(n1.uid_, n2.uid_);
  }

  /** Return the edge with nodes i and j
   * @pre i, j are valid node uids
   * @pre edge(i,j) already exists
   * Complexity: O(d(i) + d(j)) for assertions, O(1) otherwise
   */
  Edge edge(size_type uid1, size_type uid2) const {
	assert(valid_node_uid(uid1));
    assert(valid_node_uid(uid2));
    assert(has_edge(uid1,uid2));
    assert(has_edge(uid2,uid1));
    return Edge(this,uid1,uid2);
  }

  Edge edge(const Node& n1, const Node& n2) const {
	  // TODO: figure out the invariant for ordering of UIDs in edge construct
	  // This is the only invariant that gives the correct behavior when
	  // Accessing the value of the edge and modifying the data in that 
	  // edge data structure
	  if(n1.uid_ < n2.uid_) return Edge(this, n2.uid_, n1.uid_);
	  else return Edge(this, n1.uid_, n2.uid_);
  }

  /** Removes an edge from the graph
   * @param e The edge object that should be removed
   * complexity: O(d)
   * does not effect any RIs or invalidate any objects
   * @return the number of edges removed from the graph
   */
  size_type remove_edge(const Edge& e) {
	  return remove_edge(e.node1(), e.node2());
  }

  size_type remove_edge(const Node& n1, const Node& n2) {
	  if (!has_edge(n1.uid_, n2.uid_))
		  return 0;

	  auto it1 = std::find(adj_list_[n1.uid_].begin(), adj_list_[n1.uid_].end(), n2.uid_);
	  adj_list_[n1.uid_].erase(it1);

	  auto it2 = std::find(adj_list_[n2.uid_].begin(), adj_list_[n2.uid_].end(), n1.uid_);
	  adj_list_[n2.uid_].erase(it2);

	  --edge_count_;

	  return 1;
  }

  /* prints a graph in adj_list_ rep.
   * Utility function
   */
  void print_graph() const {
	  std::cerr << "Graph: {" << std::endl;
	  for (size_type i = 0; i < adj_list_.size(); i++)
	  {
		  std::cerr << " Node (" << i << "): ";
		  for (size_type j = 0; j < adj_list_[i].size(); j++)
		  {
			  std::cerr << adj_list_[i][j].n2_id << ", ";
		  }
		  std::cerr << std::endl;
	  }
	  std::cerr << "}" << std::endl << std::endl;

	  debug_uid(0);
  }



  ///////////////
  // Iterators //
  ///////////////

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator.
   * @invariant uid is always either a valid node index, or represents node_end() in which case uid = num_nodes()
   */
  class NodeIterator : private equality_comparable<NodeIterator> {
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

    /** Dereferences an iterator
     * @pre uid_ is a valid node index
     * */
    Node operator*() const {
    	assert(graph_->valid_node_index(idx_));
    	return graph_->node(idx_);
    }

    /** Increments UID and returns this
     * @pre uid_ is currently a valid node
     * */
    NodeIterator& operator++() {
    	assert(graph_->valid_node_index(idx_));
    	++idx_;
    	return *this;
    }

    /** Compare by graph and uid_
     * can't pass on comparison to the edge itself
     * because that would potentially dereference
     * invalid iterators
     */
    bool operator==(const NodeIterator& it) const {
      return (graph_ == it.graph_) && (idx_ == it.idx_);
    }

   private:
    friend class Graph;

    /** Private constructor
     * @pre uid is either a valid node uid, or points just to the end of the vector
     */
    NodeIterator(const Graph* graph, size_type idx)
        : graph_(const_cast<Graph*>(graph)), idx_(idx) {
    	assert(idx_ >= 0 && idx_ <= graph_->num_nodes());
    }

    // data members
    Graph* graph_;
    size_type idx_;
  };



  /** returns iterator to node with uid 0 */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** returns iterator to node with uid num_nodes() */
  node_iterator node_end() const {
    return NodeIterator(this, num_nodes());
  }


  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator.
   * @invariant n1_ind and n2_ind must always correspond to an element in adj_list_[n1_ind][n2_ind]
   * 	unless this represents edge_end(), in which case n1_ind = num_nodes() and n2_ind = 0
   * */
  class EdgeIterator : private equality_comparable<EdgeIterator>{
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


    /** Dereferences an iterator
     * @pre this edge exists in adj_list_
     * O(1)
     * */
    Edge operator*() const {
    	assert(n1_ind_ < graph_->adj_list_.size());
		assert(n2_ind_ < graph_->adj_list_[n1_ind_].size());

		size_type n1_uid = n1_ind_;
		size_type n2_uid = graph_->adj_list_[n1_ind_][n2_ind_].n2_id;
		return graph_->edge(n1_uid, n2_uid);
    }

    /** Goes to next forward edge and returns this
     * @pre 0 <= n1_ind < adj_list_.size()
     * @pre 0 <= n2_ind < adj_list_[n1_ind].size()
     *
     * Allow n1 to go one beyond for iterator.end()
     * n2 must still be in bounds
     *
     * @post (0 <= n1_ind < adj_list_.size() && 0 <= n2_ind < adj_list_[n1_ind].size()) or
     * 		(n1_ind = adj_list_.size() && n2_ind == 0 // .end() condition
     * */
    EdgeIterator& operator++() {
		assert(n1_ind_ >= 0);
		assert(n2_ind_ >= 0);
		assert(n1_ind_ < graph_->adj_list_.size());
		assert(n2_ind_ < graph_->adj_list_[n1_ind_].size());

		// only accept forward edges i.e. 2 -> 4
		do {
			++n2_ind_;
			// do we need to go to the next list?
			// be careful of nodes without any edges
			while (n2_ind_ == graph_->adj_list_[n1_ind_].size() && n1_ind_ < graph_->adj_list_.size()) {
				n2_ind_ = 0;
				++n1_ind_;
			}
		}
		while (n1_ind_ < graph_->adj_list_.size() && !graph_->forward_edge(n1_ind_, n2_ind_));

		// checking post conditions
		if (n1_ind_ == graph_->adj_list_.size())
		{
			//at end
			assert(n2_ind_ == 0);
		}
		else
		{
			// still valid
			assert(n1_ind_ < graph_->adj_list_.size());
			assert(n2_ind_ < graph_->adj_list_[n1_ind_].size());
		}

		return *this;
    }

    /** Compare by graph and indices
     * can't pass on comparison to the edge itself
     * because that would potentially dereference
     * invalid iterators
     */
    bool operator==(const EdgeIterator& it) const {
    	return (graph_ == it.graph_) && (n1_ind_ == it.n1_ind_) && (n2_ind_ == it.n2_ind_);
    }

    /* prints useful info about iterator
     *
     */
    void debug() const {
    	std::cerr << "Edge iterator: " << std::endl;
    	std::cerr << "n1: (" << n1_ind_ << " / " << graph_->adj_list_.size() << ") ";
    	std::cerr << "n2: (" << n2_ind_ << " / " << graph_->adj_list_[n1_ind_].size() << ")" << std::endl;
    }

   private:
    friend class Graph;

    /** Private constructor for Graph
     *
     * @pre graph_ is a valid graph
     * @pre (0 <= n1_ind < adj_list_.size() && 0 <= n2_ind < adj_list_[n1_ind].size()) or
     * 		(n1_ind = adj_list_.size() && n2_ind == 0 // .end() condition
     */
    EdgeIterator(const Graph* graph, size_type n1_ind, size_type n2_ind)
        : graph_(const_cast<Graph*>(graph)), n1_ind_(n1_ind), n2_ind_(n2_ind) {

    	// checking post conditions
		if (n1_ind_ == graph_->adj_list_.size())
		{
			//at end
			assert(n2_ind_ == 0);
		}
		else
		{
			// still valid
			assert(n1_ind_ < graph_->adj_list_.size());
			assert(n2_ind_ < graph_->adj_list_[n1_ind_].size());
		}
    }

    // data members
    Graph* graph_;
    size_type n1_ind_; // not uids, indices
    size_type n2_ind_; // not uids, indices
  };

  /** returns an iterator with uid 0 0 */
  edge_iterator edge_begin() const {
	  // in case the first node doesn't have an edge
	  size_type n_id = 0;
	  while (adj_list_[n_id].size() == 0)
		  ++n_id;
	  return EdgeIterator(this, n_id, 0);
  }

  /** returns an iterator with num_edges 0 */
  edge_iterator edge_end() const {
    return EdgeIterator(this, adj_list_.size(), 0);
  }


  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator.
   * @invariant n1_ind and n2_ind must always correspond to an element in adj_list_[n1_ind][n2_ind]
   * 	unless this represents edge_end(), in which case n1_ind = num_nodes() and n2_ind = 0
   * */
  class IncidentIterator : private equality_comparable<IncidentIterator> {
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

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

	/** Dereferences an iterator */
	Edge operator*() const {
		assert(n1_ind_ < graph_->adj_list_.size());
		assert(n2_ind_ < graph_->adj_list_[n1_ind_].size());

		size_type n1_uid = n1_ind_;
		size_type n2_uid = graph_->adj_list_[n1_ind_][n2_ind_].n2_id;
		return graph_->edge(n1_uid, n2_uid);
	}

	/** Goes to next forward edge and returns this
	* @pre 0 <= n1_ind < adj_list_.size()
	* @pre 0 <= n2_ind < adj_list_[n1_ind].size()
	*
	* Allow n1 to go one beyond for iterator.end()
	* n2 must still be in bounds
	*
	* @post (0 <= n1_ind < adj_list_.size() && 0 <= n2_ind < adj_list_[n1_ind].size()) or
	* 		(n1_ind = adj_list_.size() && n2_ind == 0
	* */
	IncidentIterator& operator++() {
		assert(n1_ind_ >= 0);
		assert(n2_ind_ >= 0);
		assert(n1_ind_ < graph_->adj_list_.size());
		assert(n2_ind_ < graph_->adj_list_[n1_ind_].size());

		++n2_ind_;

		assert(n1_ind_ < graph_->adj_list_.size());
		assert(n2_ind_ <= graph_->adj_list_[n1_ind_].size());


		return *this;
	}

	/** Compare by graph and indices
	* can't pass on comparison to the edge itself
	* because that would potentially dereference
	* invalid iterators
	*/
	bool operator==(const IncidentIterator& it) const {
		return (graph_ == it.graph_) && (n1_ind_ == it.n1_ind_) && (n2_ind_ == it.n2_ind_);
	}

	/* prints useful info about iterator
	*
	*/
	void debug() const {
		std::cerr << "Edge iterator: " << std::endl;
		std::cerr << "n1: (" << n1_ind_ << " / " << graph_->adj_list_.size() << ") ";
		std::cerr << "n2: (" << n2_ind_ << " / " << graph_->adj_list_[n1_ind_].size() << ")" << std::endl;
	}

   private:
    friend class Graph;
    friend class Node;

    /** private constructor
     *
     * @pre n1_ind is a valid node
     * @pre 0 <= n2_ind <= deg(n1)
     */
    IncidentIterator(const Graph* graph, size_type n1_ind, size_type n2_ind)
    	: graph_(const_cast<Graph*>(graph)), n1_ind_(n1_ind), n2_ind_(n2_ind) {

		assert(n1_ind_ < graph_->adj_list_.size());
		assert(n2_ind_ <= graph_->adj_list_[n1_ind_].size());
    }

    // data members
	Graph* graph_;
	size_type n1_ind_; // not uids, position in vectors
	size_type n2_ind_;

  };

 private:
	friend class EdgeIterator;
	friend class IncidentIterator;
	// Private functions --------------------------------------------------------
	bool valid_node_uid(const size_type i) const {
		assert(i >= 0 && i < uid2idx_.size());
		return valid_node_index(uid2idx_[i]);
	}

	/* pops the next available uid out of the list
	 * or adds space at the end if none have been created
	 * O(1)
	 * @pre any uids that have been removed have been inserted
	 * 	into available_uids_
	 * @post old size(available_uids) == 0 or
	 * 	new size(avail) = old size(avail) - 1
	 * @post uid2idx[uid] == nodes.size()
	 * @return @a uid, a uid that should be assigned to a new node
	 */
	size_type take_next_uid() {
		// just pass next spot in line
		if (available_uids_.size() == 0) {
			uid2idx_.push_back(num_nodes());
			return num_nodes();
		} else {
			size_type uid = available_uids_.back();
			available_uids_.pop_back();
			return uid;
		}
	}

	/** returns the node by uid
	 */
	NodeData& node_data_by_uid(size_type uid) {
		return nodes_[uid2idx_[uid]];
	}

	/** prints debug information about a particular node
	 */
	void debug_uid(size_type uid) const {
		std::cerr << std::endl << "Debugging Node UID" << std::endl;
		std::cerr << "uid: " << uid << std::endl;
		std::cerr << "uid2idx: " << std::endl;
		for (auto it = uid2idx_.begin(); it != uid2idx_.end(); it++) {
			std::cerr << *it << ", ";
		}
		std::cerr << std::endl;

		std::cerr << "nodes.uid: " << std::endl;
		for (auto it = nodes_.begin(); it != nodes_.end(); it++) {
			std::cerr << (*it).uid << ", ";
		}
		std::cerr << std::endl;

		std::cerr << "available uids: " << std::endl;
		for (auto it = available_uids_.begin(); it != available_uids_.end();
				++it) {
			std::cerr << *it << ", ";
		}
		std::cerr << std::endl;
	}

	/** returns whether adj_list_[i][j] is a forward edge
	 * Utility function
	 * checks whether node1 uid > node2 uid
	 */
	bool forward_edge(size_type i, size_type j) const {
		return i > adj_list_[i][j].n2_id;
	}


  // Private data -------------------------------------------------------------


  /** node data sources
   * indices are mutable
   * uids are permanent
   * */
  std::vector<NodeData> nodes_; // sorted by index, not uid
  std::vector<size_type> uid2idx_; // uid2idx_[uid] = idx. -1 = not taken
  std::vector<size_type> available_uids_;

  // edge adjacency list data source
  adj_list__type adj_list_;
  size_type edge_count_;

};

#endif
