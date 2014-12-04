#pragma once
/** @file Mesh.hpp
 * @brief A Mesh is composed of nodes, edges, and triangles such that:
 *  -- All triangles have three nodes and three edges.
 *  -- All edges belong to at least one triangle and at most two triangles.
 */

/** @class Mesh
 * @brief A template for 3D triangular meshes.
 *
 * Users can add triangles and retrieve nodes, edges, and triangles.
 */

#include <fstream>
#include <cmath>
#include <set>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"

#include "Point.hpp"
#include "Graph.hpp"
#include "CollisionDetector.hpp"

template <typename UserNodeData, typename UserEdgeData, typename UserTriData>
class Mesh {
	// HW3: YOUR CODE HERE
	// Write all typedefs, public member functions, private member functions,
	//   inner classes, private member data, etc.
public:
	/** Data members that will be inside the 2 graphs */
	struct triangle_data;
	struct link_data;
	struct vertex_data;
	struct edge_data;
	class Triangle;
	class Vertex;
	class Edge;
	class Link;

	// Iterator Clases
	template <typename IT> class AdjacentIterator;
	template <typename IT, typename RET> class TransformIterator;


	// Primitive types
	typedef unsigned idx_type;

	// inner graphs themselves
	typedef Graph<triangle_data, link_data> TriGraph;
	typedef Graph<vertex_data, edge_data> VertGraph;

	// inner nodes and edges
	typedef typename TriGraph::node_type tri_node;
	typedef typename TriGraph::edge_type tri_edge;
	typedef typename VertGraph::node_type vert_node;
	typedef typename VertGraph::edge_type vert_edge;

	// Iterator types
	typedef typename VertGraph::node_iterator vert_node_iterator;
	typedef typename VertGraph::edge_iterator vert_edge_iterator;
	typedef typename VertGraph::incident_iterator incident_iterator;
	typedef typename TriGraph::node_iterator tri_node_iterator;
	typedef typename TriGraph::edge_iterator tri_edge_iterator;
	typedef typename TriGraph::incident_iterator tri_incident_iterator;
	typedef typename std::set<Triangle>::iterator adj_triangle_iterator;

	// Define synonyms for the iterators and iterator types
	typedef TransformIterator<tri_node_iterator, Triangle> triangle_iterator;
	typedef TransformIterator<tri_edge_iterator, Link> link_iterator;
	typedef TransformIterator<vert_node_iterator, Vertex> vertex_iterator;
	typedef TransformIterator<vert_edge_iterator, Edge> edge_iterator;
	typedef TransformIterator<adj_triangle_iterator,Triangle> adjacent_iterator;
	typedef TransformIterator<tri_incident_iterator,Link> adj_link_iterator;


	// Externally visible node and triangle types, used in viewer maps
	typedef vert_node node_type;
	typedef Triangle tri_type;
	typedef Vertex vert_type;
	typedef Edge edge_type;

	// Helpful type definitions
	typedef std::set<Triangle> triangle_set;
	typedef std::vector<Edge> edge_list;

	/** triangle stores 3 vertices and user triangle data*/
	typedef struct triangle_data {
		vert_node n1_;
		vert_node n2_;
		vert_node n3_;
		UserTriData data_;
		triangle_data(vert_node n1, vert_node n2, vert_node n3, UserTriData d = UserTriData()) :
			n1_(n1), n2_(n2), n3_(n3), data_(d) {}
		triangle_data() {}
	} triangle_data;

	/** link stores its dual edge, and user edge data
	 * link must have dual
	 * */
	typedef struct link_data {
		vert_edge dual_;
		UserEdgeData data_;
		link_data(vert_edge dual, UserEdgeData data = UserEdgeData()) : dual_(dual), data_(data) {}
		link_data() {}
	} link_data;

	/** vertex stores its triangles, node data */
	typedef struct vertex_data {
		std::set<Triangle> triangles_;
		UserNodeData data_;	
		vertex_data(UserNodeData data = UserNodeData()) : data_(data) {/* triangles starts empty*/}
	} vertex_data;

	/** edge just stores its dual link
	 * edges don't necessarily have dual
	 * */
	typedef struct edge_data {
		tri_edge dual_;
		UserEdgeData data_;
		edge_data(tri_edge dual) : dual_(dual) {}
		edge_data() {};
	} edge_data;


	/** Type of indexes and sizes. Return type of Mesh::num_nodes(). */
	typedef unsigned size_type;

	/** Return the number of nodes in the mesh. */
	size_type num_nodes() const {
		return vertex_graph_.num_nodes();
	}

	/** Return the number of edges in the mesh. */
	size_type num_edges() const {
		return vertex_graph_.num_edges();
	}

	/** Return the number of triangles in the mesh. */
	size_type num_triangles() const {
		return triangle_graph_.num_nodes();
	}

	/** Return number of triangle-triangle connections in the mesh. */
	size_type num_links() const {
		return triangle_graph_.num_edges();
	}

	/** Vertex Type
	 * Used as a thin wrapper around Node
	 */
	class Vertex : private totally_ordered<Vertex> {
		public:
			Vertex() {};
			const Point& position() const {return v_.position();}
			idx_type index() const { v_.index();}
			bool operator==(const Vertex& other) const {return v_ == other.v_;}
			bool operator<(const Vertex& other) const { return v_ < other.v_;}
			UserNodeData& value() {return v_.value().data_;}
			size_type degree() const {return v_.degree();}
			incident_iterator edge_begin() {return v_.edge_begin();}
			incident_iterator edge_end() {return v_.edge_end();}

			/** Return an iterator to the first element in the adjacent 
			 * 	triangles set
			 */
			adjacent_iterator triangles_begin() {
				return adjacent_iterator(mesh_, v_.value().triangles_.begin());
			}

			/** Return an iterator to the last element in the adjacent 
			 * 	triangles set
			 */
			adjacent_iterator triangles_end() {
				return adjacent_iterator(mesh_, v_.value().triangles_.end());
			}

		private:
			friend class Mesh;
			const Mesh* mesh_;
			vert_node v_;
			Vertex(const Mesh* mesh, const vert_node v) : 
					mesh_(mesh), v_(v) {}
	};

	/** Link type
	 * Used as a thin wrapper around Graph edges
	 */
	class Link : private totally_ordered<Link> {
		public:
			Link() {};
			Triangle node1() const {return Triangle(mesh_, e_.node1());}
			Triangle node2() const {return Triangle(mesh_, e_.node2());}
			Triangle triangle1() const {return node1();};
			Triangle triangle2() const {return node2();};
			UserEdgeData& value() {
				return e_.value().data_;
			}
			bool operator==(const Link& other) const {
				return e_ == other.e_;
			}
			bool operator<(const Link& other) const {
				return e_ < other.e_;
			}
		private: 
			friend class Mesh;
			const Mesh* mesh_;
			tri_edge e_;
			Link(const Mesh* mesh, tri_edge e) : mesh_(mesh), e_(e) {}
	};

	/** Edge type
	 * Used as a thin wrapper around Graph edges
	 */
	class Edge : private totally_ordered<Edge> {
		public:
			Edge() {};
			Vertex node1() const {return Vertex(mesh_, e_.node1());}
			Vertex node2() const {return Vertex(mesh_, e_.node2());}
			Vertex vertex1() const {return node1();}
			Vertex vertex2() const {return node2();}
			UserEdgeData& value() { return e_.value().data_; }
			bool operator==(const Edge& other) const { return e_ == other.e_; }
			bool operator<(const Edge& other) const { return e_ < other.e_; }

			triangle_set adjacent_triangles() const {
				return mesh_->common_triangles(
							e_.node1(), e_.node2());
			}
			/** returns set of triangle in this edge that is not
			 * equal to @a t
			 * */
			triangle_set other_triangle(Triangle t) const {
				triangle_set group = adjacent_triangles();
				triangle_set other;
				for (auto it = group.begin(); it != group.end(); ++it)
					if (*it != t)
						other.insert(*it);
				return other;
			}

			double length() const {
				return e_.length();
			}

		private: 
			friend class Mesh;
			const Mesh* mesh_;
			vert_edge e_;
			Edge(const Mesh* mesh, vert_edge e) : mesh_(mesh), e_(e) {}
	};

	/** Triangle Type
	 * used as a proxy pattern to access nodes in our triangle_graph_
	 *
	 */
	class Triangle : private totally_ordered<Triangle> {
	public:
		/** construct invalid node
		 * use Mesh.add_triangle or Mesh.get_triangle instead
		 */
		Triangle() {
			// this is a hack for comparisons. replace once we have iterators
			uid_ = 10000;
			mesh_ = nullptr;
		}

		/** returns one of the 3 vertices of the triangle */
		vert_node vertex(size_type i) const {
			switch(i) {
			case 1:
				return mesh_->triangle_graph_.node(uid_).value().n1_;
			case 2:
				return mesh_->triangle_graph_.node(uid_).value().n2_;
			case 3:
				return mesh_->triangle_graph_.node(uid_).value().n3_;
			default:
				// invalid vertex option
				assert(false);
				return vert_node();
			}
		}

		/** calculates area of triangle from nodes using .5 |ab x ac | */
		double area() const {
			return 0.5 * norm(cross(vertex(2).position() - vertex(1).position(),
					vertex(3).position() - vertex(1).position()));
		}

		/** returns the user's value for the triangle. NOT our value! */
		UserTriData& value() {
			return mesh_->triangle_graph_.node(uid_).value().data_;
		}
		const UserTriData& value() const {
			return mesh_->triangle_graph_.node(uid_).value().data_;
		}

		/** returns how many connections this triangle has */
		size_type degree() const {
			return mesh_->triangle_graph_.node(uid_).degree();
		}

		adj_link_iterator link_begin() const {
			return adj_link_iterator(mesh_, 
					mesh_->get_tri_node(*this).edge_begin());
		}

		adj_link_iterator link_end() const {
			return adj_link_iterator(mesh_, 
					mesh_->get_tri_node(*this).edge_end());
		}

		/** comparison operators forward to underlying graph */
		bool operator==(const Triangle& x) const {
			return uid_ == x.uid_;
		}
		bool operator<(const Triangle& x) const {
			return uid_ < x.uid_;
		}

		Point position() const {
			double x_mid = (vertex(1).position().x + vertex(2).position().x
						 	+ vertex(3).position().x) / 3.0;
			double y_mid = (vertex(1).position().y + vertex(2).position().y
						 	+ vertex(3).position().y) / 3.0;
			double z_mid = (vertex(1).position().z + vertex(2).position().z
						 	+ vertex(3).position().z) / 3.0;
			return Point(x_mid, y_mid, z_mid);
		}

	private:
		friend class Mesh;

		// private data
		const Mesh* mesh_;
		size_type uid_;

		// private constructor, for Mesh functions
		Triangle(const Mesh* mesh, const size_type uid) : mesh_(mesh), uid_(uid) {}
		Triangle(const Mesh* mesh, const tri_node tn) : mesh_(mesh), uid_(tn.index()) {}
		Triangle(const Mesh* mesh, const Triangle t) : mesh_(mesh), uid_(t.uid_)
			{}

	};

	/* Returns a vector of edges adjacent to the triangle 
	 */ 
	edge_list edges(Triangle t) const {
		edge_list ret;
		ret.push_back(Edge(this, vertex_graph_.edge(t.vertex(1), t.vertex(2))));
		ret.push_back(Edge(this, vertex_graph_.edge(t.vertex(1), t.vertex(3))));
		ret.push_back(Edge(this, vertex_graph_.edge(t.vertex(2), t.vertex(3))));
		return ret;
	}
	/** add node
	 * function for the user to add nodes into the mesh
	 * should later be added into triangles
	 *
	 * @ post new num_nodes() = old num_nodes() + 1
	 */
	vert_node add_node(Point p, UserNodeData d = UserNodeData()) {
		return vertex_graph_.add_node(p, d);

	}

	/** add triangle
	 * turn existing nodes into a triangle
	 * @pre all nodes come from same graph
	 * @pre n1 != n2 != n3
	 * @pre no edge defined by n1, n2, n3 is already part of 2 triangles
	 *
	 * @RI each node should list this triangle as one of its triangles
	 * @RI any two triangles which share two points have a link between them
	 * @RI any link contains a reference to the edge defined by those 2 points
	 * @RI all nodes in a triangle have edges between them
	 *
	 * @post new num_triangles() = old num_triangles() + 1 if successful
	 * @post new num_links() = old num_links() + (0-3)
	 * @post new num_edges() = old num_edges() + (0-3)
	 * @post new num_edges() + new num_links() = old num_edges() + old num_links() + 3
	 *
	 * the sum of edges created + links created always == 3
	 */
	Triangle add_triangle(vert_node n1, vert_node n2, vert_node n3, UserTriData d = UserTriData()) {


		// check for existing triangles
		Triangle existing = get_triangle(n1, n2, n3);
		if (existing != Triangle()) {
			existing.value() = d;
			return existing;
		}
		// creating in private data
		tri_node tn = triangle_graph_.add_node(Point(),triangle_data(n1,n2,n3,d));

		//adding edges: some might already exist
		vert_edge e12 = vertex_graph_.add_edge(n2, n1);
		vert_edge e23 = vertex_graph_.add_edge(n2, n3);
		vert_edge e31 = vertex_graph_.add_edge(n3, n1);

		// find neighbors, create links and maintain duals
		std::set<Triangle> adj12 = common_triangles(n1, n2);
		if (adj12.size() != 0) {
			tri_edge temp = triangle_graph_.add_edge(tn, get_tri_node(*adj12.begin()), link_data(e12));
			e12.value().dual_ = temp;
		}
		std::set<Triangle> adj23 = common_triangles(n2, n3);
		if (adj23.size() != 0) {
			tri_edge temp = triangle_graph_.add_edge(tn, get_tri_node(*adj23.begin()), link_data(e23));
			e23.value().dual_ = temp;
		}
		std::set<Triangle> adj31 = common_triangles(n3, n1);
		if (adj31.size() != 0) {
			tri_edge temp = triangle_graph_.add_edge(tn, get_tri_node(*adj31.begin()), link_data(e31));
			e31.value().dual_ = temp;
		}

		// creating Triangle object
		Triangle tri = Triangle(this, tn);

		// attaching this triangle to each node
		n1.value().triangles_.insert(tri);
		n2.value().triangles_.insert(tri);
		n3.value().triangles_.insert(tri);

		return tri;
	}

	/** Returns all of the triangles that the vertices have in common */
	std::set<Triangle> common_triangles(vert_node n1, vert_node n2) const {
		std::set<Triangle> overlap;
		std::set_intersection(
				n1.value().triangles_.begin(), n1.value().triangles_.end(),
				n2.value().triangles_.begin(), n2.value().triangles_.end(),
				std::inserter(overlap, overlap.begin()));
		return overlap;
	}

	/** Return triangle that all 3 vertices have in common
	 * @pre A triangle exists that has all 3 nodes
	 * */
	Triangle get_triangle(vert_node n1, vert_node n2, vert_node n3) const {
		std::set<Triangle> overlap1;
		std::set_intersection(
			n1.value().triangles_.begin(), n1.value().triangles_.end(),
			n2.value().triangles_.begin(), n2.value().triangles_.end(),
			std::inserter(overlap1, overlap1.begin()));

		// finding intersection of intersection
		std::set<Triangle> overlap2;
		std::set_intersection(
			overlap1.begin(), overlap1.end(),
			n3.value().triangles_.begin(), n3.value().triangles_.end(),
			std::inserter(overlap2, overlap2.begin()));

		// 1 or 0
		if (overlap2.begin() != overlap2.end())
			return *overlap2.begin();
		else
			return Triangle();
	}

	/** Return a triangle that is bounded by the two edges
	 * @pre A triangle exists that contains the two edges
	 */
	Triangle get_triangle(vert_edge e1, vert_edge e2) const {
		return get_triangle(e1.node1(), e1.node2(), e2.node());
	}

	/** checks whether 3 nodes form a triangle */
	bool has_triangle(vert_node n1, vert_node n2, vert_node n3) const {
		Triangle t = get_triangle(n1, n2, n3);
		return t != Triangle();
	}


	/** Returns a vector normal to the edge, and normal to the +z direction
	 * If there is only one triangle, it points away from that triangle
	 * If the edge is bounded by two triangles, the normal points from the 
	 * 	smaller to the larger, as defined by operator<
	 */
	Point normal(Edge e) const {
		vert_edge edge = e.e_;
		triangle_set adj_triangles = e.adjacent_triangles();

		auto first = adj_triangles.begin();
		Triangle t1;
		if(adj_triangles.size() == 1) {
			t1 = *first;
		} 
		else if( adj_triangles.size() == 2) {
			auto first = adj_triangles.begin();
			t1 = std::min(*first, *(++first));	
		}
		else
			assert( false ); // every edge is bordered by a triangle
		Point edge_vec = edge.node1().position() - edge.node2().position();
		Point normal_vect = cross(edge_vec, Point(0,0,1));
		Point towards_t1 = get_unused(t1, edge).position() - edge.node1().position();
		if(dot(normal_vect, towards_t1) > 0)
			normal_vect = -1 * normal_vect;
		return normal_vect;
	}

	// other accessor for normal function
	Point normal(vert_edge e) const {
		normal(Edge(this, e));
	}

	// turns two triangles into edge
	Point normal(Triangle t1, Triangle t2) const {
		vert_edge edge = get_link(t1, t2).value().dual_;
		return normal(edge);
	}

	/** Creates a new iterator from an existing iterator IT by templating on 
	  * 	a arbitrary value type built from the value type of IT 
	  */
	template <class IT, class VALUE_TYPE>
	class TransformIterator : private totally_ordered<TransformIterator<IT,VALUE_TYPE>> {
		public:
			// These type definitions help us use STL's iterator_traits.
			/** Element type. */
			typedef VALUE_TYPE value_type;
			/** Type of pointers to elements. */
			typedef VALUE_TYPE* pointer;
			/** Type of references to elements. */
			typedef VALUE_TYPE& reference;
			/** Iterator category. */
			typedef std::input_iterator_tag iterator_category;
			/** Difference between iterators */
			typedef std::ptrdiff_t difference_type;
			/** Define the type of this iterator */
			typedef TransformIterator<IT, value_type> this_type;

			/** Construct an invalid NodeIterator. */
			TransformIterator() {
			}

			value_type operator*() const {return value_type(mesh_, *it_);}
			this_type& operator++() {++it_; return *this;}
			bool operator==(const this_type& other_) const {
				return it_ == other_.it_;
			}
			bool operator<(const this_type& other_) const {
				return it_ < other_.it_;
			}
		private: 
			friend class Mesh;
			const Mesh* mesh_;
			IT it_;
			TransformIterator(const Mesh* mesh, IT it) :
				mesh_(mesh), it_(it) {}
	};


	/** Returns an iterator to the first vertex in the graph */
	vertex_iterator vertex_begin() const {
		return vertex_iterator(this, vertex_graph_.node_begin());
	}

	/** Returns an iterator to the last vertex in the graph */
	vertex_iterator vertex_end() const {
		return vertex_iterator(this, vertex_graph_.node_end());
	}

	/** Return an iterator to the first edge in the graph */
	edge_iterator edge_begin() const {
		return edge_iterator(this, vertex_graph_.edge_begin());
	}
	
	/** Returns an iterator to the last edge in the graph */
	edge_iterator edge_end() const {
		return edge_iterator(this, vertex_graph_.edge_end());
	}

	/** Returns an iterator to the first triangle in the graph */
	triangle_iterator triangles_begin() const {
		return triangle_iterator(this, triangle_graph_.node_begin());
	}

	/** Returns an iterator to the last triangle in the graph */
	triangle_iterator triangles_end() const {
		return triangle_iterator(this, triangle_graph_.node_end());
	}

	/** Returns an iterator to the first link in the triangle graph */
	link_iterator link_begin() const {
		return link_iterator(this, triangle_graph_.edge_begin());
	}

	/** Returns an iterator to the last link in the triangle graph */
	link_iterator link_end() const {
		return link_iterator(this, triangle_graph_.edge_end());
	}

	/** Gets a Link object from two Triangles */
	Link link(Triangle t1, Triangle t2) const {
		return Link(this, get_link(t1, t2));
	}


private:
	friend class Triangle;


	/** private utility functions */

	/** returns the edge in the triangle graph between t1 and t2
	 * @pre there is a link between t1 and t2
	 */
	tri_edge get_link(Triangle t1, Triangle t2) {
		assert (triangle_graph_.has_edge(get_tri_node(t1), get_tri_node(t2)));
		return triangle_graph_.edge(get_tri_node(t1), get_tri_node(t2));
	}
	// cast from triangle to its base node
	tri_node get_tri_node(Triangle t) const {
		return triangle_graph_.node(t.uid_);
	}
	// returns the first vertex of t that is not part of edge e
	vert_node get_unused(Triangle t, vert_edge e) const {
		int i = 1;
		while (t.vertex(i) == e.node1() || t.vertex(i) == e.node2())
			++i;
		return t.vertex(i);
	}

	/** Private data members */
	TriGraph triangle_graph_;
	VertGraph vertex_graph_;

};
