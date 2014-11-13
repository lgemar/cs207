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

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"

#include "Point.hpp"
#include "Graph.hpp"

template <typename UserNodeData, typename UserEdgeData, typename UserTriData>
class Mesh {
	// HW3: YOUR CODE HERE
	// Write all typedefs, public member functions, private member functions,
	//   inner classes, private member data, etc.

	/** Data members that will be inside the 2 graphs */
	struct triangle_data;
	struct link_data;
	struct vertex_data;
	struct edge_data;
	class Triangle;

	// inner graphs themselves
	typedef Graph<triangle_data, link_data> TriGraph;
	typedef Graph<vertex_data, edge_data> VertGraph;

	// inner nodes and edges
	typedef TriGraph::node_type tri_node;
	typedef TriGraph::edge_type tri_edge;
	typedef VertGraph::node_type vert_node;
	typedef VertGraph::edge_type vert_edge;

	/** triangle stores 3 vertices and user triangle data*/
	typedef struct triangle_data {
		vert_node n1_;
		vert_node n2_;
		vert_node n3_;
		UserTriData data_;
		triangle_data(vert_node n1, vert_node n2, vert_node n3, UserTriData d) :
			n1_(n1), n2_(n2), n3_(n3), data_(d) {}
	} triangle_data;

	/** link stores its dual edge, and user edge data*/
	typedef struct link_data {
		vert_edge dual_;
		UserEdgeData data_;
	} link_data;

	/** vertex stores its triangles, node data */
	typedef struct vertex_data {
		std::vector<Triangle> triangles_;
		UserNodeData data_;
	} vertex_data;

	/** edge just stores its dual link */
	typedef struct edge_data {
		tri_edge dual_;
	} edge_data;

public:
	/** Type of indexes and sizes. Return type of Mesh::num_nodes(). */
	typedef unsigned size_type;

	/** Return the number of nodes in the mesh. */
	size_type num_nodes() const {
		return 0;
	}

	/** Return the number of edges in the mesh. */
	size_type num_edges() const {
		return 0;
	}

	/** Return the number of triangles in the mesh. */
	size_type num_triangles() const {
		return 0;
	}

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
		}

		/** returns one of the 3 vertices of the triangle */
		vert_node vertex(size_t i) {
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

		/** returns the user's value for the triangle. NOT our value! */
		UserTriData value() {
			return mesh_->triangle_graph_.node(uid_).value().data_;
		}

	private:
		friend class Mesh;

		// private data
		const Mesh* mesh_;
		size_t uid_;

		// private constructor, for Mesh functions
		Triangle(const Mesh* mesh, const size_t uid) : mesh_(mesh), uid_(uid) {}
		Triangle(const Mesh* mesh, const tri_node tn) : mesh_(mesh), uid_(tn.uid_) {}

	};

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

		// creating in private data
		tri_node tn = triangle_graph_.add_node(Point(),tri_data(n1,n2,n3,d));

		// creating Triangle object
		Triangle tri = Triangle(this, tn);

		// attaching this triangle to each node
		n1.value().triangles_.push_back(tri);
		n2.value().triangles_.push_back(tri);
		n3.value().triangles_.push_back(tri);

		// creating edges

		// creating links

	}

private:
	friend class Triangle;
	/** Private data members */
	TriGraph triangle_graph_;
	VertGraph vertex_graph_;

};
