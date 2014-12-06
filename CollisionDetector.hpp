#include<iostream>
#include<string>
#include<cmath>
#include<cassert>
#include<vector>
#include<unordered_map>

#include "Point.hpp"
#include "BoundingBox.hpp"
#include "CollisionGeometry.hpp"
#include "Graph.hpp"
#include "Mesh.hpp"
#include "debug.hpp"
#include "SpaceSearcher.hpp"

/** Simple mapping from a node to its corresponding point. */
struct NodeToPoint {
  template <typename NODE>
  Point operator()(const NODE& n) {
    return n.position();
  }
};

template <typename MeshType>
struct CollisionDetector {

	// used types ----------------------
	class Tag;
	struct Collision;
	struct Object;
	typedef typename MeshType::Triangle Triangle;
	typedef typename MeshType::Node Node;
	typedef typename MeshType::Edge Edge;
	typedef typename std::vector<Collision>::iterator CollIter;
	typedef class Graph<Object,int> object_graph;
	typedef typename object_graph::Node object_node;
	typedef SpaceSearcher<Node, NodeToPoint> space_searcher;

	CollisionDetector() : next_tag_id_(1) {
	}

	/** 3D object that can be checked for collisions
	 * @a mesh the closed mesh that defines the boundary of the object
	 * @a tag the tag representing how we should check it
	 *
	 * these are what we operate over
	 */
	struct Object {
		MeshType* mesh;
		Tag tag;

		// explicitly pass in a tag
		Object(MeshType& m, Tag& t) : mesh(&m), tag(t) {};

	};

	/** struct to store collision information
	 * @a n_ the node that is inside another mesh
	 * @a t_ the triangle in the mesh that is closest to the node
	 * 	(most likely guess for the collision)
	 */
	struct Collision {

		Node n;
		Triangle t;

		Collision(Node n, Triangle t) : n(n), t(t) {}
		Collision(Node n) : n(n), t() {}
	};

	/** Tag type to specify what gets checked for what
	 * @a id_ the unique id to represent this tag
	 * @a white_ whether this is a white list or a black list
	 * @a list_ a list of other Tags to specify relationship with this tag
	 *
	 * @RI no two tags have the same id
	 *
	 * Default tag is (0, false)
	 *
	 * White list tags will only check collisions against tags on their list
	 * Black list tags will check collisions against any tag not on their list
	 *
	 * Decision are made by the more conservative tag. i.e. a tag that checks everything
	 * will not be checked against a tag that checks nothing.
	 * All checking must be bidirectional.
	 */
	class Tag {
	//private:
	public:
		//friend struct CollisionDetector;
		int id_;
		bool white_;
		std::vector<size_t> list_;

		/** private constructor for collision detector */
		Tag(int id, bool white) : id_(id), white_(white) {}

	public:
		/** create invalid tag */
		Tag() : id_(0), white_(false) {}
		void add(Tag t) {
			list_.push_back(t.id_);
		}
	};

	/** tag has an empty blacklist, will check against everything */
	Tag getAllTag() {
		return Tag();
	}

	/** tag has empty white list, will check against nothing */
	Tag getNoneTag() {
		return get_tag(true);
	}

	/** tag will check only against self */
	Tag getSelfTag() {
		Tag t = get_tag(true);
		t.list_.push_back(t.id_);
		return t;
	}

	/** tag will check anything but self */
	Tag getOtherTag() {
		Tag t = get_tag(false);
		t.list_.push_back(t.id_);
		return t;
	}

	Tag get_tag(bool white) {
		size_t id = next_tag_id_;
		++next_tag_id_;
		return Tag(id, white);

	}

	/** Adds an object to the world of objects */
	void add_object(MeshType& m, Tag tag = Tag()) {

		// create a node for this mesh
		Point approx_pos = (*(m.vertex_begin())).position();
		Object o = Object(m, tag);
		object_node n = object_graph_.add_node(approx_pos, o);

		// saving link to this mesh in hash table
		mesh2node[&m] = n;

		// create edges for the graph
		for(auto it = object_graph_.node_begin(); 
				 it != object_graph_.node_end(); ++it) {

			Tag tag2 = (*it).value().tag;

			// check if other tag in list
			auto optr = std::find(tag.list_.begin(), tag.list_.end(), tag2.id_);
			bool tag_found = (optr != tag.list_.end());

			// white listed and not found
			if (tag.white_ && !tag_found)
				continue;

			// blacklisted and found
			if (!tag.white_ && tag_found)
				continue;

			// checking if it's in other tag's list
			optr = std::find(tag2.list_.begin(), tag2.list_.end(), tag.id_);
			tag_found = (optr != tag2.list_.end());

			// white listed and not found
			if (tag2.white_ && !tag_found)
				continue;

			// blacklisted and found
			if (!tag2.white_ && tag_found)
				continue;

			// self loop
			if (*it == n)
				continue;

			// no conflicts found, adding edge
			object_graph_.add_edge(n, (*it));
		}

	}

	/** removes a mesh from our collision detection
	 *
	 */
	void remove_object(const MeshType& m) {
		// finding node from mesh
		object_node on = mesh2node[&m];
		object_graph_.remove_node(on);
		mesh2node.erase(&m);
	}

	/** Finds all collisions within the meshes defined by the range
	 * store them in our internal collisions array
	 * @pre @a first and @a last must define a valid iterator range
	 */
	void check_collisions() {
		collisions_.clear();
		// Iterate over all edges in the graph and find intersections
		for(auto it = object_graph_.edge_begin(); 
				it != object_graph_.edge_end(); ++it) {
			// Get the two meshes that are a part of this edge
			auto e = (*it);
			auto m1 = *e.node1().value().mesh;
			auto m2 = *e.node2().value().mesh;
			// Build spatial search objects
			space_searcher s1 = space_searcher(m1.vertex_begin(), 
											m1.vertex_end(),
											NodeToPoint());
			space_searcher s2 = space_searcher(m2.vertex_begin(), 
											m2.vertex_end(),
											NodeToPoint());
			// Find the bounding boxes corresponding to spaces
			BoundingBox bb1 = s1.bounding_box();
			BoundingBox bb2 = s2.bounding_box();
			// Find the collisions
			find_collisions(s1.begin(bb2), s1.end(bb2), m2);
			find_collisions(s2.begin(bb1), s2.end(bb1), m1);
		}
	}

	/** Add the collision information to the collisions array */
	template<typename IT, typename MESH>
	int find_collisions(IT first, IT last, MESH m) {
		int num_collisions = 0;
		for(IT it1 = first; it1 != last; ++it1) {
			int num_intersections = 0;
			Node n = (*it1);
			
			// Find the two points that define a line
			Point p0 = n.position();
			Point p1 = 2 * p0;
			for(auto it2 = m.triangles_begin(); 
					it2 != m.triangles_end(); ++it2) {
					Triangle t = (*it2);
					// Find the three points that make up triangle
					db("got here");
					Point t1 = t.vertex(1).position();
					Point t2 = t.vertex(2).position();
					Point t3 = t.vertex(3).position();

					// Determine intersection point
					Point p;
					if( is_plane_line_intersect(t1, t2, t3, p0, p1)) {
						p = plane_line_intersect(t1, t2, t3, 
														p0, p1);

						// Check if intersection points same direction
						// as the outgoing ray
						bool check1 = dot(p-p0, p1-p0) > 0;

						// Check if intersection point is inside of
						// the triangle being checked against
						bool check2 = is_inside_triangle(t1, t2, t3, p);

						// Increase the num_intersections is the two
						// check are true
						if( check1 && check2 ) {
							++num_intersections;
						}
					}
			}

			// If the number of intersections is odd, add to collisions
			if( num_intersections % 2 != 0 ) {
				Collision c = Collision(n);
				collisions_.push_back(c);
				++num_collisions;
			}
		}
		return num_collisions;
	}

	/** returns iterator to beginning of our found collisions
	 */
	CollIter begin() {
		return collisions_.begin();
	}

	/** returns iterator to end of our vector of collisions
	 */
	CollIter end() {
		return collisions_.end();
	}

	void print_graph() const {
		object_graph_.print_graph();
	}

	private: 
		std::vector<std::pair<BoundingBox, Object>> bounding_boxes_;
		std::vector<Collision> collisions_;
		size_t next_tag_id_;
		object_graph object_graph_;
		std::unordered_map<const MeshType*, object_node> mesh2node;

};
