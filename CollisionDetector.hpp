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

	/** 3D object that can be checked for collisions
	 * @a mesh the closed mesh that defines the boundary of the object
	 * @a tag the tag representing how we should check it
	 *
	 * these are what we operate over
	 */
	struct Object {
		MeshType& mesh;
		Tag& tag;

		// need this for some reason
		const Object& operator=(const Object& o) {
			return o;
		}

		// default to all Tag
		Object(MeshType& m) : mesh(m), tag(Tag()) {};

		// explicitly pass in a tag
		Object(MeshType& m, Tag& t) : mesh(m), tag(t) {};

	};

	/** struct to store collision information
	 * @a n_ the node that is inside another mesh
	 * @a t_ the triangle in the mesh that is closest to the node
	 * 	(most likely guess for the collision)
	 */
	struct Collision {
		Node n_;
		Triangle t_;

		Collision(Node n, Triangle t) : n_(n), t_(t) {}
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
	private:
		friend struct CollisionDetector;
		int id_;
		bool white_;

		/** private constructor for collision detector */
		Tag(int id, bool white) : id_(id), white_(white) {}

	public:
		/** create invalid tag */
		Tag() : id_(0), white_(false) {}
		std::vector<size_t> list_;
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
	void add_object(MeshType m, Tag tag) {

		// create a node for this mesh
		Point approx_pos = (*(m.vertex_begin())).position();
		Object o = Object(m, tag);
		object_node n = object_graph_.add_node(approx_pos, o);

		// saving link to this mesh in hash table
		mesh2node[&m] = n;

		// create edges for the graph
		for(auto it = object_graph_.node_begin(); 
				 it != object_graph_.node_end(); ++it) {

			// check if other tag in list
			Tag tag2 = (*it).value().tag;
			auto optr = std::find(tag.list_.begin(), tag.list_.end(), tag2.id_);
			bool tag_found = (optr != tag.list_.end());

			// white listed and not found
			if (tag.white_ && !tag_found	)
				continue;

			// blacklisted and found
			if (!tag.white_ && tag_found)
				continue;


			// checking if it's in other tag's list
			optr = std::find(tag2.list_.begin(), tag2.list_.end(), tag.id_);
			tag_found = (optr != tag.list_.end());

			// white listed and not found
			if (tag2.white_ && !tag_found	)
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
	void remove_object(MeshType& m) {
		db("in remove");
		// finding node from mesh
		object_node on = mesh2node[&m];
		db("got node");
		object_graph_.remove_node(on);
		db("removed node");
		mesh2node.erase(&m);
		db("erased mesh");
	}

	/** Finds all collisions within the meshes defined by the range
	 * store them in our internal collisions array
	 * @pre @a first and @a last must define a valid iterator range
	 */
	template <typename IT>
	void check_collisions() {
		// Build bounding boxes out of all of the meshes
		for(auto it = object_graph_.node_begin(); it != object_graph_.node_end(); ++it) {
			Object obj = (*it).value();

			BoundingBox b = build_bb(obj.mesh.vertex_begin(),
									 obj.mesh.vertex_end());
			bounding_boxes_.push_back(std::make_pair(b,obj));
		}

		// do other stuff
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

	private: 
		std::vector<std::pair<BoundingBox, Object>> bounding_boxes_;
		std::vector<Collision> collisions_;
		size_t next_tag_id_;
		object_graph object_graph_;
		std::unordered_map<MeshType*, object_node> mesh2node;

};
