
README.txt

<===== CollisionDetector ======>

Given a collection of meshes, our algorithm finds all "collisions"
between these meshes. A collision is defined as a point from one
mesh being inside another mesh, as checked by interior testing.

====== Example

Spheres turn red when they are in a collision

$ make demo_collision
% ./demo_collision

====== Files

CollisionDetector.hpp // main class
CollisionGeometry.hpp // 3d geometry
debug.hpp // nice wrappers for printing
MortonCoder.hpp // Thanks Vicente :) 
SpaceSearcher.hpp

demo_collision.hpp // demonstration 

====== Requirements

Collision Detector is templated on your mesh class.
Your mesh must supply the following functions / classes:

Mesh::Triangle
  Vertex({1,2,3})
Mesh::Node
Mesh::Edge

vertex_begin() // you probably named these nodes
vertex_end() 

====== Interface

// default constructor
CollisionDetector<Mesh> c;

// register meshes by reference to the collision checker
c.add_object(mesh1);
c.add_object(mesh2);

// checks for collisions between all meshes
// found collisions are stored in an internal vector
// that can be accessed by iterators. This vector is
// overwritten whenever you call check_collisions()
c.check_collisions();

// Each collision object represents a node that was inside a different mesh
auto boom = *c.begin();
Node n = boom.n1; // the node that was inside
Mesh* m1 = boom.mesh1; // the mesh that the node is a part of
Mesh* m2 = boom.mesh2; // the mesh that the node was inside

====== Fine tuned collision checking

If you want to check only certain objects against each other, you can
do so with "tags". For instance, you might want to check for collisions
between a tank and a bullet, but you don't want to check for collisions
between 2 bullets.

Tag bullet_tag = c.getOtherTag(); // checks against everything with a different tag
c.add_object(bullet1, bullet_tag);
c.add_object(bullet2, bullet_tag);

Through combined whitelists and blacklists, we can support arbitrary relationships.
The more conservative of two tags always takes precedence.

Tag tank_tag = c.getAllTag(); // checked against everything
Tag missile_tag = c.getNoneTag();
missile_tag.add(tank_tag); // missiles will only be checked against tanks
c.add_object(tank1, tank_tag);
c.add_object(missile1, missile_tag);