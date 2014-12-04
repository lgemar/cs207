#include "CollisionDetector.hpp"

int main () {
	typedef Mesh<int, int, int> MeshType;
	typedef CollisionDetector<MeshType> collider;
	collider c = collider();
	c.test();
	return 0;
}
