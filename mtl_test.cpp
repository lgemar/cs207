/**
 * @file mtl_test.cpp
 * Test script for interfacing with MTL4 and it's linear solvers.
 */

// HW3: Need to install/include Boost and MTL in Makefile
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

// Define a IdentityMatrix that interfaces with MTL
struct IdentityMatrix {
	IdentityMatrix(int m, int n) : m(m), n(n), s(m * n) {
	}

	/** Helper function to perform delayed evalutation of multiplication. 
	 * Assign::apply(a, b) resolves to an assignment such as a += b, a-=b, 
	 * 	or a = b
	 * @pre size(v) == size(w)
	 */
	template<typename VectorIn, typename VectorOut, typename Assign>
	void mult(const VectorIn& v, VectorOut& w, Assign) const {
		assert(size(v) == size(w));

		size_t highest_index = size(v);
		for (size_t i = 0; i < highest_index; i++) {
			Assign::apply(w[i], v[i]);
		}
	}

	/** Matrix-vector multiplication forwards to MTL's lazy mat_cvec_mult
	 * operator */
	template <typename Vector>
	mtl::vec::mat_cvec_multiplier<IdentityMatrix, Vector> operator*(const Vector& v) const {
		return mtl::vec::mat_cvec_multiplier
					<IdentityMatrix, Vector>(*this, v);
	}

	// number of rows / columns / size
	int m;
	int n;
	int s;
};

inline std::size_t size(const IdentityMatrix& A) { return A.s * A.s; }
inline std::size_t num_rows(const IdentityMatrix& A) { return A.m; }
inline std::size_t num_cols(const IdentityMatrix& A) { return A.n; }

/** Traits that mtl uses to determine properties of our Identity matrix */
namespace mtl {

/** Define IdentityMatrix to be a non-scalar type */
namespace ashape {
template<>
struct ashape_aux<IdentityMatrix> {
	typedef nonscal type;
};
} // end namespace ashape

/** IdentityMatrix implements Collection concept with value type 
 * and size type */
template<>
struct Collection<IdentityMatrix> {
	typedef double value_type;
	typedef unsigned size_type;
};
}


int main()
{
  // Construct an IdentityMatrix and "solve" Ix = b
  // using MTL's conjugate gradient solver

  // define sz
  const int sz = 40, N = sz * sz;

  //typedef mtl::mat::poisson2D_dirichlet matrix_type;
  typedef IdentityMatrix matrix_type;

  // Set up an identity matrix, A
  matrix_type I(N, N);

  // Create a preconditioner
  itl::pc::identity<matrix_type> P(I);

  // Set up a matrix
  mtl::dense_vector<double> x(N, 1.0), b(N);

  b = I * x;
  x = 0;

  itl::cyclic_iteration<double> iter(b, 50, 1e-10);

  cg(I, x, b, P, iter);

  // Print the results
  assert( x == b );
  std::cout << "Success" << std::endl;

  return 0;
}
