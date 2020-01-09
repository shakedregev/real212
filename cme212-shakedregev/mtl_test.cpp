/**
 * @file mtl_test.cpp
 * Test script for interfacing with MTL4 and it's linear solvers.
 */

// HW3: Need to install/include Boost and MTL in Makefile
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

// HW3: YOUR CODE HERE
// Define a IdentityMatrix that interfaces with MTL
struct IdentityMatrix {
  /* Identity Matrix constructor */
	IdentityMatrix(size_t dim)
	: dim_(dim){} 

  /* * Helper function to perform multiplication. Allows for delayed
   * evaluation of results.
   * Assign::apply(a, b) resolves to an assignment operation such as
   a += b, a -= b, or a = b.
   * @pre @a size(v) == size(w) */
  template <typename VectorIn, typename VectorOut, typename Assign>
	void mult(const VectorIn& v , VectorOut& w, Assign) const{
		for(unsigned int i = 0; i != size(v); ++i)
		{ 
			w[i] = Assign::apply(w[i], v[i]); 
		}    
	}

  /* * Matvec forwards to MTL ’s lazy mat_cvec_multiplier operator */
  template<typename Vector>
	mtl::vec::mat_cvec_multiplier <IdentityMatrix, Vector>
	operator*(const Vector& v) const {
		return {*this, v};
	}

	// dimension getter for matrix
	size_t dimen() const
	{
		return dim_;
	}
private:
	size_t dim_;

};
/** The  number  of  elements  in the  matrix. */
inline std:: size_t  size(const IdentityMatrix& A)
{
	return A.dimen()*A.dimen();
}

/** The  number  of rows in the  matrix. */
inline std:: size_t  num_rows(const IdentityMatrix& A)
{
	return A.dimen();
}
/** The  number  of columns in the  matrix. */
inline std:: size_t  num_cols(const IdentityMatrix& A)
{
	return A.dimen();
}

/* * Traits that MTL uses to determine properties of our IdentityMatrix . */
namespace mtl{
	namespace ashape{

    /* * Define IdentityMatrix to be a non-scalar type . */
    template<>
		struct ashape_aux <IdentityMatrix> {
			typedef nonscal type;
		};
  } // end namespace ashape

  /* * IdentityMatrix implements the Collection concept
   * with value_type and size_type */
  template<>
  struct Collection <IdentityMatrix>{
  	typedef double value_type;
  	typedef unsigned size_type;
  };
} // end namespace mtl

int main()
{

  // HW3: YOUR CODE HERE
  // Construct an IdentityMatrix and "solve" Ix = b
  // using MTL's conjugate gradient solver
	const unsigned int N=1000;
	const double tol=1.e-8;
	// void(os);
	IdentityMatrix eye(N); //going for that sweet matlab pun syntax
	mtl::dense_vector<double> x(N, 1.0), b(N);
	b=x;
	x=0;
	// Terminates when |Ix-b|< tol * b or N iterations
	itl::cyclic_iteration<double> stop(b, N, tol);

  // Solve Ix = b
	itl::cg(eye, x, b, stop);
	return 0;
}
