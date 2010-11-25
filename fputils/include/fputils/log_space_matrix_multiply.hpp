#ifndef FPUTILS_LOG_SPACE_MATRIX_MULTIPLY_HPP
#define FPUTILS_LOG_SPACE_MATRIX_MULTIPLY_HPP

#include <boost/numeric/ublas/matrix.hpp>

namespace fputils {
	typedef boost::numeric::ublas::matrix< double > Matrix ;
	typedef boost::numeric::ublas::scalar_matrix<double> ConstantMatrix ;

	// function: log_multiply_exp
	// Work out the product of two matrices, but interpret all quantities as being in natural log space.
	// (Thus the result at elt. i,j is the log of the i,j th entry of the product of the matrices
	// obtained by exponentiating each elt of the left and right matrices).
	Matrix log_space_matrix_multiply( Matrix const& left, Matrix const& right ) ;
}

#endif
