//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef INTEGRATION_MODIFIED_CHOLESKY_HPP
#define INTEGRATION_MODIFIED_CHOLESKY_HPP

#include <iostream>
#include <boost/noncopyable.hpp>
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>

// #define DEBUG_MODIFIED_CHOLESKY 1

namespace integration {
	//
	// Implements the 'MC' algorithm from Gill & Murray, Practical Optimisation.
	// Also expressed as Algorithm 6.5, "Modified Cholesky algorithm" in
	// Nocedal & Wright, "Numerical Optimisation".
	//
	// I have borrowed a little bit from Eigen's LDLT.h here to get the permutations and types right.
	//
	// This algorithm only touches the lower-diagonal of the matrix.
	// As described in the above references, we use the lower-diagonal of the matrix
	// to store the non-unity entries of L, and the diagonal to store the entries of D.
	// The auxiliary variables c_ij are stored in the lower diagonal too until they
	// are overwritten by entries of L and D.
	//
	// At step j the matrix looks like:
	//
	//    0 . . j . .
	// 0  d
	// .  l d
	// .  l l d
	// j  c c c c
	// .  c c c a c
	// .  c c c a a c
	//
	// Where a refers to an entry of the original matrix (after the possible permutations) and l refers to an entry of
	// the computed matrix L, d refers to an entry of D, and c to one of the auxiliary c_ijs.
	//
	// At the jth step we update this to become
	// 1: d
	// .  l d
	// .  l l d
	// j: l l l d
	// .  c c c c c
	// .  c c c c a c
	//
	// i.e. we compute the jth row of L, the jth entry of D, and the c's in
	// the jth column and on the diagonal below j.
	template< typename Matrix >
	struct ModifiedCholesky {
	public:
	    enum {
	      RowsAtCompileTime = Matrix::RowsAtCompileTime,
	      ColsAtCompileTime = Matrix::ColsAtCompileTime,
	      Options = Matrix::Options & ~Eigen::RowMajorBit, // these are the options for the TmpMatrixType, we need a ColMajor matrix here!
	      MaxRowsAtCompileTime = Matrix::MaxRowsAtCompileTime,
	      MaxColsAtCompileTime = Matrix::MaxColsAtCompileTime,
	      UpLo = Eigen::Lower
	    } ;
	    typedef typename Matrix::Scalar Scalar;
	    typedef typename Eigen::NumTraits<typename Matrix::Scalar>::Real RealScalar;
	    typedef typename Matrix::Index Index;

	    typedef Eigen::Transpositions<RowsAtCompileTime, MaxRowsAtCompileTime> Transpositions;
	    typedef Eigen::PermutationMatrix<RowsAtCompileTime, MaxRowsAtCompileTime> Permutations;
	    typedef Eigen::TriangularView< Matrix const, Eigen::UnitLower > const MatrixL ;
		typedef Eigen::Diagonal< Matrix const > Diagonal ;

	public:
		ModifiedCholesky() {}

		ModifiedCholesky( ModifiedCholesky const& other ):
			m_matrix( other.m_matrix ),
			m_transpositions( other.m_transpositions )
		{}

		ModifiedCholesky& operator=( ModifiedCholesky const& other ) {
			m_matrix = other.m_matrix ;
			m_transpositions = other.m_transpositions ;
			return *this ;
		}
		
		ModifiedCholesky& compute( Matrix const& matrix ) {
			m_matrix = matrix ;
			compute_inplace( m_matrix ) ;
			return *this ;
		}
		
		MatrixL matrixL() const {
			return m_matrix ;
		}
		
	    Diagonal vectorD() const {
			return m_matrix.diagonal();
		}

		Transpositions const matrixP() const {
			return m_transpositions ;
		}
		
		Matrix solve( Matrix const& rhs ) const
	    {
			Matrix result = rhs ;
			assert( result.rows() == m_matrix.rows() ) ;
		    // result = P rhs
		    result = m_transpositions * result ;
		    // result = L^-1 (P rhs)
		    matrixL().solveInPlace( result );
		    // result = D^-1 (L^-1 P rhs)
			for( Index i = 0; i < m_matrix.rows(); ++i ) {
				result.row(i) /= m_matrix(i,i) ;
			}
		    // result = L^-T (D^-1 L^-1 P rhs)
		    matrixL().transpose().solveInPlace( result ) ;
		    // result = P^-1 L^-T (D^-1 L^-1 P rhs)
			result = m_transpositions.transpose() * result ;
			return result ;
	    }

	private:
	    Matrix m_matrix ;
	    Transpositions m_transpositions;
		
	private:
		
		void compute_inplace( Matrix& matrix ) {
			assert( matrix.rows() == matrix.cols() ) ;
			m_transpositions.resize( matrix.rows() ) ;
			Index const size = matrix.rows() ;
			
#if DEBUG_MODIFIED_CHOLESKY
			std::cerr << "ModifiedCholesky::compute_inplace(): computing with matrix:\n" << matrix << ".\n" ;
#endif
			
			Scalar biggestOnDiagonal ;
			RealScalar beta ;
			RealScalar betaSquared ;
			RealScalar delta ;

			for( Index j = 0; j < size; ++j ) {
		        // Find largest diagonal element
		        Index indexOfBiggestOnDiagonal ;
		        biggestOnDiagonal = matrix.diagonal().tail( size - j ).cwiseAbs().maxCoeff( &indexOfBiggestOnDiagonal ) ;
				indexOfBiggestOnDiagonal += j ;

				// Initialise beta and delta if we are starting.
		        if( j == 0 ) {
					Scalar biggestOffDiagonal = 0.0 ;
					for( Index i = 0; i < size; ++i ) {
						for( Index j = i+1; j < size; ++j ) {
							biggestOffDiagonal = std::max( biggestOffDiagonal, std::abs( matrix( i, j )) ) ;
						}
					}
					delta = std::numeric_limits< Scalar >::epsilon() * std::max( biggestOnDiagonal + biggestOffDiagonal, Scalar( 1 ) ) ;
					beta = std::max(
						std::numeric_limits< Scalar >::epsilon(),
						std::max( biggestOnDiagonal, biggestOffDiagonal / std::sqrt( ( size * size ) - 1 ) )
					) ;
					betaSquared = beta * beta ;
					
#if DEBUG_MODIFIED_CHOLESKY
					std::cerr << "ModifiedCholesky::compute_inplace(): initialised with delta = "
						<< delta << ", beta^2 = " << betaSquared << ".\n" ;
#endif
					
				}
				
#if DEBUG_MODIFIED_CHOLESKY
				std::cerr << "ModifiedCholesky::compute_inplace(): at iteration "
					<< j
						<< ": largest entry on diagonal = " << biggestOnDiagonal
					<< " at index " << indexOfBiggestOnDiagonal << ".\n" ;
#endif
											
				// Swap rows and columns corresponding to the jth and largest diagonal element.
		        m_transpositions.coeffRef( j ) = indexOfBiggestOnDiagonal ;
				Index const tailSize = size - indexOfBiggestOnDiagonal - 1 ;
		        if( j != indexOfBiggestOnDiagonal ) {
					// indexOfbiggestOnDiagonal is always >= j by construction
					// we only touch the lower triangular part of the matrix.
					matrix.row( j ).head( j ).swap( matrix.row( indexOfBiggestOnDiagonal ).head( j ) ) ;
					matrix.col( j ).tail( tailSize ).swap( matrix.col( indexOfBiggestOnDiagonal ).tail( tailSize ) ) ;
					std::swap( matrix.coeffRef(j,j), matrix.coeffRef( indexOfBiggestOnDiagonal, indexOfBiggestOnDiagonal ) );
					for( int i = j+1; i < indexOfBiggestOnDiagonal; ++i ) {
						std::swap( matrix.coeffRef( i, j ), matrix.coeffRef( indexOfBiggestOnDiagonal, i ) ) ;
					}
				}
				
				// We first compute the jth row of L to produce:
				// 1: d
				// .  l d
				// .  l l d
				// j: l l l c
				// .  c c c a c
				// .  c c c a a c
				//
				// by formula: l_js = c_js / d_s for s = 0,...,j-1.
		        matrix.row( j ).head( j ).array() /= matrix.diagonal().head( j ).array() ;
#if DEBUG_MODIFIED_CHOLESKY
				std::cerr << "ModifiedCholesky::compute_inplace(): after computing ls, matrix =\n"
					<< matrix << ".\n" ;
#endif
				// We next compute jth column of c_ijs to get:
				// 1: d
				// .  l d
				// .  l l d
				// j: l l l c
				// .  c c c c c
				// .  c c c c a c
				//
				// by formula c_ij = a_ij - sum_s l_ks c_is for s=j+1...n
				for( Index i = j+1; i < size; ++i ) {
					matrix(i,j) -= ( matrix.row( j ).head( j ) * matrix.row( i ).head( j ).transpose() ) ;
				}
				Scalar const theta = ( (j+1) == size ) ? 0.0 : ( matrix.col( j ).tail( tailSize ).maxCoeff() ) ;
#if DEBUG_MODIFIED_CHOLESKY
				std::cerr << "ModifiedCholesky::compute_inplace(): after computing cs, matrix =\n"
					<< matrix << ",\n"
					<< "theta = " << theta << ".\n" ;
#endif

				// compute d_jj to produce
				// 1: d
				// .  l d
				// .  l l d
				// j: l l l d
				// .  c c c c c
				// .  c c c c a c
				double new_dj = std::max( delta, std::abs( matrix(j,j) ) ) ;
				new_dj = std::max( new_dj, (theta*theta) / betaSquared ) ;
				matrix(j,j) = new_dj ;
				//
				// Finally update the c_ii's
				matrix.diagonal().tail( tailSize ).array() -= ( matrix.col(j).tail( tailSize ).array().square() ) / matrix(j,j) ;
				
#if DEBUG_MODIFIED_CHOLESKY
				std::cerr << "ModifiedCholesky::compute_inplace(): after iteration " << j << ", matrix is:\n"
					<< matrix << ".\n" ;
#endif
			}
		}
	} ;
}

#endif
