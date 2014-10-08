//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef CORE_ASCENT_DIRECTION_PICKER_HPP
#define CORE_ASCENT_DIRECTION_PICKER_HPP

#include <iostream>
#include <boost/noncopyable.hpp>
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>

#define DEBUG_ASCENT_DIRECTION_PICKER 1

namespace integration {
	//
	// This class implements a solver which decomposes a symmetric matrix
	// using Eigen's Cholesky decomposition if that works
	// (i.e. if matrix is positive-definite )
	// or falls back to a symmetric eigenvalue decomposition if not.
	// In the latter case a positive-definite approximation to the
	// inverse is provded by setting all eigenvalues to be at least some minimal
	// positive value delta.
	// We expect the matrix passed in to usually be negative-definite; it will be
	// a second derivative of a likelihood function.
	struct CholeskyOrEigenvalueSolver: public boost::noncopyable {
	public:
		typedef Eigen::VectorXd Vector ;
		typedef Eigen::MatrixXd Matrix ;
		
		enum State { eUncomputed, eCholeskyComputed, eEigenComputed } ;

	public:
		CholeskyOrEigenvalueSolver( double const delta = 0.01 ):
			m_delta( delta )
		{}

		void compute( Matrix const& matrix ) {
#if DEBUG_ASCENT_DIRECTION_PICKER
				std::cerr << "CholeskyOrEigenvalueSolver: Solving:\n" << matrix << "...\n" ;
#endif
			m_cholesky_solver.compute( matrix ) ;
			if( m_cholesky_solver.info() == Eigen::Success && m_cholesky_solver.isNegative() ) {
				m_state = eCholeskyComputed ;
#if DEBUG_ASCENT_DIRECTION_PICKER
				std::cerr << "CholeskyOrEigenvalueSolver: LDLT decomposition computed.  Matrix was negative-definite." ;
#endif
			} else {
#if DEBUG_ASCENT_DIRECTION_PICKER
				std::cerr << "CholeskyOrEigenvalueSolver: LDLT decomposition not computed.  Matrix was positive-definite or indefinite.\n" ;
#endif
				m_eigen_solver.compute( matrix ) ;
				if( m_eigen_solver.info() == Eigen::NoConvergence ) {
					m_state = eUncomputed ;
				} else {
					m_state = eEigenComputed ;
				
					m_inverse_eigenvalues = m_eigen_solver.eigenvalues() ;
#if DEBUG_ASCENT_DIRECTION_PICKER
					std::cerr << "CholeskyOrEigenvalueSolver: Eigenvalue decomposition computed.  Eigenvalues are: "
						<< m_inverse_eigenvalues.transpose() << ".\n" ;
#endif
					double maxAbsEigenvalue = m_inverse_eigenvalues.array().abs().maxCoeff() ;
					for( int i = 0; i < m_inverse_eigenvalues.size(); ++i ) {
						if( m_inverse_eigenvalues(i) > -m_delta ) {
							m_inverse_eigenvalues(i) = -m_delta ;
						}
					}
					m_inverse_eigenvalues = m_inverse_eigenvalues.array().inverse() ;
				}
			}
		}

		Vector solve( Vector const& v ) {
			assert( m_state != eUncomputed ) ;
			switch( m_state ) {
				case eUncomputed:
					assert( 0 ) ;
					break ;
				case eCholeskyComputed:
					return m_cholesky_solver.solve( v ) ;
					break ;
				case eEigenComputed:
					return (
						m_eigen_solver.eigenvectors() * m_inverse_eigenvalues.asDiagonal() * m_eigen_solver.eigenvectors().transpose()
					) * v ;
					break ;
			}
		}

	private:
		Eigen::LDLT< Matrix > m_cholesky_solver ;
		Eigen::SelfAdjointEigenSolver< Matrix > m_eigen_solver ;
		double const m_delta ;

		// state variables
		State m_state ;
		Eigen::VectorXd m_inverse_eigenvalues ;
	} ;
}

#endif

