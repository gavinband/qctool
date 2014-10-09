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

// #define DEBUG_ASCENT_DIRECTION_PICKER 1

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
	template< typename FunctionEvaluator >
	struct CholeskyOrEigenvalueSolver: public boost::noncopyable {
	public:
		typedef Eigen::VectorXd Vector ;
		typedef Eigen::MatrixXd Matrix ;
		
	public:
		CholeskyOrEigenvalueSolver( double const delta = -0.01 ):
			m_delta( delta )
		{
			assert( delta < 0 ) ;
		}

		void set_delta( double const delta ) {
			assert( delta < 0 ) ;
			m_delta = delta ;
		}
		
		Vector compute( FunctionEvaluator& evaluator, Vector const& point ) {
			evaluator.evaluate_at( point, 2 ) ;
			Vector const& first_derivative = evaluator.get_value_of_first_derivative() ;
			Matrix const& second_derivative = evaluator.get_value_of_second_derivative() ;

#if DEBUG_ASCENT_DIRECTION_PICKER
				std::cerr << "CholeskyOrEigenvalueSolver: Solving:\n" << second_derivative << "...\n" ;
#endif
			m_cholesky_solver.compute( second_derivative ) ;
			if( m_cholesky_solver.info() == Eigen::Success && m_cholesky_solver.vectorD().array().maxCoeff() < 0 ) {
#if DEBUG_ASCENT_DIRECTION_PICKER
				std::cerr << "CholeskyOrEigenvalueSolver: LDLT decomposition computed.  Matrix was negative-definite.\n" ;
				std::cerr << "CholeskyOrEigenvalueSolver: LDLT decomposition L:\n" << m_cholesky_solver.matrixLDLT() << ".\n" ;
				std::cerr << "CholeskyOrEigenvalueSolver: LDLT decomposition D:\n" << m_cholesky_solver.vectorD().transpose() << ".\n" ;
#endif
				// Return Newton direction
				return m_cholesky_solver.solve( -first_derivative ) ;
			} else {
#if DEBUG_ASCENT_DIRECTION_PICKER
				std::cerr << "CholeskyOrEigenvalueSolver: LDLT decomposition not computed.  Matrix was positive-definite or indefinite.\n" ;
#endif
				m_eigen_solver.compute( second_derivative ) ;
				if( m_eigen_solver.info() == Eigen::NoConvergence ) {
					throw NumericalError( "integration::CholeskyOrEigenvalueSolver::operator()", "Eigenvalue decomposition did not converge" ) ;
				} 
			
				m_d = m_eigen_solver.eigenvalues() ;
#if DEBUG_ASCENT_DIRECTION_PICKER
				std::cerr << "CholeskyOrEigenvalueSolver: Eigenvalue decomposition computed.  Eigenvalues are: "
					<< m_d.transpose() << ".\n" ;
#endif
				double maxAbsEigenvalue = m_d.array().abs().maxCoeff() ;
				for( int i = 0; i < m_d.size(); ++i ) {
					if( m_d(i) > m_delta ) {
						m_d(i) = m_delta ;
					}
				}
				m_d = m_d.array().inverse() ;
				
				return -(
					m_eigen_solver.eigenvectors() * m_d.asDiagonal() * m_eigen_solver.eigenvectors().transpose()
				) * first_derivative ;
			}
		}

	private:
		Eigen::LDLT< Matrix > m_cholesky_solver ;
		Eigen::SelfAdjointEigenSolver< Matrix > m_eigen_solver ;
		double m_delta ;
		Eigen::VectorXd m_d ;
	} ;
}

#endif

