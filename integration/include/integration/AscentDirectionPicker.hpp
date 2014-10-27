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
#include "integration/ModifiedCholesky.hpp"

// #define DEBUG_ASCENT_DIRECTION_PICKER 1

namespace integration {
	template< typename Function >
	struct CholeskyOrEigenvalueSolver: public boost::noncopyable {
	public:
		typedef Eigen::VectorXd Vector ;
		typedef Eigen::MatrixXd Matrix ;
		
	public:
		CholeskyOrEigenvalueSolver(
			double const delta = -0.1
		):
			m_delta( delta )
		{
			assert( delta < 0 ) ;
		}

		void set_delta( double const delta ) {
			assert( delta < 0 ) ;
			m_delta = delta ;
		}
		
		Vector compute( Function& function, Vector const& point ) {
			function.evaluate_at( point, 2 ) ;
			Matrix const matrix = function.get_value_of_second_derivative() ;

			



			Vector const v = -function.get_value_of_first_derivative() ;
			m_cholesky_solver.compute( matrix ) ;
			if( m_cholesky_solver.info() == Eigen::Success && m_cholesky_solver.vectorD().array().maxCoeff() < 0 ) {
				return m_cholesky_solver.solve( v ) ;
			} else {
				m_eigen_solver.compute( matrix ) ;
				if( m_eigen_solver.info() == Eigen::NoConvergence ) {
					throw NumericalError( "integration::CholeskyOrEigenvalueSolver::solve()", "Eigenvalue decomposition did not converge" ) ;
				} 
			
				m_d = m_eigen_solver.eigenvalues() ;
				for( int i = 0; i < m_d.size(); ++i ) {
					if( m_d(i) > m_delta ) {
						m_d(i) = m_delta ;
					}
				}
				m_d = m_d.array().inverse() ;
				
				return (
					m_eigen_solver.eigenvectors() * m_d.asDiagonal() * m_eigen_solver.eigenvectors().transpose()
				) * v ;
			}
		}

	private:
		double m_delta ;
		Eigen::LDLT< Matrix > m_cholesky_solver ;
		Eigen::SelfAdjointEigenSolver< Matrix > m_eigen_solver ;
		Eigen::VectorXd m_d ;
	} ;
	
}

#endif
