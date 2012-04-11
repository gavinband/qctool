
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_MULTIVARIATENORMAL_DISTRIBUTION_HPP
#define QCTOOL_MULTIVARIATENORMAL_DISTRIBUTION_HPP

#include <Eigen/Core>
#include <Eigen/Cholesky>
#include "Distribution.hpp"

template< typename Vector, typename Matrix >
struct MultivariateNormalDistribution: public Distribution< Vector, Matrix > {
	MultivariateNormalDistribution( Vector const& mean, Matrix const& variance ):
		m_mean( mean ),
		m_variance( variance ),
		m_eigenvalue_tolerance( 0.000000001 )
	{
		assert( m_mean.cols() == 1 ) ;
		assert( m_variance.rows() == m_variance.cols() == m_mean.rows() ) ;
	}


	typedef boost::function< void( std::string const& ) > NameSetter ;

	double get_log_density_at( Vector const& v ) const {
		// decomposition is
		//    variance = P^t L D L^t P
		// we need 
		Vector a = m_decomposition.matrixL().solve( m_decomposition.transpositionsP() * ( v - m_mean ) ) ;
		double result = a.transpose() * m_decomposition.vectorD().asDiagonal() * a ;
		double const PI = 3.14159265358979323846264338327950288 ;

		double rank = 0.0 ;
		for( int i = 0; i < m_decomposition.vectorD().size(); ++i ) {
			if( m_decomposition.vectorD()(i) > 0.0 ) {
				++rank ;
			}
		}
		double constant = -0.5 * ( rank * std::log( 2.0 * PI ) + std::log( m_decomposition.vectorD().prod() ) ) ;
		return constant + result ;
	}

private:
	Vector m_mean ;
	Matrix m_variance ;	
	Eigen::LDLT< Matrix > m_decomposition ;
	double const m_eigenvalue_tolerance ;
	double m_determinant ;
	double m_rank ;
} ;

#endif
