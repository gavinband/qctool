
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <Eigen/Core>
#include <boost/math/distributions/hypergeometric.hpp>
#include "metro/FishersExactTest.hpp"

namespace metro {
	FishersExactTest::FishersExactTest( Eigen::Matrix2d const& matrix ):
		m_matrix( matrix ),
		m_distribution( m_matrix.col( 0 ).sum(), m_matrix.row( 0 ).sum(), m_matrix.sum() )
	{}

	double FishersExactTest::get_OR() const {
		return ( m_matrix(0,0) * m_matrix(1,1) ) / ( m_matrix(0,1)*m_matrix(1,0) ) ;
	}

	std::pair< double, double > FishersExactTest::get_confidence_interval() const {
		assert(0) ; // not implemented yet
	}

	double FishersExactTest::get_pvalue( Alternative const alternative ) const {
		using boost::math::cdf ;
		using boost::math::complement ;

		double result = 0 ;
		switch( alternative ) {
			case eGreater:
				if( m_matrix(0,0) == 0 || m_matrix(1,1) == 0 ) {
					result = 1 ;
				} else {
					result = cdf( complement( m_distribution, m_matrix( 0, 0 ) - 1.0 )) ;
				}
				break ;
			case eLess:
				result = cdf( m_distribution, m_matrix( 0, 0 ) ) ;
				break ;
			case eTwoSided:
				// As in R's fisher.test(), we take as "at least as extreme" all tables
				// with equal or lower probability under hypergeometric distribution.
				
				double const p = pdf( m_distribution, m_matrix( 0, 0 ) ) ;
				result = p ;

				// computation over all tables can be quite slow.
				// speed it up by inspecting which tail of the distribution our value is in
				// and then using the appropriate cdf.
				int const min_value = std::max( 0.0, m_matrix(0,0) - m_matrix(1,1)) ;
				int const lower_mode = std::max(
					int( ( m_matrix.col( 0 ).sum() + 1 ) * ( m_matrix.row( 0 ).sum() + 1 ) / ( m_matrix.sum() + 2 ) ),
					min_value
				) ;
				if( m_matrix( 0, 0 ) == lower_mode ) {
					result = 1 ;
				}
				else if( m_matrix( 0, 0 ) > lower_mode ) {
					double const min_a = std::max( 0.0, m_matrix(0,0) - m_matrix(1,1)) ;
					double const max_a = lower_mode ;
					result = cdf( complement( m_distribution, m_matrix( 0, 0 ) - 1.0 )) ;

					for(
						double a = min_a ;
						a <= max_a;
						++a
					) {
						double q = pdf( m_distribution, a ) ;
						if( q > p ) {
							// we are only summing lower tail, so shortcut if possible
							break ;
						}
						result += q ;
					}
				}
				else {
					double const min_a = lower_mode ;
					double const max_a = std::min( m_matrix.row(0).sum(), m_matrix.col(0).sum() ) ;
					result = cdf( m_distribution, m_matrix( 0, 0 ) ) ;

					for(
						double a = max_a ;
						a >= min_a;
						--a
					) {
						double q = pdf( m_distribution, a ) ;
						if( q > p ) {
							// we are only summing upper tail, so shortcut if possible
							break ;
						}
						result += q ;
					}
				}
				
			break ;
		}
		return result ;
	}
}
