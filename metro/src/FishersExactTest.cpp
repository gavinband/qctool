
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
				else {
					// compute one tail using boost::cdf() and the other using
					// a recursion trick.
					double begin_a = 0 ;
					double end_a = 0 ;
					double direction = 0 ;
					if( m_matrix( 0, 0 ) > lower_mode ) {
						// compute upper tail
						result = cdf( complement( m_distribution, m_matrix( 0, 0 ) - 1.0 )) ;
						// compute lower tail, walking out from mode towards tail.
						begin_a = lower_mode ;
						end_a = std::max( 0.0, m_matrix(0,0) - m_matrix(1,1)) - 1 ;
						direction = -1 ;
					}
					else {
						// compute lower tail
						result = cdf( m_distribution, m_matrix( 0, 0 ) ) ;
						// compute upper tail, walking out from mode towards tail.
						begin_a = lower_mode ;
						end_a = std::min( m_matrix.row(0).sum(), m_matrix.col(0).sum() ) + 1 ;
						direction = 1 ;
					}
				
					double const R0 = m_matrix.row(0).sum() ;
					double const C0 = m_matrix.col(0).sum() ;
					double const C1 = m_matrix.col(1).sum() ;
					for(
						double a = begin_a, q = pdf( m_distribution, begin_a ) ;
						(direction * a ) < (direction * end_a) ;
						a += direction
					) {
						if( q == 0 ) {
							break ;
						}
						else if( ( q / p ) < 10 && ( q / p ) > 0.1 ) {
							// Within an order of magnitude of p.  Recompute q to avoid accuracy loss.
							q = pdf( m_distribution, a ) ;
						}
						if( q <= p ) {
							result += q ;
						}
						double const b = R0 - a ;
						double const c = C0 - a ;
						double const d = C1 - b ;
						if( direction == 1 ) {
							q *= ( b / ( a + 1 ) ) * ( c / ( d + 1 )) ;
						} else {
							q /= ( ( b + 1 ) / a ) * ( ( c + 1 ) / d ) ;
						}
					}
				}
			break ;
		}
		return result ;
	}
}
