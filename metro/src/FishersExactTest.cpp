
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

		double result ;
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
				double const p = pdf( m_distribution, m_matrix( 0, 0 ) ) ;
				double const min_margin = std::min( m_matrix.row(0).sum(), m_matrix.col(0).sum() ) ;
				result = 0.0 ;
				for(
					double a = std::max( 0.0, m_matrix(0,0) - m_matrix(1,1));
					a <= min_margin;
					++a
				) {
					double q = pdf( m_distribution, a ) ;
					if( q <= p ) {
						result += q ;
					}
				}
			break ;
		}
		return result ;
	}
}
