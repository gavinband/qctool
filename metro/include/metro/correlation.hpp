
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_CORRELATION
#define METRO_CORRELATION

#include <limits>

namespace metro {
	template< typename Vector1, typename Vector2, typename NonMissingVector >
	double compute_correlation( Vector1 const& v1, Vector2 const& v2, NonMissingVector const& non_missingness_indicator ) {
		assert( v1.size() == v2.size() ) ;
		double non_missingness = non_missingness_indicator.sum() ;
		double mean1 = 0.0 ;
		double mean2 = 0.0 ;
		for( int i = 0; i < v1.size(); ++i ) {
			if( non_missingness_indicator( i )) {
				mean1 += v1(i) / non_missingness ;
				mean2 += v2(i) / non_missingness ;
			}
		}

		double covariance = 0.0 ;
		double variance1 = 0.0 ;
		double variance2 = 0.0 ;
		for( int i = 0; i < v1.size(); ++i ) {
			if( non_missingness_indicator( i )) {
				covariance += ( v1(i) - mean1 ) * ( v2(i) - mean2 ) ;
				variance1 += ( v1(i) - mean1 ) * ( v1(i) - mean1 ) ;
				variance2 += ( v2(i) - mean2 ) * ( v2(i) - mean2 ) ;
			}
		}
		
		// We should divide the covariance by N-1 and also
		// divide each variance by the same quantity.
		// But this washes out in the ratio.
		
		return covariance / std::sqrt( variance1 * variance2 ) ;
	}
}

#endif
