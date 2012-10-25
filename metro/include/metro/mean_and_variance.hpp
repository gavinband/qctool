
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_MEAN_AND_VARIANCE
#define METRO_MEAN_AND_VARIANCE

#include <limits>

namespace metro {
	template< typename Data >
	std::pair< double, double > compute_mean_and_variance( Data const& data ) {
		double const mean = data.sum() / data.size() ;
		double variance = std::numeric_limits< double >::quiet_NaN() ;
		if( data.size() > 1 ) {
			variance = ( data.array() - mean ).square().sum() / ( data.size() - 1  ) ;
		}
		return std::make_pair( mean, variance ) ;
	}
}

#endif
