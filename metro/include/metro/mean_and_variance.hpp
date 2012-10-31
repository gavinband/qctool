
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_MEAN_AND_VARIANCE
#define METRO_MEAN_AND_VARIANCE

#include <limits>

namespace metro {
	// Compute mean and variance for a vector.
	// All values must be non-missing (i.e. not NaN.)
	template< typename Data >
	std::pair< double, double > compute_mean_and_variance( Data const& data ) {
		double const mean = data.sum() / data.size() ;
		double variance = std::numeric_limits< double >::quiet_NaN() ;
		if( data.size() > 1 ) {
			variance = ( data.array() - mean ).square().sum() / ( data.size() - 1  ) ;
		}
		return std::make_pair( mean, variance ) ;
	}

	// Compute mean and variance for a vector ignoring missing values.
	template< typename Data >
	std::pair< double, double > compute_mean_and_variance( Data const& data, Data const& nonmissingness ) {
		assert( data.size() == nonmissingness.size() ) ;
		// Ensure all non-missing values are zero.
		double const mean = ( data.array() * nonmissingness.array() ).sum() / nonmissingness.sum() ;
		double variance = std::numeric_limits< double >::quiet_NaN() ;
		if( data.size() > 1 ) {
			variance = (( data - mean ).array() * nonmissingness.array() ).square().sum() / ( nonmissingness.sum() - 1 ) ;
		}
		return std::make_pair( mean, variance ) ;
	}
}

#endif
