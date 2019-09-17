
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_SUMMATION_HPP
#define METRO_SUMMATION_HPP

namespace metro {
	template< typename Vector >
	double naive_sum( Vector const& data ) {
		return data.sum() ;
	}

	template< typename Vector >
	double neumaier_sum( Vector const& data ) {
		double sum = 0.0 ;
		double c = 0.0 ;
		for( int i = 0; i < data.size(); ++i ) {
			double t = sum + data[i] ;
			if( std::abs( sum ) > std::abs( data[i] ) ) {
				c += (sum - t) + data[i] ;
			} else {
				c += (data[i] - t) + sum ;
			}
			sum = t ;
		}
		return( sum + c ) ;
	}
}

#endif
