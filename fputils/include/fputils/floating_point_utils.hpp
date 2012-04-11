
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef FLOATING_POINT_UTILS_HPP
#define FLOATING_POINT_UTILS_HPP

#include <utility>
#include <iterator>
#include <cmath>
#include <limits>
#include <vector>

namespace fputils {
	template< typename F >
	bool floats_are_equal( F a, F b, F tolerance ) {
		//See http://www.boost.org/doc/libs/1_34_1/libs/test/doc/components/test_tools/floating_point_comparison.html
	    F difference = std::abs( a - b ) ;
	    return (difference <= tolerance * std::abs( a ))
	        && (difference <= tolerance * std::abs( b )) ;
	}

	template< typename F >
	bool floats_are_equal_to_within_epsilon( F a, F b, F epsilon ) {
		if( std::abs( a ) == std::numeric_limits< F >::infinity() || std::abs( b ) == std::numeric_limits< F >::infinity() )  {
			return a == b ;
		}
		else {
	    	return (std::abs( a - b ) < epsilon ) ;
		}
	}

	template< typename F >
	bool check_if_less_than_with_tolerance( F a, F b, F tolerance ) {
		return a < (b + tolerance) ;
	}

	template< typename F >
	bool check_if_even( F a, F tolerance = 0.0 ) {
		return floats_are_equal( std::fmod( a, 2.0 ), 0.0, tolerance ) ;
	}

	template< typename F >
	bool check_if_odd( F a, F tolerance = 0.0 ) {
		return floats_are_equal( std::fmod( a, 2.0 ), 1.0, tolerance ) ;
	}

	template< typename F >
	F round_to_nearest_integer( F a ) {
		F result = std::floor( a ) ;
		if( a - result > 0.5 ) {
			result += 1.0 ;
		}
		return result ;
	}

	//
	// log-sum-exp implementations
	//

	namespace impl {
		template< typename Iterator, typename F >
		F log_sum_exp( Iterator begin, Iterator const& end, F const max_value ) {
			F running_sum = 0.0 ;
			for( ; begin != end; ++begin ) {
				running_sum += std::exp( *begin - max_value ) ;
			}
			return max_value + std::log( running_sum ) ;
		}

		// class value_traits is a proxy for std::iterator_traits.
		// This is provided so that functions with similar signatures can be used
		// for value types (e.g. double) as well as iterators.
		template< typename Iterator >
		struct value_traits {
			typedef typename std::iterator_traits< Iterator >::value_type value_type ;
		} ;
	
		template<>
		struct value_traits< double >
		{
			typedef double value_type ;
		} ;
	}

	double log_sum_exp( double a, double b ) ;

	template< typename Iterator >
	typename impl::value_traits< Iterator >::value_type log_sum_exp( Iterator const& begin, Iterator const& end ) {
		typedef typename impl::value_traits< Iterator >::value_type F ;
		F max_value = -std::numeric_limits< F >::infinity() ;
		for( Iterator i = begin; i != end; ++i ) {
			max_value = std::max( max_value, *i ) ;
		}
		if( max_value == -std::numeric_limits< F >::infinity() ) {
			return -std::numeric_limits< F >::infinity() ;
		}
		else {
			return impl::log_sum_exp( begin, end, max_value ) ;
		}
	}

	template< typename F >
	double log_sum_exp( std::pair< F, F > const& values ) {
		return log_sum_exp( values.first, values.second ) ;
	}

	template< typename F >
	F log_sum_exp( std::vector< F > const& values ) {
		return log_sum_exp( &values[0], &values[0] + values.size() ) ;
	}

	template< typename F >
	double log_diff_exp( F a, F b ) {
		assert( a >= b ) ;
		if( b == -std::numeric_limits< F >::infinity() ) {
			return a ;
		}
		else {
			return a + std::log( 1.0 - std::exp( b - a )) ;
		}
	}

	// Calculate log( exp(value1) - exp(value2) ).
	template < typename F >
	double log_diff_exp( std::pair< F, F > const& values ) {
		return log_diff_exp( values.first, values.second ) ;
	}
}

#endif

