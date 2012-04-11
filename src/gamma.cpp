
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


#include <vector>
#include <cmath>
#include <iostream>
#include <cassert>
#include "gamma.hpp"

namespace impl {
	// coefficients from GNU scientific library
	// G = 7
	double const coefficient_array1[] = {
		0.99999999999980993,
		676.5203681218851,
		-1259.1392167224028,
		771.32342877765313,
		-176.61502916214059,
		12.507343278686905,
		-0.13857109526572012,
		9.9843695780195716E-6,
		1.5056327351493116E-7
	} ;

	std::size_t const coefficient_array1_size = 9 ;

	// coefficients from http://home.att.net/~numericana/answer/info/godfrey.htm
	// G = 9
	double const coefficient_array2[] = {
		1.000000000000000174663,
		5716.400188274341379136,
		-14815.30426768413909044,
		14291.49277657478554025,
		-6348.160217641458813289,
		1301.608286058321874105,
		-108.1767053514369634679,
		2.605696505611755827729,
		-0.7423452510201416151527e-2,
		0.5384136432509564062961e-7,
		-0.4023533141268236372067e-8
	} ;

	std::size_t const coefficient_array2_size = 11 ;

	double const pi = 3.1415926535897932384626433832795 ;
	double const sqrt_of_2_pi = 2.506628274631000502415765284810 ;
	
	std::vector< double > log_of_factorial_table ;
}

double log_of_gamma_function( double x, double const * const coefficient_array, std::size_t const number_of_coefficients ) {
	// Implementation based on Lanczos approximation, see the following sources:
	// Numerical Recipes in C++, 3rd edition.
	// Wikipedia page for Lanczos approximation
	// http://home.att.net/~numericana/answer/info/godfrey.htm
	// http://en.literateprograms.org/Gamma_function_with_the_Lanczos_approximation_(Ruby)
	if( x < 0.5 ) {
		return std::log( impl::pi ) - std::log( std::sin( impl::pi * x ) * gamma( 1-x )) ;
	}
	else {
		std::size_t G = number_of_coefficients - 2;
		x -= 1.0; // Lanczos approximation is for \Gamma(x-1).
		double counting_x = x ;
		double tmp = coefficient_array[0] ;
		for( std::size_t i = 1; i < number_of_coefficients; ++i ) {
			tmp += coefficient_array[i] / (++counting_x) ;
		}
		double another_tmp = x + G + 0.5 ;
		return ((x+0.5) * std::log( another_tmp )) - another_tmp + std::log( impl::sqrt_of_2_pi * tmp );
	}
}

double log_of_factorial( unsigned long x ) {
	// Lookup the log-of-factorial value in the table.
	// Calculate if not yet calculated.
	std::size_t table_size = impl::log_of_factorial_table.size() ;
	
	if( table_size < x + 1 ) {
		impl::log_of_factorial_table.resize( x + 1, -1.0 ) ; // fill with -1s, indicating value not calculated.
	}

	double result = impl::log_of_factorial_table[ x ] ;
	if( result < 0.0 ) {
		// result not calculated yet
		result = log_of_gamma_function(x + 1) ;
	}
	impl::log_of_factorial_table[ x ] = result ;

	return result ;
}
