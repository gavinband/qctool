#ifndef __GTOOL__GAMMA_HPP__
#define __GTOOL__GAMMA_HPP__

#include <vector>
#include <cstddef>

namespace impl {
	extern double const coefficient_array1[] ;
	extern std::size_t const coefficient_array1_size ;
	extern double const coefficient_array2[] ;
	extern std::size_t const coefficient_array2_size ;
	extern double const pi ;
	extern double const sqrt_of_2_pi ;
	extern std::vector< double > log_of_factorial_table ;
}

double log_of_gamma_function( double x, double const* const coefficient_array = impl::coefficient_array2, std::size_t const number_of_coefficients = impl::coefficient_array2_size ) ;
double log_of_factorial( unsigned long x ) ;

#endif
