
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


#include <string>
#include <deque>
#include <iostream>
#include <algorithm>
#include <boost/math/distributions/chi_squared.hpp>
#include "genfile/string_utils.hpp"

int main() {
	std::cerr << "Welcome to inflation.cpp\n" ;
	std::cerr << "(c) Gavin Band (2013)\n" ;

	std::cerr << "inflation.cpp: reading..." ;
	float d ;
	std::string line ;
	std::deque< float > numbers ;
	std::size_t missing_count = 0 ;
	while( std::getline( std::cin, line )) {
		if( line.empty() || line.find_first_not_of( "0123456789.E+-" ) != std::string::npos ) {
			++missing_count ;
		} else {
			numbers.push_back( genfile::string_utils::to_repr< float >( line ) ) ;
		}
		
		if( numbers.size() % 1000000 == 0 ) {
			std::cerr << ( numbers.size() / 1000000 ) << ".." ;
		}
	}
	
	std::cerr << ( numbers.size() / 1000000 ) << ".\n" ;
	std::cerr << "inflation.cpp: read " << numbers.size() << " P-values from std::cin.\n" ;
	std::cerr << "inflation.cpp: converting to chi-squared statistics...\n" ;
	using namespace boost::math ;
	chi_squared_distribution< float > chi_square( 1 ) ;
	for( std::size_t i = 0; i < numbers.size(); ++i ) {
		std::cerr << numbers[i] << " -> " ;
		numbers[i] = quantile( complement( chi_square, numbers[i] )) ;
		std::cerr << numbers[i] << "\n" ;
	}
	
	std::cerr << "inflation.cpp: finding median...\n" ;

	std::cout << "{ \"N\"=" << numbers.size() << ", \"missing\"=" << missing_count ;
	// Compute basic lambda
	{
		std::size_t const mid = numbers.size() / 2 ;
		std::nth_element( numbers.begin(), numbers.begin() + mid, numbers.end() ) ;
		double const median = numbers[ mid ] ;
		std::cout << ", \"median\"=" << median << ", \"lambda\"=" << median / quantile( complement( chi_square, 0.5 )) ;
	}
	
	// Compute lambda after remove 1% of signal
	{
		std::size_t const one_percent = numbers.size() * 0.01 ;
		std::deque< float >::iterator one_percent_i = numbers.begin() + one_percent ;
		numbers.erase( numbers.begin(), one_percent_i ) ;
	}
	
	{
		std::size_t const mid = numbers.size() / 2 ;
		std::nth_element( numbers.begin(), numbers.begin() + mid, numbers.end() ) ;
		double const median = numbers[ mid ] ;
		std::cout << ", \"1pc_adjusted_median\"=" << median << ", \"1pc_adjusted_lambda\"=" << median / quantile( complement( chi_square, 0.5 )) ;
	}
	std::cout << " }\n" ;

	return 0 ;
}
