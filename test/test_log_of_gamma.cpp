
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "gamma.hpp"
#include "floating_point_utils.hpp"
#include "OptionProcessor.hpp"

double naive_log_of_factorial( double x ) {
	assert( x > 0.0 ) ;
	if( x > 2.0 ) {
		return std::log( x ) + naive_log_of_factorial( x - 1.0 ) ;
	}
	else {
		return std::log( x ) ;
	}
}

void test_log_of_gamma( double tolerance ) {
	for( long x = 2; x < 10000; ++x ) {
		double a = naive_log_of_factorial( static_cast< double >( x - 1.0 )) ;
		double b = log_of_gamma_function( static_cast< double >( x )) ;
		assert( floats_are_equal( a, b, tolerance )) ;
	}
}

int main( int argc, char** argv ) {
	OptionProcessor options ;
    try {
		options[ "--floating_point_tolerance" ]
			.set_description( "Tolerance to use in floating point comparisons" )
			.set_takes_value()
			.set_default_value( 0.0000000001 ) ;

		options.process( argc, argv ) ;
    }
    catch( ArgumentProcessingException const& exception ) {
        std::cerr << "!! Error: " << exception.message() << ".\n";
        std::cerr << "Usage: gen-select [options]\n"
                << options
                << "\n" ;
        return -1 ;
    }

	double tolerance = options.get_value<double>( "--floating_point_tolerance" ) ;
	
	test_log_of_gamma( tolerance ) ;

	return 0 ;
}
