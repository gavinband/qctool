
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include "LikelihoodRatioTestStatistic.hpp"
#include "floating_point_utils.hpp"
#include "OptionProcessor.hpp"
#include "GenotypeAssayStatistics.hpp"
#include "GenRowStatistics.hpp"
#include "GenotypeAssayStatisticFactory.hpp"

std::map< std::string, double > get_test_data() {
	std::map< std::string, double > data ;
	// List of test configurations
	// The number on the right is the maximum possible probability (under the assumption that genotypes are drawn independently
	// according to fixed probabilities p0, p1, p2 summing to 1) of observing the data.
	//                            |  1  |  2  |  3  |  4  |  5  |  6  |  7  |  8  |
	data ["SA01 rs001 10000001 A G 0 1 0 1 0 0 0 0 1 0 0 1 0 0 1"] =  (1.0/5.0) * (1.0/5.0) * std::pow(3.0/5.0, 3.0);
	data ["SA01 rs001 10000001 A G 0 1 0 0 1 0 0 1 0 0 0 1 0 0 1"] =  (1.0) * std::pow(3.0/5.0, 3.0) * std::pow(2.0/5.0, 2.0);
	data ["SA01 rs001 10000001 A G 1 0 0 1 0 0 0 1 0 0 1 0 0 1 0 0 1 0 0 0 1 0 0 1"] = std::pow( 0.25, 2.0 ) * std::pow(0.5, 4.0 ) * std::pow( 0.25, 2.0 ) ;
	data ["SA01 rs001 10000001 A G 1 0 0 1 0 0 1 0 0 0 1 0 0 1 0 0 1 0 0 0 1 0 0 1"] = std::pow( 3./8., 3.0 ) * std::pow( 3./8., 3.0 ) * std::pow( 2./8., 2.0 ) ;
	data ["SA01 rs001 10000001 A G 1 0 0 1 0 0 1 0 0 1 0 0 0 1 0 0 1 0 0 0 1 0 0 1"] = std::pow( 4./8., 4.0 ) * std::pow( 2./8., 2.0 ) * std::pow( 2./8., 2.0 ) ;

	return data ;
}


std::map< std::string, double > get_test_data_for_HW() {
	std::map< std::string, double > data ;
	// List of test configurations
	// The number on the right is the maximum possible probability of observing the data, under the assumption
	// that the genotypes are drawn independently according to fixed probabilities p0, p1, p2 summing to 1,
	// which are in hardy-weinberg proportions.
	//                            |  1  |  2  |  3  |  4  |  5  |  6  |  7  |  8  |
	data ["SA01 rs001 10000001 A G 0 1 0 1 0 0 0 0 1 0 0 1 0 0 1"] =  std::pow( 0.3, 2.0 * 1.0 ) * std::pow( 2.*0.3*0.7, 1.0 ) * std::pow( 0.7, 2.0 * 3.0 ) ;
	data ["SA01 rs001 10000001 A G 0 1 0 0 1 0 0 1 0 0 0 1 0 0 1"] =  1.0 * std::pow( 2.*0.3*0.7, 3.0 ) * std::pow( 0.7, 2.*2. ) ;
	data ["SA01 rs001 10000001 A G 1 0 0 1 0 0 0 1 0 0 1 0 0 1 0 0 1 0 0 0 1 0 0 1"] = std::pow( (8./16.), 2.0 * 2.0 ) * std::pow( 2.0 * (8./16.) * (8./16.), 4.0 ) * std::pow( (8./16.), 2.0 * 2.0 ) ;
	data ["SA01 rs001 10000001 A G 1 0 0 1 0 0 1 0 0 0 1 0 0 1 0 0 1 0 0 0 1 0 0 1"] = std::pow( (9./16.), 2.0 * 3.0 ) * std::pow( 2.0 * (9./16.) * (7./16.), 3.0 ) * std::pow( (7./16.), 2.0 * 2.0 ) ;
	data ["SA01 rs001 10000001 A G 1 0 0 1 0 0 1 0 0 1 0 0 0 1 0 0 1 0 0 0 1 0 0 1"] = std::pow( (10./16.), 2.0 * 4.0 ) * std::pow( 2.0 * (10./16.) * (6./16.), 2.0 ) * std::pow( (6./16.), 2.0 * 2.0 ) ;

	return data ;
}

void test_maximum_likelihood( std::map< std::string, double > const& data, std::string const& statistic_spec, double const epsilon ) {
	std::map< std::string, double >::const_iterator
		i( data.begin() ),
		end_i( data.end() ) ;
		
	int count = 0 ;

	std::cout << "         " ;
	GenRowStatistics row_statistics ;
	row_statistics.add_statistic( statistic_spec, GenotypeAssayStatisticFactory::create_statistic( statistic_spec )) ;
	row_statistics.format_column_headers( std::cout ) << "\n" ;

	for( ; i != end_i; ++i ) {
		std::istringstream inStream( i->first ) ;
		++ count ;
		GenRow row ;
		inStream >> row ;
		row_statistics.process( row ) ;
		std::cout << "row " << std::setw(3) << count << ": " << row_statistics << "\n" ;
		
		assert( floats_are_equal_to_within_epsilon( row_statistics.get_statistic_value< double >( statistic_spec ), i->second, epsilon )) ;
	}
}


int main( int argc, char** argv ) {
	OptionProcessor options ;
    try {
		options[ "--epsilon" ]
			.set_description( "Epsilon to use when comparing floating point numbers" )
			.set_takes_value()
			.set_default_value( 0.00001 ) ;

		options.process( argc, argv ) ;
    }
    catch( ArgumentProcessingException const& exception ) {
        std::cerr << "!! Error: " << exception.message() << ".\n";
        std::cerr << "Usage: gen-select [options]\n"
                << options
                << "\n" ;
        return -1 ;
    }

	double epsilon = options.get_value<double>( "--epsilon" ) ;
	test_maximum_likelihood( get_test_data(), "MLIG", epsilon ) ;
	test_maximum_likelihood( get_test_data_for_HW(), "MLIGHW", epsilon ) ;

	return 0 ;
}
