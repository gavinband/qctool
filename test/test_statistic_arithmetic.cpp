#include <iostream>
#include <iomanip>
#include <sstream>
#include <cassert>
#include <map>
#include <string>
#include "GenRow.hpp"
#include "GenRowStatistics.hpp"
#include "GenotypeAssayStatisticFactory.hpp"
#include "floating_point_utils.hpp"
#include "test_case.hpp"

namespace {
	std::map< std::string, double > get_product_test_data() {
		std::map< std::string, double > data ;
		// List of test configurations
		// The number on the right is the product of:
		// 1. the maximum possible probability (under the assumption that genotypes are drawn independently
		// according to fixed probabilities p0, p1, p2 summing to 1) of observing the data
		// AND
		// 2. the minor allele frequency
		//                            |  1  |  2  |  3  |  4  |  5  |  6  |  7  |  8  |
		data ["SA01 rs001 10000001 A G 0 1 0 1 0 0 0 0 1 0 0 1 0 0 1"] =  std::log( (1.0/5.0) * (1.0/5.0) * std::pow(3.0/5.0, 3.0) ) * ( 3. / 10.0 );
		data ["SA01 rs001 10000001 A G 0 1 0 0 1 0 0 1 0 0 0 1 0 0 1"] =  std::log( (1.0) * std::pow(3.0/5.0, 3.0) * std::pow(2.0/5.0, 2.0) ) * (3. / 10.0);
		data ["SA01 rs001 10000001 A G 1 0 0 1 0 0 0 1 0 0 1 0 0 1 0 0 1 0 0 0 1 0 0 1"] = std::log( std::pow( 0.25, 2.0 ) * std::pow(0.5, 4.0 ) * std::pow( 0.25, 2.0 ) ) * ( 8. / 16.);
		data ["SA01 rs001 10000001 A G 1 0 0 1 0 0 1 0 0 0 1 0 0 1 0 0 1 0 0 0 1 0 0 1"] = std::log( std::pow( 3./8., 3.0 ) * std::pow( 3./8., 3.0 ) * std::pow( 2./8., 2.0 ) ) * ( 7./16.);
		data ["SA01 rs001 10000001 A G 1 0 0 1 0 0 1 0 0 1 0 0 0 1 0 0 1 0 0 0 1 0 0 1"] = std::log( std::pow( 4./8., 4.0 ) * std::pow( 2./8., 2.0 ) * std::pow( 2./8., 2.0 ) ) * ( 6./16.);

		return data ;
	}

	std::map< std::string, double > get_difference_test_data() {
		std::map< std::string, double > data ;
		// List of test configurations
		// The number on the right is the product of:
		// 1. the maximum possible probability (under the assumption that genotypes are drawn independently
		// according to fixed probabilities p0, p1, p2 summing to 1) of observing the data
		// AND
		// 2. the minor allele frequency
		//                            |  1  |  2  |  3  |  4  |  5  |  6  |  7  |  8  |
		data ["SA01 rs001 10000001 A G 0 1 0 1 0 0 0 0 1 0 0 1 0 0 1"] =  std::log( (1.0/5.0) * (1.0/5.0) * std::pow(3.0/5.0, 3.0) ) - ( 3. / 10.0 );
		data ["SA01 rs001 10000001 A G 0 1 0 0 1 0 0 1 0 0 0 1 0 0 1"] =  std::log( (1.0) * std::pow(3.0/5.0, 3.0) * std::pow(2.0/5.0, 2.0) ) - (3. / 10.0);
		data ["SA01 rs001 10000001 A G 1 0 0 1 0 0 0 1 0 0 1 0 0 1 0 0 1 0 0 0 1 0 0 1"] = std::log( std::pow( 0.25, 2.0 ) * std::pow(0.5, 4.0 ) * std::pow( 0.25, 2.0 ) ) - ( 8. / 16.);
		data ["SA01 rs001 10000001 A G 1 0 0 1 0 0 1 0 0 0 1 0 0 1 0 0 1 0 0 0 1 0 0 1"] = std::log( std::pow( 3./8., 3.0 ) * std::pow( 3./8., 3.0 ) * std::pow( 2./8., 2.0 ) ) - ( 7./16.);
		data ["SA01 rs001 10000001 A G 1 0 0 1 0 0 1 0 0 1 0 0 0 1 0 0 1 0 0 0 1 0 0 1"] = std::log( std::pow( 4./8., 4.0 ) * std::pow( 2./8., 2.0 ) * std::pow( 2./8., 2.0 ) ) - ( 6./16.);

		return data ;
	}

	std::map< std::string, double > get_ratio_test_data() {
		std::map< std::string, double > data ;
		// List of test configurations
		// The number on the right is the product of:
		// 1. the maximum possible probability (under the assumption that genotypes are drawn independently
		// according to fixed probabilities p0, p1, p2 summing to 1) of observing the data
		// AND
		// 2. the minor allele frequency
		//                            |  1  |  2  |  3  |  4  |  5  |  6  |  7  |  8  |
		data ["SA01 rs001 10000001 A G 0 1 0 1 0 0 0 0 1 0 0 1 0 0 1"] =  std::log( (1.0/5.0) * (1.0/5.0) * std::pow(3.0/5.0, 3.0) ) / ( 0.3 );
		data ["SA01 rs001 10000001 A G 0 1 0 0 1 0 0 1 0 0 0 1 0 0 1"] =  std::log( (1.0) * std::pow(3.0/5.0, 3.0) * std::pow(2.0/5.0, 2.0) ) / (  0.3 );
		data ["SA01 rs001 10000001 A G 1 0 0 1 0 0 0 1 0 0 1 0 0 1 0 0 1 0 0 0 1 0 0 1"] = std::log( std::pow( 0.25, 2.0 ) * std::pow(0.5, 4.0 ) * std::pow( 0.25, 2.0 ) ) / ( 0.5 );
		data ["SA01 rs001 10000001 A G 1 0 0 1 0 0 1 0 0 0 1 0 0 1 0 0 1 0 0 0 1 0 0 1"] = std::log( std::pow( 3./8., 3.0 ) * std::pow( 3./8., 3.0 ) * std::pow( 2./8., 2.0 ) ) / ( 7./16.);
		data ["SA01 rs001 10000001 A G 1 0 0 1 0 0 1 0 0 1 0 0 0 1 0 0 1 0 0 0 1 0 0 1"] = std::log( std::pow( 4./8., 4.0 ) * std::pow( 2./8., 2.0 ) * std::pow( 2./8., 2.0 ) ) / ( 6./16.);

		return data ;
	}

	std::map< std::string, double > get_sum_test_data() {
		std::map< std::string, double > data ;
		// List of test configurations
		// The number on the right is the product of:
		// 1. the maximum possible probability (under the assumption that genotypes are drawn independently
		// according to fixed probabilities p0, p1, p2 summing to 1) of observing the data
		// AND
		// 2. the minor allele frequency
		//                            |  1  |  2  |  3  |  4  |  5  |  6  |  7  |  8  |
		data ["SA01 rs001 10000001 A G 0 1 0 1 0 0 0 0 1 0 0 1 0 0 1"] =  std::log( (1.0/5.0) * (1.0/5.0) * std::pow(3.0/5.0, 3.0) ) + ( 3. / 10.0 );
		data ["SA01 rs001 10000001 A G 0 1 0 0 1 0 0 1 0 0 0 1 0 0 1"] =  std::log( (1.0) * std::pow(3.0/5.0, 3.0) * std::pow(2.0/5.0, 2.0) ) + (3. / 10.0);
		data ["SA01 rs001 10000001 A G 1 0 0 1 0 0 0 1 0 0 1 0 0 1 0 0 1 0 0 0 1 0 0 1"] = std::log( std::pow( 0.25, 2.0 ) * std::pow(0.5, 4.0 ) * std::pow( 0.25, 2.0 ) ) + ( 8. / 16.);
		data ["SA01 rs001 10000001 A G 1 0 0 1 0 0 1 0 0 0 1 0 0 1 0 0 1 0 0 0 1 0 0 1"] = std::log( std::pow( 3./8., 3.0 ) * std::pow( 3./8., 3.0 ) * std::pow( 2./8., 2.0 ) ) + ( 7./16.);
		data ["SA01 rs001 10000001 A G 1 0 0 1 0 0 1 0 0 1 0 0 0 1 0 0 1 0 0 0 1 0 0 1"] = std::log( std::pow( 4./8., 4.0 ) * std::pow( 2./8., 2.0 ) * std::pow( 2./8., 2.0 ) ) + ( 6./16.);

		return data ;
	}

	void do_test( std::string statistic_spec, std::map< std::string, double > const& data ) {
		std::map< std::string, double >::const_iterator
			i( data.begin() ),
			end_i( data.end() ) ;

		int count = 0 ;

		double epsilon = 0.1 ;

		GenRowStatistics row_statistics ;
		row_statistics.add_statistic( "MLIG", GenotypeAssayStatisticFactory::create_statistic( "MLIG" )) ;
		row_statistics.add_statistic( "MAF", GenotypeAssayStatisticFactory::create_statistic( "MAF" )) ;
		row_statistics.add_statistic( statistic_spec, GenotypeAssayStatisticFactory::create_statistic( statistic_spec )) ;
		std::cout << "          " ;
		row_statistics.format_column_headers( std::cout ) << "\n" ;

		for( ; i != end_i; ++i ) {
			std::istringstream inStream( i->first ) ;
			InternalStorageGenRow row ;
			inStream >> row ;
			row_statistics.process( row ) ;
			std::cout << "row " << std::setw(3) << ++count << " : " << row_statistics << " -- expected: " << i->second << ".\n" ;

			TEST_ASSERT( floats_are_equal_to_within_epsilon( row_statistics.get_value< double >( statistic_spec ), i->second, epsilon )) ;
		}
	}
}

AUTO_TEST_CASE( test_statistic_ratio ) {
	do_test( "MLIG / MAF", get_ratio_test_data() ) ;
}	

AUTO_TEST_CASE( test_statistic_sum ) {
	do_test( "MLIG + MAF", get_sum_test_data() ) ;
}	

AUTO_TEST_CASE( test_statistic_product ) {
	do_test( "MLIG * MAF", get_product_test_data() ) ;
}	

AUTO_TEST_CASE( test_statistic_difference ) {
	do_test( "MLIG - MAF", get_difference_test_data() ) ;
}	
