#include <iostream>
#include <iomanip>
#include <sstream>
#include <cassert>
#include <map>
#include <string>
#include "GenRow.hpp"
#include "GenRowStatistics.hpp"
#include "InformationStatistic.hpp"
#include "SimpleGenotypeAssayStatistics.hpp"
#include "floating_point_utils.hpp"
#include "flip_alleles.hpp"
#include "../config.hpp"
#if HAVE_BOOST_UNIT_TEST
	#define BOOST_AUTO_TEST_MAIN
	#include "boost/test/auto_unit_test.hpp"
	#define AUTO_TEST_CASE( param ) BOOST_AUTO_TEST_CASE(param)
	#define TEST_ASSERT( param ) BOOST_ASSERT( param )
#else
	#define AUTO_TEST_CASE( param ) void param()
	#define TEST_ASSERT( param ) assert( param )
#endif

std::map< std::string, double > get_data() {
	std::map< std::string, double > data ;
	// Miscellaneous tests
	data[ "SA1 rs001 10000000 A G" ] = 0.0 ;
	data[ "SA1 rs001 10000000 A G 0 0 0" ] = 0.0 ;
	data[ "SA1 rs001 10000000 A G 0 0 0 0 0 0" ] = 0.0 ;
	data[ "SA1 rs001 10000000 A G 0.33333333333 0.333333333333 0.33333333333" ] = 0 ;
	data[ "SA1 rs001 10000000 A G 0.5 0.5 0.0" ] = 0.33333333333333 ;
	data[ "SA1 rs001 10000000 A G 0.5 0.0 0.5" ] = 0 ;
	data[ "SA1 rs001 10000000 A G 0.0 0.5 0.5" ] = 0.33333333333333 ;
	// Information = 1 when there is no uncertainty.
	data[ "SA1 rs001 10000000 A G 1 0 0" ] = 1 ;
	data[ "SA1 rs001 10000000 A G 0 1 0" ] = 1 ;
	data[ "SA1 rs001 10000000 A G 0 0 1" ] = 1 ;
	data[ "SA1 rs001 10000000 A G 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0" ] = 1 ;
	data[ "SA1 rs001 10000000 A G 1 0 0 0 1 0 0 0 1 0 1 0 1 0 0" ] = 1 ;
	// Information = 0 when p_i = .25 .5 .25 for all i
	data[ "SA1 rs001 10000000 A G 0.25 0.5 0.25" ] = 0 ;
	data[ "SA1 rs001 10000000 A G 0.25 0.5 0.25 0.25 0.5 0.25" ] = 0 ;
	data[ "SA1 rs001 10000000 A G 0.25 0.5 0.25 0.25 0.5 0.25 0.25 0.5 0.25" ] = 0 ;
	data[ "SA1 rs001 10000000 A G 0.25 0.5 0.25 0.25 0.5 0.25 0.25 0.5 0.25 0.25 0.5 0.25 0.25 0.5 0.25 0.25 0.5 0.25" ] = 0 ;
	return data ;
}

std::vector< std::string > get_flipping_data() {
	std::vector< std::string > data ;
	// Tests with various amounts of missing data.
	// One individual
	data.push_back( "SA1 rs001 10000000 A G 1 0 0" ) ;
	data.push_back( "SA1 rs001 10000000 A G 0.9 0 0" ) ;
	data.push_back( "SA1 rs001 10000000 A G 0.5 0 0" ) ;
	data.push_back( "SA1 rs001 10000000 A G 0.1 0 0" ) ;
	data.push_back( "SA1 rs001 10000000 A G 0.0 0 0" ) ;

	data.push_back( "SA1 rs001 10000000 A G 0.5 0.5 0" ) ;
	data.push_back( "SA1 rs001 10000000 A G 0.5 0.2 0" ) ;
	data.push_back( "SA1 rs001 10000000 A G 0.5 0.1 0" ) ;

	// One individual
	data.push_back( "SA1 rs001 10000000 A G 0.5 0 0 0.5 0 0" ) ;
	data.push_back( "SA1 rs001 10000000 A G 0.5 0 0 0.5 0 0.5" ) ;
	data.push_back( "SA1 rs001 10000000 A G 0.5 0 0 0.0 0 0.5" ) ;
	data.push_back( "SA1 rs001 10000000 A G 0.5 0 0.5 0.0 0 0.5" ) ;

	return data ;
}

AUTO_TEST_CASE( test_information ) {
	std::cerr << "test_information()..." << "\n" ;
	std::map< std::string, double >
		data = get_data() ;
	
	std::map< std::string, double >::const_iterator
		i( data.begin() ),
		end_i( data.end() ) ;

	int count = 0 ;

	double epsilon = 0.000001 ;

	GenRowStatistics row_statistics ;
	std::auto_ptr< GenotypeAssayStatistic > statistic( new FillingInformationStatistic ) ;
	statistic->set_precision( 5 ) ;
	row_statistics.add_statistic( "Information", statistic ) ;
	row_statistics.format_column_headers( std::cout ) << "\n" ;

	for( ; i != end_i; ++i ) {
		std::istringstream inStream( i->first ) ;
		InternalStorageGenRow row ;
		inStream >> row ;
		row_statistics.process( row ) ;
		std::cout << "row " << std::setw(3) << ++count << " : " << row << " : " << row_statistics << "\n" ;
		
		TEST_ASSERT( floats_are_equal_to_within_epsilon( row_statistics.get_value< double >( "Information" ), i->second, epsilon )) ;
	}
}

AUTO_TEST_CASE( test_allele_flipping ) {
	std::cerr << "test_allele_flipping()..." << "\n" ;
	std::vector< std::string > data = get_flipping_data() ;
	
	GenRowStatistics row_statistics ;
	std::auto_ptr< GenotypeAssayStatistic > statistic( new FillingInformationStatistic ) ;
	statistic->set_precision( 5 ) ;
	row_statistics.add_statistic( "Information", statistic ) ;
	row_statistics.format_column_headers( std::cout ) << "\n" ;
	
	double epsilon = 0.000001 ;
	
	for( std::size_t i = 0; i < data.size(); ++i ) {
		std::istringstream inStream( data[i] ) ;
		InternalStorageGenRow row ;
		inStream >> row ;
		row_statistics.process( row ) ;
		double info = row_statistics.get_value< double >( "Information" ) ;
		std::cout << "row " << std::setw(3) << i << " : " << row << " : " << info << " : " ;
		flip_alleles( &row ) ;
		row_statistics.process( row ) ;
		double flipped_info = row_statistics.get_value< double >( "Information" ) ;
		std::cout << row << " : " << flipped_info << "\n" ;
		TEST_ASSERT( floats_are_equal_to_within_epsilon( info, flipped_info, epsilon )) ;
		
	}
}

#ifndef HAVE_BOOST_UNIT_TEST
	int main( int argc, char** argv ) {
		test_information() ;
		test_allele_flipping() ;
	}
#endif
