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
	data[ "SA1 rs001 10000000 A G 0 0 0" ] = 1.0 ;
	data[ "SA1 rs001 10000000 A G 0 0 0 0 0 0" ] = 1.0 ;
	data[ "SA1 rs001 10000000 A G 0.33333333333 0.333333333333 0.33333333333" ] = 0 ;
	data[ "SA1 rs001 10000000 A G 0.5 0.5 0.0" ] = 0.33333333333333 ;
	data[ "SA1 rs001 10000000 A G 0.5 0.0 0.5" ] = 0 ;
	data[ "SA1 rs001 10000000 A G 0.0 0.5 0.5" ] = 0.33333333333333 ;
	// Imnformation = 1 when there is no uncertainty.
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

AUTO_TEST_CASE( test_information ) {
	std::map< std::string, double >
		data = get_data() ;
	
	std::map< std::string, double >::const_iterator
		i( data.begin() ),
		end_i( data.end() ) ;

	int count = 0 ;

	double epsilon = 0.000001 ;

	GenRowStatistics row_statistics ;
	std::auto_ptr< GenotypeAssayStatistic > statistic( new InformationStatistic ) ;
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

#ifndef HAVE_BOOST_UNIT_TEST
	int main( int argc, char** argv ) {
		test_information() ;
	}
#endif
