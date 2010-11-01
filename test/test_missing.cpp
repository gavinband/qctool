#include <iostream>
#include <iomanip>
#include <sstream>
#include <cassert>
#include <map>
#include <string>
#include "GenRow.hpp"
#include "GenRowStatistics.hpp"
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
	data[ "SA1 rs001 10000000 A G" ] = std::numeric_limits< double >::quiet_NaN() ;
	data[ "SA1 rs001 10000000 A G 1 0 0" ] = 0 ;
	data[ "SA1 rs001 10000000 A G 0 1 0" ] = 0 ;
	data[ "SA1 rs001 10000000 A G 0 0 1" ] = 0 ;
	data[ "SA1 rs001 10000000 A G 0 0 0 1 0 0" ] = 0.5 ;
	data[ "SA1 rs001 10000000 A G 0 0 0 0 1 0" ] = 0.5 ;
	data[ "SA1 rs001 10000000 A G 0 0 0 0 0 1" ] = 0.5 ;
	data[ "SA1 rs001 10000000 A G 0 0 0 0 0 0" ] = 1 ;
	data[ "SA1 rs001 10000000 A G 0 0 0 0 0 0 0 0 0 0 0 0.5721 0 0.0207 0.9792" ] = .6856 ;
	data[ "SA2 rs002 10010000 A G 0 0 1 0 1 0 1 0 0 0 1 0 1 0 0" ] = 0.0 ;
	data[ "SA3 rs003 10020000 C T 1 0 0 0 1 0 0 0 1 0 0.9967 0 0 0 1" ] = .00066 ;
	data[ "SA4 rs004 10030000 G T 1 0 0 0 1 0 0 0 1 0 1 0 0 0 1" ] = 0 ;
	data[ "SA5 rs005 10040000 C G 0 0 .5 0 .5 0 .5 0 0 0 .5 0 .5 0 0" ] = 0.5 ;
	data[ "SA6 rs006 10050000 A G 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0" ] = 1 ;
	data[ "SA7 rs007 10060000 G T 0.5 0.5 0 0 0.5 0 1 0 0 0.5 0 0 0.5 0 0" ] = 1.5 / 5.0 ;

	return data ;
}

bool is_NaN( double value ) {
	return value != value ;
}

// Test that we can read in InternalStorageGenRows, output them, and that this gives the same results.
AUTO_TEST_CASE( test_missing ) {
	std::map< std::string, double >
		data = get_data() ;
	
	std::map< std::string, double >::const_iterator
		i( data.begin() ),
		end_i( data.end() ) ;

	int count = 0 ;

	double epsilon = 0.000001 ;

	GenRowStatistics row_statistics ;
	std::auto_ptr< GenotypeAssayStatistic > statistic( new MissingDataProportionStatistic ) ;
	statistic->set_precision( 5 ) ;
	row_statistics.add_statistic( "missing", statistic ) ;
	row_statistics.format_column_headers( std::cout ) << "\n" ;

	for( ; i != end_i; ++i ) {
		std::istringstream inStream( i->first ) ;
		InternalStorageGenRow row ;
		inStream >> row ;
		row_statistics.process( row ) ;
		std::cout << "row " << std::setw(3) << ++count << " : " << row_statistics << "\n" ;
		
		if( is_NaN( i->second )) {
			TEST_ASSERT( is_NaN( row_statistics.get_value< double >( "missing" ) )) ;
		}
		else {
			TEST_ASSERT( floats_are_equal_to_within_epsilon( row_statistics.get_value< double >( "missing" ), i->second, epsilon )) ;
		}
	}
}

#ifndef HAVE_BOOST_UNIT_TEST
	int main( int argc, char** argv ) {
		test_missing() ;
	}
#endif
