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
#include "RowCondition.hpp"
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

std::map< std::string, bool > get_data() {
	// rightmost value is whether the missing data proportion is > 0.5
	std::map< std::string, bool > data ;
	data[ "SA1 rs001 10000000 A G 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0" ] = true ;
	data[ "SA1 rs001 10000000 A G 0.5 0 0 0 0 0 0 0 0 0 0 0 0 0 0" ] = true ;
	data[ "SA1 rs001 10000000 A G 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0" ] = true ;
	data[ "SA1 rs001 10000000 A G 1 0 0 0.5 0 0 0 0 0 0 0 0 0 0 0" ] = true ;
	data[ "SA1 rs001 10000000 A G 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0" ] = true ;
	data[ "SA1 rs001 10000000 A G 1 0 0 1 0 0 0 0.4 0 0 0 0 0 0 0" ] = true ;
	data[ "SA1 rs001 10000000 A G 1 0 0 1 0 0 0 0.5 0 0 0 0 0 0 0" ] = false ;
	data[ "SA1 rs001 10000000 A G 0 0 0 0 0 0 0 0 0 0 0 0.5721 0 0.0207 0.9792" ] = true ;
	data[ "SA2 rs002 10010000 A G 0 0 1 0 1 0 1 0 0 0 1 0 1 0 0" ] = false ;
	data[ "SA3 rs003 10020000 C T 1 0 0 0 1 0 0 0 1 0 0.9967 0 0 0 1" ] = false ;
	data[ "SA4 rs004 10030000 G T 1 0 0 0 1 0 0 0 1 0 1 0 0 0 1" ] = false ;
	data[ "SA5 rs005 10040000 C G 0 0 .5 0 .5 0 .5 0 0 0 .5 0 .5 0 0" ] = false ;
	data[ "SA6 rs006 10050000 A G 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0" ] = true ;
	data[ "SA7 rs007 10060000 G T 0.5 0.5 0 0 0.5 0 1 0 0 0.5 0 0 0.5 0 0" ] = false ;

	return data ;
}

AUTO_TEST_CASE( test_less_than_greater_than ) {
	std::map< std::string, bool >
		data = get_data() ;
	
	std::map< std::string, bool >::const_iterator
		i( data.begin() ),
		end_i( data.end() ) ;

	int count = 0 ;

	GenRowStatistics row_statistics ;
	std::auto_ptr< GenotypeAssayStatistic > statistic( new MissingDataProportionStatistic ) ;
	statistic->set_precision( 5 ) ;
	row_statistics.add_statistic( "missing", statistic ) ;
	row_statistics.format_column_headers( std::cout ) << "\n" ;

	StatisticGreaterThan gt_condition( "missing", 0.5 ) ;
	StatisticLessThan lt_condition( "missing", 0.5 ) ;

	for( ; i != end_i; ++i ) {
		std::istringstream inStream( i->first ) ;
		InternalStorageGenRow row ;
		inStream >> row ;
		row_statistics.process( row ) ;
		std::cout << "row " << std::setw(3) << ++count << " : " << row_statistics << ": " << std::boolalpha << i->second << "\n" ;
		TEST_ASSERT( gt_condition.check_if_satisfied( row_statistics ) == i->second ) ;
		TEST_ASSERT(( row_statistics.get_value< double >( "missing" ) == 0.5 ) || ( lt_condition.check_if_satisfied( row_statistics ) != i->second )) ;
	}
}

#ifndef HAVE_BOOST_UNIT_TEST
	int main( int argc, char** argv ) {
		test_less_than_greater_than() ;
	}
#endif
