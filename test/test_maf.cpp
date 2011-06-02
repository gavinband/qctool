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
#include "test_case.hpp"

namespace {
	std::map< std::string, double > get_data() {
		std::map< std::string, double > data ;
		data[ "SA1 rs001 10000000 A G 0 0 0 0 0 0 0 0 0 0 0 0.5721 0 0.0207 0.9792" ] = .00658396946564885496 ;
		data[ "SA2 rs002 10010000 A G 0 0 1 0 1 0 1 0 0 0 1 0 1 0 0" ] = 0.4 ;
		data[ "SA3 rs003 10020000 C T 1 0 0 0 1 0 0 0 1 0 0.9967 0 0 0 1" ] = 0.39993395651123141273 ;
		data[ "SA4 rs004 10030000 G T 1 0 0 0 1 0 0 0 1 0 1 0 0 0 1" ] = 0.4 ;
		data[ "SA5 rs005 10040000 C G 0 0 1 0 1 0 1 0 0 0 1 0 1 0 0" ] = 0.4 ;
		data[ "SA6 rs006 10050000 A G 0 1 0 0 1 0 1 0 0 0 1 0 1 0 0" ] = 0.3 ;
		data[ "SA7 rs007 10060000 G T 1 0 0 0 1 0 1 0 0 1 0 0 1 0 0" ] = 0.1 ;
		data[ "SA8 rs008 10070000 G T 1 0 0 0 1 0 0 0 1 0 1 0 0 0 1" ] = 0.4 ;
		data[ "SA9 rs009 10080000 G T 0 0 1 0 1 0 0 0 1 0 1 0 0 0 1" ] = 0.2 ;
		data[ "SA10 rs010 10090000 A G 0 0 1 0 1 0 1 0 0 0 1 0 1 0 0" ] = 0.4 ;
		data[ "SA11 rs011 10100000 C T 0 0 0 0 1 0 1 0 0 1 0 0 0.67 0 0.23" ] = .18717948717948717948 ;

		return data ;
	}
}

// Test that we can read in InternalStorageGenRows, output them, and that this gives the same results.
AUTO_TEST_CASE( test_maf ) {
	std::map< std::string, double >
		data = get_data() ;
	
	std::map< std::string, double >::const_iterator
		i( data.begin() ),
		end_i( data.end() ) ;

	int count = 0 ;

	double epsilon = 0.000001 ;

	GenRowStatistics row_statistics ;
	std::auto_ptr< GenotypeAssayStatistic > statistic( new AlleleProportionStatistic( AlleleProportionStatistic::minor_allele ) ) ;
	statistic->set_precision( 5 ) ;
	row_statistics.add_statistic( "MAF", statistic ) ;
	row_statistics.format_column_headers( std::cout ) << "\n" ;

	for( ; i != end_i; ++i ) {
		std::istringstream inStream( i->first ) ;
		InternalStorageGenRow row ;
		inStream >> row ;
		row_statistics.process( row ) ;
		std::cout << "row " << std::setw(3) << ++count << " : " << row_statistics << "\n" ;
		
		TEST_ASSERT( floats_are_equal_to_within_epsilon( row_statistics.get_value< double >( "MAF" ), i->second, epsilon )) ;
	}
}

