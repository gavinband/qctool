#include <iostream>
#include <iomanip>
#include <sstream>
#include <cassert>
#include <map>
#include <string>
#include "GenRow.hpp"
#include "GenRowStatistics.hpp"
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

std::map< std::string, std::string > get_data() {
	std::map< std::string, std::string > data ;
	data[ "S1 S1 0 A G 1 0 0"  ] = "G" ;
	data[ "S1 S1 0 A G 0 0 1"  ] = "A" ;
	data[ "S1 S1 0 A G 0.8 0.2 0"  ] = "G" ;
	data[ "S1 S1 0 A G 0 0.2 0.8"  ] = "A" ;
	
	return data ;
}

AUTO_TEST_CASE( test_alleles ) {
	std::map< std::string, std::string >
		data = get_data() ;
	
	std::map< std::string, std::string >::const_iterator
		i( data.begin() ),
		end_i( data.end() ) ;

	int count = 0 ;

	GenRowStatistics row_statistics ;
	std::auto_ptr< GenotypeAssayStatistic > minor_allele_statistic( new GenRowAllele( GenRowAllele::minor ) ) ;
	std::auto_ptr< GenotypeAssayStatistic > major_allele_statistic( new GenRowAllele( GenRowAllele::major ) ) ;
	row_statistics.add_statistic( "minor allele", minor_allele_statistic ) ;
	row_statistics.add_statistic( "major allele", major_allele_statistic ) ;
	row_statistics.format_column_headers( std::cout ) << "\n" ;

	for( ; i != end_i; ++i ) {
		std::istringstream inStream( i->first ) ;
		InternalStorageGenRow row ;
		inStream >> row ;
		row_statistics.process( row ) ;
		std::cout << "row " << std::setw(3) << ++count << " : " << row_statistics << "\n" ;
		TEST_ASSERT( row_statistics.get_value< std::string >( "minor allele" ) == i->second ) ;
		TEST_ASSERT( row_statistics.get_value< std::string >( "major allele" ) == ((i->second == "A") ? std::string("G") : std::string("A") )) ;
	}
}

#ifndef HAVE_BOOST_UNIT_TEST
	int main( int argc, char** argv ) {
		test_alleles() ;
	}
#endif
