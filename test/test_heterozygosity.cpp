
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

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
#include "test_case.hpp"

namespace {
	std::map< std::string, double > get_data() {
		std::map< std::string, double > data ;
		data[ "S1 S1 0 A G 0 1 0 0 1 0 0 1 0 0 1 0 0 1 0" ] = 1.0 ;
		data[ "S1 S1 0 A G 1 0 0 0 1 0 0 1 0 0 1 0 0 1 0" ] = 4.0 / 5.0 ;
		data[ "S1 S1 0 A G 1 0 0 1 0 0 0 1 0 0 1 0 0 1 0" ] = 3.0 / 5.0 ;
		data[ "S1 S1 0 A G 1 0 0 1 0 0 1 0 0 0 1 0 0 1 0" ] = 2.0 / 5.0 ;
		data[ "S1 S1 0 A G 1 0 0 1 0 0 1 0 0 1 0 0 0 1 0" ] = 1.0 / 5.0 ;
		data[ "S1 S1 0 A G 1 0 0 1 0 0 1 0 0 1 0 0 1 0 0" ] = 0.0 ;
                 
		data[ "S1 S1 0 A G 0 0 1 0 1 0 0 1 0 0 1 0 0 1 0" ] = 4.0/5.0 ;
		data[ "S1 S1 0 A G 0 0 1 0 0 1 0 1 0 0 1 0 0 1 0" ] = 3.0/5.0 ;
		data[ "S1 S1 0 A G 0 0 1 0 0 1 0 0 1 0 1 0 0 1 0" ] = 2.0/5.0 ;
		data[ "S1 S1 0 A G 0 0 1 0 0 1 0 0 1 0 0 1 0 1 0" ] = 1.0/5.0 ;
		data[ "S1 S1 0 A G 0 0 1 0 0 1 0 0 1 0 0 1 0 0 1" ] = 0.0 ;
                 
		data[ "S1 S1 0 A G 0 0.5 0 0 0.2 0 0 0.9 0 0 0.1 0 0 0.2 0" ] = 1.0 ;
		data[ "S1 S1 0 A G 0.5 0.0 0 0 0.2 0 0 0.9 0 0 0.1 0 0 0.2 0" ] = 1.4 / 1.9 ;

		data[ "SA1 rs001 10000000 A G 0 0 0 0 0 0 0 0 0 0 0 0.5721 0 0.0207 0.9792" ] = 0.0207 / (0.0207 + 0.5721 + 0.9792) ;
		data[ "SA2 rs002 10010000 A G 0 0 1 0 1 0 1 0 0 0 1 0 1 0 0" ] = 2.0 / 5.0 ;
		data[ "SA3 rs003 10020000 C T 1 0 0 0 1 0 0 0 1 0 0.9967 0 0 0 1" ] = 1.9967 / 4.9967 ;
		data[ "SA4 rs004 10030000 G T 1 0 0 0 1 0 0 0 1 0 1 0 0 0 1" ] = 2.0 / 5.0 ;
		data[ "SA5 rs005 10040000 C G 0 0 .5 0 .5 0 .5 0 0 0 .5 0 .5 0 0" ] = 1.0 / 2.5 ;
		data[ "SA7 rs007 10060000 G T 0.5 0.5 0 0 0.5 0 1 0 0 0.5 0 0 0.5 0 0" ] = 1.0 / 3.5 ;


		return data ;
	}
}

AUTO_TEST_CASE( test_heterozygosity ) {
	std::map< std::string, double >
		data = get_data() ;
	
	std::map< std::string, double >::const_iterator
		i( data.begin() ),
		end_i( data.end() ) ;

	int count = 0 ;
	double epsilon = 0.000001 ;

	GenotypeAssayStatistics row_statistics ;
	std::auto_ptr< GenotypeAssayStatistic > statistic( new HeterozygosityStatistic ) ;
	statistic->set_precision( 5 ) ;
	row_statistics.add_statistic( "test stat", statistic ) ;
	row_statistics.format_column_headers( std::cout ) << "\n" ;

	for( ; i != end_i; ++i ) {
		std::istringstream inStream( i->first ) ;
		InternalStorageGenRow row ;
		inStream >> row ;
		row_statistics.process( row.begin_genotype_proportions(), row.end_genotype_proportions() ) ;
		std::cout << "row " << std::setw(3) << ++count << " : " << row_statistics << "\n" ;
		
		TEST_ASSERT( floats_are_equal_to_within_epsilon( row_statistics.get_value< double >( "test stat" ), i->second, epsilon )) ;
	}
}

#ifndef HAVE_BOOST_UNIT_TEST_FRAMEWORK
	int main( int argc, char** argv ) {
		test_heterozygosity() ;
	}
#endif
