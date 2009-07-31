#include <iostream>
#include <sstream>
#include <cassert>
#include <boost/bind.hpp>
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
#include "GenRow.hpp"
#include "endianness_utils.hpp"
#include "stdint.h"
#include "bgen.hpp"

namespace data {
	std::string data =
		"SA1 rs001 10000000 A G 0 0 0 0 0 0 0 0 0 0 0 0.5721 0 0.0207 0.9792\n"
		"SA2 rs002 10010000 A G 0 0 1 0 1 0 1 0 0 0 1 0 1 0 0\n"
		"SA3 rs003 10020000 C T 1 0 0 0 1 0 0 0 1 0 0.9967 0 0 0 1\n"
		"SA4 rs004 10030000 G T 1 0 0 0 1 0 0 0 1 0 1 0 0 0 1\n"
		"SA5 rs005 10040000 C G 0 0 1 0 1 0 1 0 0 0 1 0 1 0 0\n"
		"SA6 rs006 10050000 A G 0 1 0 0 1 0 1 0 0 0 1 0 1 0 0\n"
		"SA7 rs007 10060000 G T 1 0 0 0 1 0 1 0 0 1 0 0 1 0 0\n"
		"SA8 rs008 10070000 G T 1 0 0 0 1 0 0 0 1 0 1 0 0 0 1\n"
		"SA9 rs009 10080000 G T 0 0 1 0 1 0 0 0 1 0 1 0 0 0 1\n"
		"SA10 rs010 10090000 A G 0 0 1 0 1 0 1 0 0 0 1 0 1 0 0\n"
		"SA11 rs011 10100000 C T 0 0 0 0 1 0 1 0 0 1 0 0 0.67 0 0.23\n" ;

	std::string data_missing_0th_and_3rd_samples =
		"SA1 rs001 10000000 A G 0 0 0 0 0 0 0 0.0207 0.9792\n"
		"SA2 rs002 10010000 A G 0 1 0 1 0 0 1 0 0\n"
		"SA3 rs003 10020000 C T 0 1 0 0 0 1 0 0 1\n"
		"SA4 rs004 10030000 G T 0 1 0 0 0 1 0 0 1\n"
		"SA5 rs005 10040000 C G 0 1 0 1 0 0 1 0 0\n"
		"SA6 rs006 10050000 A G 0 1 0 1 0 0 1 0 0\n"
		"SA7 rs007 10060000 G T 0 1 0 1 0 0 1 0 0\n"
		"SA8 rs008 10070000 G T 0 1 0 0 0 1 0 0 1\n"
		"SA9 rs009 10080000 G T 0 1 0 0 0 1 0 0 1\n"
		"SA10 rs010 10090000 A G 0 1 0 1 0 0 1 0 0\n"
		"SA11 rs011 10100000 C T 0 1 0 1 0 0 0.67 0 0.23\n" ;

	std::string data_missing_3rd_and_4th_samples =
		"SA1 rs001 10000000 A G 0 0 0 0 0 0 0 0 0\n"
		"SA2 rs002 10010000 A G 0 0 1 0 1 0 1 0 0\n"
		"SA3 rs003 10020000 C T 1 0 0 0 1 0 0 0 1\n"
		"SA4 rs004 10030000 G T 1 0 0 0 1 0 0 0 1\n"
		"SA5 rs005 10040000 C G 0 0 1 0 1 0 1 0 0\n"
		"SA6 rs006 10050000 A G 0 1 0 0 1 0 1 0 0\n"
		"SA7 rs007 10060000 G T 1 0 0 0 1 0 1 0 0\n"
		"SA8 rs008 10070000 G T 1 0 0 0 1 0 0 0 1\n"
		"SA9 rs009 10080000 G T 0 0 1 0 1 0 0 0 1\n"
		"SA10 rs010 10090000 A G 0 0 1 0 1 0 1 0 0\n"
		"SA11 rs011 10100000 C T 0 0 0 0 1 0 1 0 0\n" ;
}

// Test that we can read in InternalStorageGenRows, output them, and that this gives the same results.
AUTO_TEST_CASE( test_genrow_sample_filtering ) {
	std::string const& data = data::data ;
	std::string const& filtered_data_1 = data::data_missing_0th_and_3rd_samples ;
	std::string const& filtered_data_2 = data::data_missing_3rd_and_4th_samples ;

	std::istringstream inStream( data ) ;
	std::istringstream filteredInStream1( filtered_data_1 ) ;
	std::istringstream filteredInStream2( filtered_data_2 ) ;

	InternalStorageGenRow row, row2, row3, row4 ;
	int count = 0;
	
	std::vector< std::size_t > zero_and_three(2) ;
	zero_and_three[0] = 0 ;
	zero_and_three[1] = 3 ;
	std::vector< std::size_t > three_and_four(2) ;
	three_and_four[0] = 3 ;
	three_and_four[1] = 4 ;

	while( inStream >> row && filteredInStream1 >> row2 && filteredInStream2 >> row4 ) {
		++count ;
		row3 = row ;
		row.filter_out_samples_with_indices( zero_and_three ) ;
		TEST_ASSERT( row == row2 ) ;
		row3.filter_out_samples_with_indices( three_and_four ) ;
		TEST_ASSERT( row3 == row4 ) ;
	}
	
	TEST_ASSERT( count == 11 ) ;
}

#ifndef HAVE_BOOST_UNIT_TEST

int main( int argc, char** argv ) {
	test_genrow_sample_filtering() ;
}

#endif

