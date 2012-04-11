
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <sstream>
#include <cassert>
#include <boost/bind.hpp>
#include "../config.hpp"
#include "test_case.hpp"
#include "GenRow.hpp"
#include "genfile/endianness_utils.hpp"
#include "stdint.h"
#include "genfile/bgen.hpp"

namespace data {
	namespace {
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

		std::string data_missing_all_samples =
			"SA1 rs001 10000000 A G\n"
			"SA2 rs002 10010000 A G\n"
			"SA3 rs003 10020000 C T\n"
			"SA4 rs004 10030000 G T\n"
			"SA5 rs005 10040000 C G\n"
			"SA6 rs006 10050000 A G\n"
			"SA7 rs007 10060000 G T\n"
			"SA8 rs008 10070000 G T\n"
			"SA9 rs009 10080000 G T\n"
			"SA10 rs010 10090000 A G\n"
			"SA11 rs011 10100000 C T\n" ;
	}
}

// Test that we can read in InternalStorageGenRows, output them, and that this gives the same results.
AUTO_TEST_CASE( test_genrow_sample_filtering ) {
	std::string const& data = data::data ;
	std::string const& filtered_data_1 = data::data_missing_0th_and_3rd_samples ;
	std::string const& filtered_data_2 = data::data_missing_3rd_and_4th_samples ;
	std::string const& filtered_data_3 = data::data_missing_all_samples ;

	std::istringstream inStream( data ) ;
	std::istringstream filteredInStream1( filtered_data_1 ) ;
	std::istringstream filteredInStream2( filtered_data_2 ) ;
	std::istringstream filteredInStream3( filtered_data_3 ) ;

	InternalStorageGenRow row, row2, row3, row4, row5, row6 ;
	int count = 0;
	
	std::vector< std::size_t > zero_and_three(2) ;
	zero_and_three[0] = 0 ;
	zero_and_three[1] = 3 ;
	std::vector< std::size_t > three_and_four(2) ;
	three_and_four[0] = 3 ;
	three_and_four[1] = 4 ;
	std::vector< std::size_t > all_five(5) ;
	all_five[0] = 0 ;
	all_five[1] = 1 ;
	all_five[2] = 2 ;
	all_five[3] = 3 ;
	all_five[4] = 4 ;

	while( inStream >> row && filteredInStream1 >> row2 && filteredInStream2 >> row4 && filteredInStream3 >> row6 ) {
		++count ;
		row3 = row5 = row ;
		row.filter_out_samples_with_indices( zero_and_three ) ;
		TEST_ASSERT( row == row2 ) ;
		row3.filter_out_samples_with_indices( three_and_four ) ;
		TEST_ASSERT( row3 == row4 ) ;
		row5.filter_out_samples_with_indices( all_five ) ;
		TEST_ASSERT( row5 == row6 ) ;
	}
	
	TEST_ASSERT( count == 11 ) ;
}

