
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <utility>
#include <map>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "test_case.hpp"
#include "statfile/DelimitedStatSource.hpp"
#include "statfile/FilteringStatSource.hpp"

BOOST_AUTO_TEST_SUITE( test_filtering_stat_source )

namespace globals {
	namespace {
		std::string data =
		"# This is a file\n"
		"# These lines are comments, which should be ignored.\n"
		"#\n"
		"index Column1 Column2 Column3 v0 v1 v2 v3 v4 v5 v6 v7 v8 v9\n"
		"0 Hello hello hello 0 0 0 0 0 0 0 0 0 0\n"
		"1 row1 A A 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n"
		"2 row2 B B 100.0 100000.254 576.22 -0.1 -100.001 100.001 1000000 4294967296 9007199254740992 -9007199254740992\n"
		"3 row3 C C 0.00001 9999.99999 0.00001 9999.99999 0.00001 9999.99999 0.00001 9999.99999 0.00001 9999.99999\n" ;

		std::size_t number_of_rows = 4 ;
	}
}

AUTO_TEST_CASE( test_equal_to ) {
	TEST_ASSERT( statfile::Constraint::equal_to( "hello" ).test( "hello" ) ) ;
	TEST_ASSERT( !statfile::Constraint::equal_to( "hello" ).test( "goodbye" ) ) ;
	TEST_ASSERT( !statfile::Constraint::equal_to( "hello" ).test( 5 ) ) ;

	TEST_ASSERT( statfile::Constraint::equal_to( 10 ).test( 10 ) ) ;
	TEST_ASSERT( statfile::Constraint::equal_to( 10 ).test( "10" ) ) ;
	TEST_ASSERT( statfile::Constraint::equal_to( 10 ).test( 10.0 ) ) ;
	TEST_ASSERT( statfile::Constraint::equal_to( 10.0 ).test( 10 ) ) ;

	TEST_ASSERT( !statfile::Constraint::equal_to( 10 ).test( 11 ) ) ;
	TEST_ASSERT( !statfile::Constraint::equal_to( 10 ).test( 10.1 ) ) ;
}

AUTO_TEST_CASE( test_less_than ) {
	TEST_ASSERT( statfile::Constraint::less_than( "abcd" ).test( "aacd" ) ) ;
	TEST_ASSERT( statfile::Constraint::less_than( "abcd" ).test( "abbd" ) ) ;
	TEST_ASSERT( statfile::Constraint::less_than( "abcd" ).test( "abcc" ) ) ;
	TEST_ASSERT( !statfile::Constraint::less_than( "abcd" ).test( "abcd" ) ) ;
	TEST_ASSERT( !statfile::Constraint::less_than( "abcd" ).test( "bbcd" ) ) ;
	TEST_ASSERT( !statfile::Constraint::less_than( "abcd" ).test( "accd" ) ) ;

	TEST_ASSERT( statfile::Constraint::less_than( 10 ).test( 9 ) ) ;
	TEST_ASSERT( statfile::Constraint::less_than( 10 ).test( "9" ) ) ;
	TEST_ASSERT( statfile::Constraint::less_than( 10 ).test( 9.0 ) ) ;

	TEST_ASSERT( !statfile::Constraint::less_than( 10 ).test( 10 ) ) ;
	TEST_ASSERT( !statfile::Constraint::less_than( 10 ).test( 10.0 ) ) ;
	TEST_ASSERT( !statfile::Constraint::less_than( 10 ).test( 10.0 ) ) ;
	TEST_ASSERT( !statfile::Constraint::less_than( 10 ).test( "10.0" ) ) ;

	TEST_ASSERT( !statfile::Constraint::less_than( 10 ).test( 11 ) ) ;
	TEST_ASSERT( !statfile::Constraint::less_than( 10 ).test( 10.1 ) ) ;
	TEST_ASSERT( !statfile::Constraint::less_than( 10 ).test( "11" ) ) ;
	TEST_ASSERT( !statfile::Constraint::less_than( 10 ).test( "10.1" ) ) ;
}

AUTO_TEST_CASE( test_less_than_or_equal_to ) {
	TEST_ASSERT( statfile::Constraint::less_than_or_equal_to( "abcd" ).test( "aacd" ) ) ;
	TEST_ASSERT( statfile::Constraint::less_than_or_equal_to( "abcd" ).test( "abbd" ) ) ;
	TEST_ASSERT( statfile::Constraint::less_than_or_equal_to( "abcd" ).test( "abcc" ) ) ;
	TEST_ASSERT( statfile::Constraint::less_than_or_equal_to( "abcd" ).test( "abcd" ) ) ;
	TEST_ASSERT( !statfile::Constraint::less_than_or_equal_to( "abcd" ).test( "bbcd" ) ) ;
	TEST_ASSERT( !statfile::Constraint::less_than_or_equal_to( "abcd" ).test( "accd" ) ) ;

	TEST_ASSERT( statfile::Constraint::less_than_or_equal_to( 10 ).test( 9 ) ) ;
	TEST_ASSERT( statfile::Constraint::less_than_or_equal_to( 10 ).test( "9" ) ) ;
	TEST_ASSERT( statfile::Constraint::less_than_or_equal_to( 10 ).test( 9.0 ) ) ;

	TEST_ASSERT( statfile::Constraint::less_than_or_equal_to( 10 ).test( 10 ) ) ;
	TEST_ASSERT( statfile::Constraint::less_than_or_equal_to( 10 ).test( "10" ) ) ;
	TEST_ASSERT( statfile::Constraint::less_than_or_equal_to( 10 ).test( 10.0 ) ) ;
	TEST_ASSERT( statfile::Constraint::less_than_or_equal_to( 10 ).test( "10.0" ) ) ;

	TEST_ASSERT( !statfile::Constraint::less_than_or_equal_to( 10 ).test( 11 ) ) ;
	TEST_ASSERT( !statfile::Constraint::less_than_or_equal_to( 10 ).test( 10.1 ) ) ;
	TEST_ASSERT( !statfile::Constraint::less_than_or_equal_to( 10 ).test( "11" ) ) ;
	TEST_ASSERT( !statfile::Constraint::less_than_or_equal_to( 10 ).test( "10.1" ) ) ;
}

AUTO_TEST_CASE( test_greater_than ) {
	TEST_ASSERT( !statfile::Constraint::greater_than( "abcd" ).test( "abcd" ) ) ;
	TEST_ASSERT( !statfile::Constraint::greater_than( "abcd" ).test( "aacd" ) ) ;
	TEST_ASSERT( !statfile::Constraint::greater_than( "abcd" ).test( "abbd" ) ) ;
	TEST_ASSERT( !statfile::Constraint::greater_than( "abcd" ).test( "abcc" ) ) ;
	TEST_ASSERT( statfile::Constraint::greater_than( "abcd" ).test( "bbcd" ) ) ;
	TEST_ASSERT( statfile::Constraint::greater_than( "abcd" ).test( "accd" ) ) ;

	TEST_ASSERT( !statfile::Constraint::greater_than( 10 ).test( 9 ) ) ;
	TEST_ASSERT( !statfile::Constraint::greater_than( 10 ).test( "9" ) ) ;
	TEST_ASSERT( !statfile::Constraint::greater_than( 10 ).test( 9.0 ) ) ;

	TEST_ASSERT( !statfile::Constraint::greater_than( 10 ).test( 10 ) ) ;
	TEST_ASSERT( !statfile::Constraint::greater_than( 10 ).test( "10" ) ) ;
	TEST_ASSERT( !statfile::Constraint::greater_than( 10 ).test( 10.0 ) ) ;
	TEST_ASSERT( !statfile::Constraint::greater_than( 10 ).test( "10.0" ) ) ;

	TEST_ASSERT( statfile::Constraint::greater_than( 10 ).test( 11 ) ) ;
	TEST_ASSERT( statfile::Constraint::greater_than( 10 ).test( 10.1 ) ) ;
	TEST_ASSERT( statfile::Constraint::greater_than( 10 ).test( "11" ) ) ;
	TEST_ASSERT( statfile::Constraint::greater_than( 10 ).test( "10.1" ) ) ;
}

AUTO_TEST_CASE( test_greater_than_or_equal_to ) {
	TEST_ASSERT( statfile::Constraint::greater_than_or_equal_to( "abcd" ).test( "abcd" ) ) ;
	TEST_ASSERT( !statfile::Constraint::greater_than_or_equal_to( "abcd" ).test( "aacd" ) ) ;
	TEST_ASSERT( !statfile::Constraint::greater_than_or_equal_to( "abcd" ).test( "abbd" ) ) ;
	TEST_ASSERT( !statfile::Constraint::greater_than_or_equal_to( "abcd" ).test( "abcc" ) ) ;
	TEST_ASSERT( statfile::Constraint::greater_than_or_equal_to( "abcd" ).test( "bbcd" ) ) ;
	TEST_ASSERT( statfile::Constraint::greater_than_or_equal_to( "abcd" ).test( "accd" ) ) ;

	TEST_ASSERT( !statfile::Constraint::greater_than_or_equal_to( 10 ).test( 9 ) ) ;
	TEST_ASSERT( !statfile::Constraint::greater_than_or_equal_to( 10 ).test( "9" ) ) ;
	TEST_ASSERT( !statfile::Constraint::greater_than_or_equal_to( 10 ).test( 9.0 ) ) ;

	TEST_ASSERT( statfile::Constraint::greater_than_or_equal_to( 10 ).test( 10 ) ) ;
	TEST_ASSERT( statfile::Constraint::greater_than_or_equal_to( 10 ).test( "10" ) ) ;
	TEST_ASSERT( statfile::Constraint::greater_than_or_equal_to( 10 ).test( 10.0 ) ) ;
	TEST_ASSERT( statfile::Constraint::greater_than_or_equal_to( 10 ).test( "10.0" ) ) ;

	TEST_ASSERT( statfile::Constraint::greater_than_or_equal_to( 10 ).test( 11 ) ) ;
	TEST_ASSERT( statfile::Constraint::greater_than_or_equal_to( 10 ).test( 10.1 ) ) ;
	TEST_ASSERT( statfile::Constraint::greater_than_or_equal_to( 10 ).test( "11" ) ) ;
	TEST_ASSERT( statfile::Constraint::greater_than_or_equal_to( 10 ).test( "10.1" ) ) ;
}

AUTO_TEST_CASE( test_between ) {
	TEST_ASSERT( !statfile::Constraint::between( "B", "D" ).test( "A" ) ) ;
	TEST_ASSERT( statfile::Constraint::between( "B", "D" ).test( "B" ) ) ;
	TEST_ASSERT( statfile::Constraint::between( "B", "D" ).test( "C" ) ) ;
	TEST_ASSERT( statfile::Constraint::between( "B", "D" ).test( "D" ) ) ;
	TEST_ASSERT( !statfile::Constraint::between( "B", "D" ).test( "E" ) ) ;

	TEST_ASSERT( !statfile::Constraint::between( 1, 10 ).test( 0 ) ) ;
	TEST_ASSERT( !statfile::Constraint::between( 1, 10 ).test( 0.0 ) ) ;
	TEST_ASSERT( statfile::Constraint::between( 1, 10 ).test( 1 ) ) ;
	TEST_ASSERT( statfile::Constraint::between( 1, 10 ).test( 5 ) ) ;
	TEST_ASSERT( statfile::Constraint::between( 1, 10 ).test( 10 ) ) ;
	TEST_ASSERT( statfile::Constraint::between( 1, 10 ).test( 10.0 ) ) ;
	TEST_ASSERT( !statfile::Constraint::between( 1, 10 ).test( 11 ) ) ;
	TEST_ASSERT( !statfile::Constraint::between( 1, 10 ).test( 10.1 ) ) ;
}

namespace {
	statfile::FilteringStatSource::UniquePtr
		build_filtering_source(
			std::string column,
			statfile::Constraint const& constraint
		) {
			statfile::FilteringStatSource::UniquePtr source(
				new statfile::FilteringStatSource(
					statfile::DelimitedStatSource::UniquePtr(
						new statfile::DelimitedStatSource(
							std::auto_ptr< std::istream >( new std::istringstream( globals::data ) ),
							" "
						)
					),
					column,
					constraint
				)
			) ;
			return source ;
		}
}

AUTO_TEST_CASE( test_filter ) 
{
	std::cerr << "Testing filtering stat source..." ;
	std::string value ;
	
	{
		statfile::FilteringStatSource::UniquePtr source = build_filtering_source(
			"Column2",
			statfile::Constraint::equal_to( "hello" )
		) ;
		(*source) >> value >> statfile::ignore_all() ;
		TEST_ASSERT( value == "0" ) ;
		TEST_ASSERT( !(*source) ) ;
	}

	{
		statfile::FilteringStatSource::UniquePtr source = build_filtering_source(
			"Column2",
			statfile::Constraint::equal_to( "hello2" )
		) ;

		TEST_ASSERT( !(*source) ) ;
	}

	{
		statfile::FilteringStatSource::UniquePtr source = build_filtering_source(
			"Column2",
			statfile::Constraint::greater_than( "A" )
		) ;

		(*source) >> value >> statfile::ignore_all() ;
		TEST_ASSERT( value == "0" ) ;
		(*source) >> value >> statfile::ignore_all() ;
		TEST_ASSERT( value == "2" ) ;
		(*source) >> value >> statfile::ignore_all() ;
		TEST_ASSERT( value == "3" ) ;
		TEST_ASSERT( !(*source) ) ;
	}

	{
		statfile::FilteringStatSource::UniquePtr source = build_filtering_source(
			"v8",
			statfile::Constraint::greater_than( 1000000 )
		) ;

		(*source) >> value >> statfile::ignore_all() ;
		TEST_ASSERT( value == "2" ) ;
		TEST_ASSERT( !(*source) ) ;
	}
}

BOOST_AUTO_TEST_SUITE_END()
