#include <iostream>
#include "test_case.hpp"
#include "db/fill_SQL.hpp"

AUTO_TEST_CASE( test_fill_SQL ) {
	std::cerr << "test_fill_SQL()..." ;
	
	TEST_ASSERT( db::fill_SQL( "" ) == "" ) ;
	TEST_ASSERT( db::fill_SQL( "A string" ) == "A string" ) ;
	TEST_ASSERT( db::fill_SQL( "%s", "hello" ) == "\"hello\"" ) ;
	TEST_ASSERT( db::fill_SQL( "%s", 5 ) == "5" ) ;
	TEST_ASSERT( db::fill_SQL( "%s", 5.56 ) == "5.56" ) ;

	TEST_ASSERT( db::fill_SQL( "%s %s", "hello", "goodbye" ) == "\"hello\" \"goodbye\"" ) ;
	TEST_ASSERT( db::fill_SQL( "%s %s", "hello", 0 ) == "\"hello\" 0" ) ;
	TEST_ASSERT( db::fill_SQL( "%s %s", 4, 0 ) == "4 0" ) ;
	TEST_ASSERT( db::fill_SQL( "%s%s", 4, "doof" ) == "4\"doof\"" ) ;
	std::cerr << "success.\n" ;
}

AUTO_TEST_MAIN {
	test_fill_SQL() ;
}