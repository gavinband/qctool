#include "test_case.hpp"
#include "genfile/string_utils/string_utils.hpp"

using namespace genfile::string_utils ;

AUTO_TEST_CASE( test_split ) {
	std::cerr << "test_split()..." ;
	std::vector< std::string > expected ;
	TEST_ASSERT( expected == split( "", "," )) ;

	expected.push_back( "hello" ) ;
	TEST_ASSERT( expected == split( "hello", "," )) ;

	expected.push_back( "there" ) ;
	TEST_ASSERT( expected == split( "hello,there", "," )) ;

	expected.push_back( " you" ) ;
	TEST_ASSERT( expected == split( "hello,there, you", "," )) ;

	expected.push_back( " beautiful\t#" ) ;
	TEST_ASSERT( expected == split( "hello,there, you, beautiful\t#", "," )) ;

	expected.push_back( std::string( 1, 0 ) + "creature" ) ;
	TEST_ASSERT( expected == split( "hello,there, you, beautiful\t#," + std::string( 1, 0 ) + "creature", "," )) ;

	std::cerr << "ok.\n" ;
}

AUTO_TEST_MAIN {
	test_split() ;
}
