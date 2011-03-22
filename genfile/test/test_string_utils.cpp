#include "test_case.hpp"
#include "genfile/string_utils/string_utils.hpp"
#include "genfile/string_utils/slice.hpp"

using namespace genfile::string_utils ;

AUTO_TEST_CASE( test_split ) {
	std::cerr << "test_split()..." ;

	{
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
	}
	std::cerr << "ok.\n" ;
}

AUTO_TEST_CASE( test_slice_split ) {
	std::cerr << "test_slice_split()..." ;

	{
		std::string data = "" ;
		std::vector< slice > expected ;
		expected.push_back( data ) ;
		TEST_ASSERT( expected == slice( data ).split( "," )) ;
		TEST_ASSERT( expected == slice( data ).split( ":" )) ;
	}
	{
		std::string data = "hello" ;
		std::vector< slice > expected ;
		expected.push_back( data ) ;
		TEST_ASSERT( expected == slice( data ).split( "," )) ;
	}
	{
		std::string data = "hello,there" ;
		std::vector< slice > expected ;
		expected.push_back( slice( data, 0, 5 )) ;
		expected.push_back( slice( data, 6, 11 )) ;
		TEST_ASSERT( expected == slice( data ).split( "," )) ;
	}
	{
		std::string data = "hello,there, you" ;
		std::vector< slice > expected ;
		expected.push_back( slice( data, 0, 5 ) ) ;
		expected.push_back( slice( data, 6, 11 ) ) ;
		expected.push_back( slice( data, 12, 16 ) ) ;
		TEST_ASSERT( expected == slice( data ).split( "," )) ;
	}
	{
		std::string data = "hello,there, you, beautiful\t#" ;
		std::vector< slice > expected ;
		expected.push_back( slice( data, 0, 5 ) ) ;
		expected.push_back( slice( data, 6, 11 ) ) ;
		expected.push_back( slice( data, 12, 16 ) ) ;
		expected.push_back( slice( data, 17, 29 ) ) ;
		TEST_ASSERT( expected == slice( data ).split( "," )) ;
	}
	{
		std::string data = "hello,there, you, beautiful\t#" + std::string( 1, 0 ) + "creature" ;
		std::vector< slice > expected ;
		expected.push_back( slice( data, 0, 5 ) ) ;
		expected.push_back( slice( data, 6, 11 ) ) ;
		expected.push_back( slice( data, 12, 16 ) ) ;
		expected.push_back( slice( data, 17, data.size() ) ) ;
		TEST_ASSERT( expected == slice( data ).split( "," )) ;
	}

	{
		std::string data = "hello,there, you, beautiful\t#" + std::string( 1, 0 ) + "creature" ;
		std::vector< slice > expected ;
		expected.push_back( slice( data, 0, 5 ) ) ;
		expected.push_back( slice( data, 6, 11 ) ) ;
		expected.push_back( slice( data, 12, 16 ) ) ;
		expected.push_back( slice( data, 17, 29 ) ) ;
		expected.push_back( slice( data, 30, data.size() ) ) ;
		TEST_ASSERT( expected == slice( data ).split( std::string( "," ) + std::string( 1, 0 ) )) ;
	}

	{
		std::string data = "hello,there, you, beautiful\t#" + std::string( 1, 0 ) + "creature" ;
		std::vector< slice > expected ;
		expected.push_back( slice( data, 0, 5 ) ) ;
		expected.push_back( slice( data, 6, 11 ) ) ;
		expected.push_back( slice( data, 12, 16 ) ) ;
		expected.push_back( slice( data, 17, 27 ) ) ;
		expected.push_back( slice( data, 28, 29 ) ) ;
		expected.push_back( slice( data, 30, data.size() ) ) ;
		TEST_ASSERT( expected == slice( data ).split( std::string( ",\t" ) + std::string( 1, 0 ) )) ;
	}

	std::cerr << "ok.\n" ;
}

AUTO_TEST_MAIN {
	test_split() ;
	test_slice_split() ;
}