
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "test_case.hpp"
#include "genfile/string_utils/slice.hpp"

using namespace genfile::string_utils ;

BOOST_AUTO_TEST_SUITE( test_slice ) ;

AUTO_TEST_CASE( test_constructors ) {
	std::cerr << "test_constructors()..." ;

	{
		char const* s = "Hello" ;
		TEST_ASSERT( std::string( slice( s )) == "Hello" ) ;
		TEST_ASSERT( slice( s ).size() == 5 ) ;
		TEST_ASSERT( slice( s )[0] == 'H' ) ;
		TEST_ASSERT( slice( s )[2] == 'l' ) ;
		TEST_ASSERT( slice( s )[4] == 'o' ) ;
		TEST_ASSERT( std::string( slice( s, 1, 4 ) ) == "ell" ) ;
		TEST_ASSERT( slice( s, 1, 4 ).size() == 3 ) ;
		TEST_ASSERT( slice( s, 1, 4 )[0] == 'e' ) ;
		TEST_ASSERT( slice( s, 1, 4 )[2] == 'l' ) ;
		TEST_ASSERT( std::string( slice( s, 3, 5 ) ) == "lo" ) ;
		TEST_ASSERT( slice( s, 3, 5 ).size() == 2 ) ;
		TEST_ASSERT( slice( s, 3, 5 )[0] == 'l' ) ;
		TEST_ASSERT( slice( s, 3, 5 )[1] == 'o' ) ;
		TEST_ASSERT( std::string( slice( s, 0, 2 ) ) == "He" ) ;
		TEST_ASSERT( slice( s, 0, 2 ).size() == 2 ) ;
		TEST_ASSERT( std::string( slice( s, 0, 0 ) ) == "" ) ;
		TEST_ASSERT( slice( s, 0, 0 ).size() == 0 ) ;
		TEST_ASSERT( std::string( slice( s, 5, 5 ) ) == "" ) ;
		TEST_ASSERT( slice( s, 5, 5 ).size() == 0 ) ;
	}

	{
		char const* s = "Hello" ;
		TEST_ASSERT( std::string( slice( s, s+5 )) == "Hello" ) ;
		TEST_ASSERT( slice( s, s+5 ).size() == 5 ) ;
		TEST_ASSERT( slice( s, s+5 )[0] == 'H' ) ;
		TEST_ASSERT( slice( s, s+5 )[2] == 'l' ) ;
		TEST_ASSERT( slice( s, s+5 )[4] == 'o' ) ;
		TEST_ASSERT( std::string( slice( s, s+5, 1, 4 ) ) == "ell" ) ;
		TEST_ASSERT( slice( s, s+5, 1, 4 ).size() == 3 ) ;
		TEST_ASSERT( slice( s, s+5, 1, 4 )[0] == 'e' ) ;
		TEST_ASSERT( slice( s, s+5, 1, 4 )[2] == 'l' ) ;
		TEST_ASSERT( std::string( slice( s, s+5, 3, 5 ) ) == "lo" ) ;
		TEST_ASSERT( slice( s, s+5, 3, 5 ).size() == 2 ) ;
		TEST_ASSERT( slice( s, s+5, 3, 5 )[0] == 'l' ) ;
		TEST_ASSERT( slice( s, s+5, 3, 5 )[1] == 'o' ) ;
		TEST_ASSERT( std::string( slice( s, s+5, 0, 2 ) ) == "He" ) ;
		TEST_ASSERT( slice( s, s+5, 0, 2 ).size() == 2 ) ;
		TEST_ASSERT( std::string( slice( s, s+5, 0, 0 ) ) == "" ) ;
		TEST_ASSERT( slice( s, s+5, 0, 0 ).size() == 0 ) ;
		TEST_ASSERT( std::string( slice( s, s+5, 5, 5 ) ) == "" ) ;
		TEST_ASSERT( slice( s, s+5, 5, 5 ).size() == 0 ) ;
	}

	{
		std::string const s = "Hello" ;
		TEST_ASSERT( std::string( slice( s )) == "Hello" ) ;
		TEST_ASSERT( slice( s ).size() == 5 ) ;
		TEST_ASSERT( slice( s )[0] == 'H' ) ;
		TEST_ASSERT( slice( s )[2] == 'l' ) ;
		TEST_ASSERT( slice( s )[4] == 'o' ) ;
		TEST_ASSERT( std::string( slice( s, 1, 4 ) ) == "ell" ) ;
		TEST_ASSERT( slice( s, 1, 4 ).size() == 3 ) ;
		TEST_ASSERT( slice( s, 1, 4 )[0] == 'e' ) ;
		TEST_ASSERT( slice( s, 1, 4 )[2] == 'l' ) ;
		TEST_ASSERT( std::string( slice( s, 3, 5 ) ) == "lo" ) ;
		TEST_ASSERT( slice( s, 3, 5 ).size() == 2 ) ;
		TEST_ASSERT( slice( s, 3, 5 )[0] == 'l' ) ;
		TEST_ASSERT( slice( s, 3, 5 )[1] == 'o' ) ;
		TEST_ASSERT( std::string( slice( s, 0, 2 ) ) == "He" ) ;
		TEST_ASSERT( slice( s, 0, 2 ).size() == 2 ) ;
		TEST_ASSERT( std::string( slice( s, 0, 0 ) ) == "" ) ;
		TEST_ASSERT( slice( s, 0, 0 ).size() == 0 ) ;
		TEST_ASSERT( std::string( slice( s, 5, 5 ) ) == "" ) ;
		TEST_ASSERT( slice( s, 5, 5 ).size() == 0 ) ;
	}

	{
		std::string const line = "Hello" ;
		slice const s( line ) ;
		TEST_ASSERT( std::string( s ) == "Hello" ) ;
		TEST_ASSERT( std::string( slice( s )) == "Hello" ) ;
		TEST_ASSERT( slice( s ).size() == 5 ) ;
		TEST_ASSERT( slice( s )[0] == 'H' ) ;
		TEST_ASSERT( slice( s )[2] == 'l' ) ;
		TEST_ASSERT( slice( s )[4] == 'o' ) ;
		TEST_ASSERT( std::string( slice( s, 1, 4 ) ) == "ell" ) ;
		TEST_ASSERT( std::string( s ) == "Hello" ) ;
		TEST_ASSERT( slice( s, 1, 4 ).size() == 3 ) ;
		TEST_ASSERT( slice( s, 1, 4 )[0] == 'e' ) ;
		TEST_ASSERT( slice( s, 1, 4 )[2] == 'l' ) ;
		TEST_ASSERT( std::string( slice( s, 3, 5 ) ) == "lo" ) ;
		TEST_ASSERT( std::string( s ) == "Hello" ) ;
		TEST_ASSERT( slice( s, 3, 5 ).size() == 2 ) ;
		TEST_ASSERT( slice( s, 3, 5 )[0] == 'l' ) ;
		TEST_ASSERT( slice( s, 3, 5 )[1] == 'o' ) ;
		TEST_ASSERT( std::string( slice( s, 0, 2 ) ) == "He" ) ;
		TEST_ASSERT( std::string( s ) == "Hello" ) ;
		TEST_ASSERT( slice( s, 0, 2 ).size() == 2 ) ;
		TEST_ASSERT( std::string( slice( s, 0, 0 ) ) == "" ) ;
		TEST_ASSERT( slice( s, 0, 0 ).size() == 0 ) ;
		TEST_ASSERT( std::string( slice( s, 5, 5 ) ) == "" ) ;
		TEST_ASSERT( slice( s, 5, 5 ).size() == 0 ) ;
	}
	std::cerr << "ok.\n" ;
}

AUTO_TEST_CASE( test_find ) {
	std::cerr << "test_find()..." ;

	{
		std::string const test_string = "Hello, I am a string.$$**88  \t" ;
		slice s( test_string ) ;
		for( std::size_t pos = 0; pos < test_string.size(); ++pos ) {
			for( int i = 0; i < 256; ++i ) {
				char c = char( i ) ;
				TEST_ASSERT( test_string.find( c ) == s.find( c )) ;
				TEST_ASSERT( test_string.find( c, pos ) == s.find( c, pos )) ;
			}
		}
	}
	{
		std::string const test_string = "Hello, I am a string.$$**88  \t" ;
		std::string const big_string = "PREFIXHello, I am a string.$$**88  \tSUFFIX" ;
		slice s( big_string, 6, 36 ) ;
		for( std::size_t pos = 0; pos < test_string.size(); ++pos ) {
			for( int i = 0; i < 256; ++i ) {
				char c = char( i ) ;
				TEST_ASSERT( s.get_start() == 6 ) ;
				TEST_ASSERT( s.get_end() == 36 ) ;
				TEST_ASSERT( test_string.find( c ) == s.find( c ) ) ;
				TEST_ASSERT( test_string.find( c, pos ) == s.find( c, pos ) ) ;
			}
		}
	}
	std::cerr << "ok.\n" ;
}

AUTO_TEST_CASE( test_find_first_of ) {
	std::cerr << "test_find_first_of()..." ;

	{
		std::string test_string = "Hello, I am a string.$$**88  \t" ;
		slice s( test_string ) ;
		std::string chars( 3, ' ' ) ;

		for( std::size_t pos = 0; pos < test_string.size(); ++pos ) {
			for( int i = 3; i < 256; ++i ) {
				char c = char( i ) ;
				chars[0] = c - 2 ;
				chars[1] = c - 1 ;
				chars[2] = c ;
				TEST_ASSERT( test_string.find_first_of( chars ) == s.find_first_of( chars )) ;
				TEST_ASSERT( test_string.find_first_of( chars, pos ) == s.find_first_of( chars, pos )) ;
			}
		}
	}

	{
		std::string const test_string = "Hello, I am a string.$$**88  \t" ;
		std::string const big_string = "PREFIXHello, I am a string.$$**88  \tSUFFIX" ;
		slice s( big_string, 6, 36 ) ;
		std::string chars( 3, ' ' ) ;
	
		for( std::size_t pos = 0; pos < test_string.size(); ++pos ) {
			for( int i = 3; i < 256; ++i ) {
				char c = char( i ) ;
				chars[0] = c - 2 ;
				chars[1] = c - 1 ;
				chars[2] = c ;
				TEST_ASSERT( test_string.find_first_of( chars ) == s.find_first_of( chars ) ) ;
				TEST_ASSERT( test_string.find_first_of( chars, pos ) == s.find_first_of( chars, pos ) ) ;
			}
		}
	}
	std::cerr << "ok.\n" ;
}

AUTO_TEST_CASE( test_split_slice ) {
	std::cerr << "test_split_slice()..." ;

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

BOOST_AUTO_TEST_SUITE_END() ;

