#include "test_case.hpp"
#include "genfile/string_utils/slice.hpp"

using namespace genfile::string_utils ;

AUTO_TEST_CASE( test_constructors ) {
	std::cerr << "test_constructors()..." ;

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

AUTO_TEST_CASE( test_split ) {
	std::cerr << "test_split()..." ;

	{
		std::string const test_string = "" ;
		std::vector< slice > elts = slice( test_string ).split( "," ) ;
		TEST_ASSERT( elts.size() == 0 ) ;
	}	

	{
		std::string const test_string = "" ;
		std::vector< slice > elts = slice( test_string ).split( ",\t" ) ;
		TEST_ASSERT( elts.size() == 0 ) ;
	}	

	{
		std::string const test_string = "," ;
		std::vector< slice > elts = slice( test_string ).split( "," ) ;
		TEST_ASSERT( elts.size() == 2 ) ;
		TEST_ASSERT( elts[0] == "" ) ;
		TEST_ASSERT( elts[1] == "" ) ;
	}	

	{
		std::string const test_string = "Hello,This, is me, and\tI, am A" ;
		std::vector< slice > elts = slice( test_string ).split( "," ) ;
		TEST_ASSERT( elts.size() == 5 ) ;
		TEST_ASSERT( elts[0] == "Hello" ) ;
		TEST_ASSERT( elts[1] == "This" ) ;
		TEST_ASSERT( elts[2] == " is me" ) ;
		TEST_ASSERT( elts[3] == " and\tI" ) ;
		TEST_ASSERT( elts[4] == " am A" ) ;
	}	

	{
		std::string const test_string = "Hello,This, is me, and\tI, am A" ;
		std::vector< slice > elts = slice( test_string ).split( ",\t" ) ;
		TEST_ASSERT( elts.size() == 6 ) ;
		TEST_ASSERT( elts[0] == "Hello" ) ;
		TEST_ASSERT( elts[1] == "This" ) ;
		TEST_ASSERT( elts[2] == " is me" ) ;
		TEST_ASSERT( elts[3] == " and" ) ;
		TEST_ASSERT( elts[4] == "I" ) ;
		TEST_ASSERT( elts[5] == " am A" ) ;
	}	

	{
		std::string const test_string = "Hello,This, is me, and\tI, am A" ;
		std::vector< slice > elts = slice( test_string ).split( ",\tl" ) ;
		TEST_ASSERT( elts.size() == 8 ) ;
		TEST_ASSERT( elts[0] == "He" ) ;
		TEST_ASSERT( elts[1] == "" ) ;
		TEST_ASSERT( elts[2] == "o" ) ;
		TEST_ASSERT( elts[3] == "This" ) ;
		TEST_ASSERT( elts[4] == " is me" ) ;
		TEST_ASSERT( elts[5] == " and" ) ;
		TEST_ASSERT( elts[6] == "I" ) ;
		TEST_ASSERT( elts[7] == " am A" ) ;
	}
	std::cerr << "ok.\n" ;
}

AUTO_TEST_MAIN {
	test_constructors() ;
	test_find() ;
	test_find_first_of() ;
	test_split() ;
}
