#include <iostream>
#include <iomanip>
#include <sstream>
#include <cassert>
#include <string>
#include <vector>
#include "FileUtil.hpp"
#include "wildcard.hpp"

#if HAVE_BOOST_UNIT_TEST
	#define BOOST_AUTO_TEST_MAIN
	#include "boost/test/auto_unit_test.hpp"
	#define AUTO_TEST_CASE( param ) BOOST_AUTO_TEST_CASE(param)
	#define TEST_ASSERT( param ) BOOST_ASSERT( param )
#else
	#define AUTO_TEST_CASE( param ) void param()
	#define TEST_ASSERT( param ) assert( param )
#endif

void do_test( std::string path_with_wildcards, bool expect ) {
	std::cout << "Searching for files matching \"" << path_with_wildcards << "\"...\n" ;
	std::vector< wildcard::FilenameMatch > filenames = wildcard::find_files_matching_path_with_wildcard( path_with_wildcards ) ;
	if( expect )
		TEST_ASSERT( filenames.size() > 0 ) ;
	else
		TEST_ASSERT( filenames.size() == 0 ) ;

	for( std::size_t i = 0; i < filenames.size(); ++i ) {
		TEST_ASSERT( filenames[i].filename().size() >= filenames[i].match().size() ) ;
		TEST_ASSERT( filenames[i].filename().find( filenames[i].match() ) != std::string::npos ) ;
		std::cout << "file: " << filenames[i].filename() << ", matching part: " << filenames[i].match() << ".\n" ;
	}
}

void do_numerical_test( std::string path_with_wildcards, bool expect ) {
	std::cout << "Searching for files matching \"" << path_with_wildcards << "\"...\n" ;
	std::vector< wildcard::FilenameMatch > filenames = wildcard::find_matches_for_path_with_integer_wildcard( path_with_wildcards ) ;
	if( expect )
		TEST_ASSERT( filenames.size() > 0 ) ;
	else
		TEST_ASSERT( filenames.size() == 0 ) ;

	for( std::size_t i = 0; i < filenames.size(); ++i ) {
		TEST_ASSERT( filenames[i].filename().size() >= filenames[i].match().size() ) ;
		TEST_ASSERT( filenames[i].filename().find( filenames[i].match() ) != std::string::npos ) ;
		std::cout << "file: " << filenames[i].filename() << ", matching part: " << filenames[i].match() << ".\n" ;
	}
}

AUTO_TEST_CASE( test_find_files_matching_path_with_wildcard ) {
	// There should be a file in the top-level dir called ./waf-1.5.8
	// (since that's what we use to build).
	do_test( "waf-1*", true ) ;
	do_test( "w*", true ) ;
	do_test( "waf*.5.8", true ) ;
	do_test( "./waf-1*", true ) ;
	do_test( "./w*", true ) ;
	do_test( "./waf*.5.8", true ) ;
	do_test( "./waf-1.5.8", true ) ;
	do_test( "waf-1.5.8", true ) ;
	do_test( "asdasdsad454545waf-1.5.8", false ) ;
	do_test( "./asdasdsad454545waf-1.5.8", false ) ;

	do_numerical_test( "waf-#.5.8", true ) ;
	do_numerical_test( "./waf-#.5.8", true ) ;
	do_numerical_test( "waf-1.#.8", true ) ;
	do_numerical_test( "./waf-1.#.8", true ) ;
	do_numerical_test( "waf-1.5.#", true ) ;
	do_numerical_test( "./waf-1.5.#", true ) ;
	do_numerical_test( "waf-1.#", false ) ;
	do_numerical_test( "./waf-1.#", false ) ;
	do_numerical_test( "waf-#.8", false ) ;
	do_numerical_test( "./waf-#.8", false ) ;
	do_numerical_test( "waf-1.5.8", true ) ;
	do_numerical_test( "./waf-1.5.8", true ) ;
	do_numerical_test( "acacackakcak", false ) ;
	do_numerical_test( "./acacackakcak", false ) ;
	do_numerical_test( "acacackakcak#", false ) ;
	do_numerical_test( "./acacackakcak#", false ) ;
}

#ifndef HAVE_BOOST_UNIT_TEST
	int main( int argc, char** argv ) {
		test_find_files_matching_path_with_wildcard() ;
	}
#endif
