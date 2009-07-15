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
	std::pair< std::vector< std::string >, std::vector< std::string > > filenames = find_files_matching_path_with_wildcard( path_with_wildcards ) ;
	TEST_ASSERT( filenames.first.size() == filenames.second.size() ) ;
	if( expect )
		TEST_ASSERT( filenames.first.size() > 0 ) ;
	else
		TEST_ASSERT( filenames.first.size() == 0 ) ;

	for( std::size_t i = 0; i < filenames.first.size(); ++i ) {
		std::cout << "file: " << filenames.first[i] << ", matching part: " << filenames.second[i] << ".\n" ;
	}
}

AUTO_TEST_CASE( test_find_files_matching_path_with_wildcard ) {
	do_test( "waf-1*", true ) ;
	do_test( "w*", true ) ;
	do_test( "waf*.5.8", true ) ;
	do_test( "./waf-1*", true ) ;
	do_test( "./w*", true ) ;
	do_test( "./waf*.5.8", true ) ;
	do_test( "asdasdsad454545waf-1.5.8", false ) ;
	do_test( "./asdasdsad454545waf-1.5.8", false ) ;
}

#ifndef HAVE_BOOST_UNIT_TEST
	int main( int argc, char** argv ) {
		test_find_files_matching_path_with_wildcard() ;
	}
#endif
