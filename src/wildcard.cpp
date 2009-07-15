#include <iostream>
#include <string>
#include <cassert>
#include "../config.hpp"
#if HAVE_BOOST_FILESYSTEM
	#include <boost/filesystem.hpp>
	namespace BFS = boost::filesystem ;
#endif
#include "parse_utils.hpp"
#include "wildcard.hpp"

namespace impl {
	bool check_if_filename_matches_filename_with_wildcard( std::string filename, std::string filename_before_wildcard, std::string filename_after_wildcard ) {
		std::string matching_segment ;
		if( string_has_prefix_and_suffix( filename, filename_before_wildcard, filename_after_wildcard, &matching_segment )) {
			return true ;
		}
		return false ;
	}
}

std::pair< std::vector< std::string >, std::vector< std::string > > find_files_matching_path_with_wildcard( std::string filename_with_wildcard, char wildcard_char ) {
	std::pair< std::vector< std::string >, std::vector< std::string > > result ;
#if HAVE_BOOST_FILESYSTEM
	BFS::path dir = BFS::path( filename_with_wildcard ).parent_path() ;
  	filename_with_wildcard = BFS::path( filename_with_wildcard ).filename() ;
	std::size_t wildcard_pos = filename_with_wildcard.find( wildcard_char ) ;
	if( dir.empty() ) {
		dir = "." ;
	}
	BFS::directory_iterator dir_i( dir ), end_i ;

	std::string filename_before_wildcard = filename_with_wildcard.substr( 0, wildcard_pos ) ;
	std::string filename_after_wildcard = "" ;
	bool have_wildcard = ( wildcard_pos != std::string::npos ) ;
	if( have_wildcard ) {
		filename_after_wildcard = filename_with_wildcard.substr( wildcard_pos + 1, filename_with_wildcard.size()) ;
	}

	for( ; dir_i != end_i; ++dir_i ) {
		if( BFS::is_regular_file( *dir_i )) {
			std::string filename = dir_i->filename();
			std::string matching_part ;
			if( have_wildcard && impl::check_if_filename_matches_filename_with_wildcard( filename, filename_before_wildcard, filename_after_wildcard )) {
				result.first.push_back(( dir / filename ).string()) ;
			}
			else if( filename == filename_before_wildcard ){
				result.first.push_back(( dir / filename ).string()) ;
			}
		}
	}
#else
	// Oh dear, no boost filesystem.  Just return the file as is.
	result.push_back( filename_with_wildcard ) ;
	result.push_back ( std::string( std::size_t(1), wildcard_char )) ;
#endif

	// Now construct the list of matching segments in the second part of our return value
	for( std::size_t i = 0; i < result.first.size(); ++i ) {
		std::string const& ith_filename = result.first[i] ;
		result.second.push_back( ith_filename.substr( filename_before_wildcard.size(), ith_filename.size() - filename_before_wildcard.size() - filename_after_wildcard.size() )) ;
	}

	return result ;
}
