#include <iostream>
#include <string>
#include <map>
#include <cassert>
#include "../config.hpp"
#if HAVE_BOOST_FILESYSTEM
	#include <boost/filesystem.hpp>
	namespace BFS = boost::filesystem ;
#endif
#include "parse_utils.hpp"
#include "wildcard.hpp"

namespace wildcard {
	namespace impl {
	#if HAVE_BOOST_FILESYSTEM
		std::vector< FilenameMatch > find_files_with_prefix_and_suffix( BFS::path const& directory, std::string const& prefix, std::string const& suffix ) {
			std::vector< FilenameMatch > result ;
			BFS::directory_iterator dir_i( directory ), end_dir_i ;
			for( ; dir_i != end_dir_i; ++dir_i ) {
				std::string filename = dir_i->filename() ;
				std::string matching_part ;
				if( string_has_prefix_and_suffix( filename, prefix, suffix, &matching_part )) {
					result.push_back( FilenameMatch( (directory / filename).string(), matching_part )) ;
				}
			}
			return result ;
		}

		BFS::path get_dir_part( std::string const& path ) {
			BFS::path dir = BFS::path( path ).parent_path() ;
			if( dir.empty() ) {
				dir = "." ;
			}
			return dir ;
		}
	#endif

		bool has_wildcard( std::string const& a_string, char const wildcard ) {
			return a_string.find( wildcard ) != std::string::npos ;
		}
	}

	std::vector< FilenameMatch > find_files_matching_path_with_wildcard( std::string filename_with_wildcard, char wildcard_char ) {
		std::vector< FilenameMatch > result ;
	#if HAVE_BOOST_FILESYSTEM
		BFS::path dir = impl::get_dir_part( filename_with_wildcard ) ;BFS::path( filename_with_wildcard ).parent_path() ;
	  	filename_with_wildcard = BFS::path( filename_with_wildcard ).filename() ;
		if( !impl::has_wildcard( filename_with_wildcard, wildcard_char )) {
			if( BFS::exists( dir / filename_with_wildcard )) {
				result.push_back( FilenameMatch( (dir / filename_with_wildcard).string(), "" )) ;
			}
		}
		else {
			std::size_t wildcard_pos = filename_with_wildcard.find( wildcard_char ) ;
			std::string filename_before_wildcard = filename_with_wildcard.substr( 0, wildcard_pos ) ;
			std::string filename_after_wildcard = filename_with_wildcard.substr( wildcard_pos + 1, filename_with_wildcard.size() ) ;

			result = impl::find_files_with_prefix_and_suffix( dir, filename_before_wildcard, filename_after_wildcard ) ;
			std::sort( result.begin(), result.end() ) ;
		}
	#else
		// Oh dear, no boost support.  Return empty list.
	#endif
		return result ;
	}

	bool operator<( FilenameMatch const& left, FilenameMatch const& right ) {
		return left.match() < right.match() ;
	}

	std::vector< FilenameMatch > find_files_matching_path_with_integer_wildcard(
		std::string filename_with_wildcard,
		char wildcard_char,
		int match_lower_bound,
		int match_upper_bound
	) {
		std::vector< FilenameMatch > result ;
	#if HAVE_BOOST_FILESYSTEM
		BFS::path dir = impl::get_dir_part( filename_with_wildcard ) ;BFS::path( filename_with_wildcard ).parent_path() ;
	  	filename_with_wildcard = BFS::path( filename_with_wildcard ).filename() ;
		if( !impl::has_wildcard( filename_with_wildcard, wildcard_char )) {
			if( BFS::exists( dir / filename_with_wildcard )) {
				result.push_back( FilenameMatch( (dir / filename_with_wildcard).string(), "" )) ;
			}
		}
		else {
			std::size_t wildcard_pos = filename_with_wildcard.find( wildcard_char ) ;
			std::string filename_before_wildcard = filename_with_wildcard.substr( 0, wildcard_pos ) ;
			std::string filename_after_wildcard = filename_with_wildcard.substr( wildcard_pos + 1, filename_with_wildcard.size() ) ;

			std::vector< FilenameMatch > candidates = impl::find_files_with_prefix_and_suffix( dir, filename_before_wildcard, filename_after_wildcard ) ;

			// Make a list of the candidates for which the matching part is an integer in the given range.
			// Order by the number itself.  Note that the actual matching string may be prefixed with
			// zeroes, so ordering by the matching string is the wrong thing to do.
			std::map< int, FilenameMatch > matching_files ;
			for( std::size_t i = 0; i < candidates.size(); ++i ) {
				int matching_number = parse_integer_in_half_open_range( candidates[i].match(), match_lower_bound, match_upper_bound + 1 ) ;
				if( matching_number <= match_upper_bound ) {
					matching_files[ matching_number ] = candidates[i] ;
				}
			}
		
			// Copy into the result.
			for( std::map< int, FilenameMatch >::const_iterator i = matching_files.begin() ;
				i != matching_files.end() ;
				++i
			) {
				result.push_back( i->second ) ;
			}
		}
	#else
		// Oh dear, no boost support.  Return empty map.
	#endif


		return result ;
	}
}