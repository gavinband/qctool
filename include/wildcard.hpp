#ifndef GTOOL_WILDCARD_HPP__
#define GTOOL_WILDCARD_HPP__

#include <iostream>
#include <algorithm>
#include <utility>
#include <string>
#include <vector>
#include <iterator>
#include <iostream>
#include <GToolException.hpp>

namespace wildcard {
	struct FilenameMatch {
		FilenameMatch( std::string const& filename = "", std::string const& match = "" ) : m_filename( filename ), m_match( match ) {}
		FilenameMatch( FilenameMatch const& other ): m_filename (other.m_filename), m_match( other.m_match ) {}
		std::string const& filename() const { return m_filename ; }
		std::string const& match() const { return m_match ; }
		private:
			std::string m_filename ;
			std::string m_match ;
	} ;

	bool operator<( FilenameMatch const& left, FilenameMatch const& right ) ;

	std::vector< FilenameMatch >
	find_files_matching_path_with_wildcard(
		std::string filename_with_wildcards,
		char wildcard_char = '*'
	) ;

	std::vector< FilenameMatch >
	find_files_matching_path_with_integer_wildcard(
		std::string filename_with_wildcard,
		char wildcard_char = '#',
		int match_lower_bound = 1,
		int match_upper_bound = 100
	) ;
}

#endif

