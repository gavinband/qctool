
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_WILDCARD_HPP
#define GENFILE_WILDCARD_HPP

#include <iostream>
#include <algorithm>
#include <utility>
#include <string>
#include <vector>
#include <iterator>
#include <iostream>
#include "snp_data_utils.hpp"

namespace genfile {
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
		bool operator==( FilenameMatch const& left, FilenameMatch const& right ) ;
		std::ostream& operator<<( std::ostream& oStream, FilenameMatch const& match ) ;
	
		enum WildcardMatchChoice { eALL_CHROMOSOMES, eAUTOSOMAL_CHROMOSOMES, eNON_SEX_CHROMOSOMES, eSEX_CHROMOSOMES } ;
	
		// Find all files matching the template.
		// If the template contains the wildcard char
		std::vector< FilenameMatch >
		find_files_by_chromosome(
			std::string path,
			WildcardMatchChoice choice = eNON_SEX_CHROMOSOMES,
			char wildcard_char = '#'
		) ;

		// Return a list of FilenameMatches consisting of filenames
		// made from the given output_filename by substituting the matches
		// from the given input filenames in place of the wildcard (if any)
		// in the output filename.  The size of the returned list is always
		// the same as the size of the input_filename_matches argument.
		std::vector< FilenameMatch >
		construct_corresponding_filenames(
			std::vector< FilenameMatch > const& filename_matches,
			std::string const& template_filename,
			char wildcard_char = '#'
		) ;
	}
}

#endif

