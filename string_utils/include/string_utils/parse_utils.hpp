#ifndef PARSE_UTILS_HPP
#define PARSE_UTILS_HPP

#include <utility>
#include <string>
#include "string_utils/string_utils.hpp"

namespace string_utils {
	// Given a string with left bracket character at the given pos,
	// Find the matching right bracket character to its right.
	// Nested brackets are taken into account.
	std::size_t find_matching_bracket( std::string const& aString, std::size_t pos, char lbracket, char rbracket ) ;
	std::string strip_brackets( std::string const& string_to_strip, char lbracket, char rbracket ) ;
	bool string_has_prefix_and_suffix( std::string const& string_to_check, std::string const& prefix, std::string const& suffix, std::string* region_outside_prefix_and_suffix = 0 ) ;

	// If the given string represents an integer in the interval [ lower_bound, upper_bound ),
	// return that integer.  Otherwise return upper_bound.
	int parse_integer_in_half_open_range( std::string const& a_string, int lower_bound, int upper_bound ) ;

	template< typename T >
	std::pair< T, T > parse_range( std::string const& range_spec ) {
		std::pair< T, T > result ;
		std::vector< std::string > bits = split( range_spec, "-" ) ;
		if( bits.size() == 1 ) {
			bits.push_back( bits[0] ) ;
		}
		assert( bits.size() == 2 ) ;
		result.first = to_repr< T >( bits[0] ) ;
		result.second = to_repr< T >( bits[1] ) ;
		return result ;
	}
}

#endif
