#include <cassert>
#include <string>
#include <sstream>
#include "string_utils/parse_utils.hpp"

namespace string_utils {
	std::size_t find_matching_bracket( std::string const& aString, std::size_t pos, char lbracket, char rbracket ) {
		assert( aString[pos] == lbracket ) ;
		int count = 1 ;
		std::size_t i = pos ;
		for( ++i; i < aString.size() && count > 0; ++i ) {
			if( aString[i] == lbracket ) {
				++count ;
			}
			else if( aString[i] == rbracket ) {
				--count ;
			}
		}
	
		if( count == 0 ) {
			return i-1 ;
		}
		else {
			return std::string::npos ;
		}
	}

	// Strip any matching brackets from the ends of the given string.
	std::string strip_brackets( std::string const& string_to_strip, char lbracket, char rbracket ) {
		if( string_to_strip.size() >= 2 && string_to_strip[0] == lbracket) {
			std::size_t rbracket_pos = find_matching_bracket( string_to_strip, 0, lbracket, rbracket ) ;
			if( rbracket_pos == string_to_strip.size() - 1 ) {
				return strip_brackets( string_to_strip.substr( 1, string_to_strip.size() - 2 ), lbracket, rbracket ) ;
			}
		}

		return string_to_strip ;
	}

	bool string_has_prefix_and_suffix( std::string const& string_to_check, std::string const& prefix, std::string const& suffix, std::string* region_outside_prefix_and_suffix ) {
		if( string_to_check.size() < ( prefix.size() + suffix.size() )) {
			return false ;
		}
		if( prefix != string_to_check.substr( 0, prefix.size()) ) {
			return false ;
		}
		if( suffix != string_to_check.substr( string_to_check.size() - suffix.size(), suffix.size() )) {
			return false ;
		}
		if( region_outside_prefix_and_suffix != 0) {
			std::size_t match_size = string_to_check.size() - prefix.size() - suffix.size() ;	
			region_outside_prefix_and_suffix->assign( string_to_check.substr( prefix.size(), match_size )) ;
		}
		return true ;
	}

	int parse_integer_in_half_open_range( std::string const& a_string, int lower_bound, int upper_bound )
	{
		assert( lower_bound <= upper_bound ) ;
		std::istringstream inStream( a_string ) ;
		int i ;
		inStream >> i ;
		inStream.peek() ;
		if( !inStream.eof()) {
			return upper_bound ;
		}
		if( i < lower_bound || i >= upper_bound ) {
			return upper_bound ;
		}
		return i ;
	}
}
