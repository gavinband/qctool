#include <cassert>
#include <string>
#include "parse_utils.hpp"

std::size_t find_matching_bracket( std::string const& aString, std::size_t pos, char lbracket, char rbracket ) {
	assert( aString[pos] == lbracket ) ;
	int count = 1 ;
	std::size_t i = 1 ;
	for( ; i < aString.size() && count > 0; ++i ) {
		if( aString[i] == lbracket ) {
			++count ;
		}
		else if( aString[i] == rbracket ) {
			--count ;
		}
	}
	
	if( count == 0 ) {
		return i ;
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

