#ifndef __GTOOL__PARSE_UTILS_HPP__
#define __GTOOL__PARSE_UTILS_HPP__

#include <string>

std::size_t find_matching_bracket( std::string const& aString, std::size_t pos, char lbracket, char rbracket ) ;
std::string strip_brackets( std::string const& string_to_strip, char lbracket, char rbracket ) ;
bool string_has_prefix_and_suffix( std::string const& string_to_check, std::string const& prefix, std::string const& suffix, std::string* region_outside_prefix_and_suffix = 0 ) ;

// If the given string represents an integer in the interval [ lower_bound, upper_bound ),
// return that integer.  Otherwise return upper_bound.
int parse_integer_in_half_open_range( std::string const& a_string, int lower_bound, int upper_bound ) ;

#endif
