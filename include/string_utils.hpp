#ifndef __GTOOL__STRING_UTILS_HPP__
#define __GTOOL__STRING_UTILS_HPP__

#include <vector>
#include <string>
#include <sstream>

std::vector< std::string > split( std::string string_to_split, std::string splitter ) ;
std::string strip( std::string string_to_strip, std::string chars = " " ) ;
std::vector< std::string > split_discarding_empty_entries( std::string string_to_split, std::string splitter ) ;
std::vector< std::string > split_and_strip( std::string string_to_split, std::string splitter, std::string strip_chars = " " ) ;
std::vector< std::string > split_and_strip_discarding_empty_entries( std::string string_to_split, std::string splitter, std::string strip_chars = " " ) ;
std::string to_lower( std::string ) ;
std::string to_upper( std::string ) ;
template< typename T > std::string to_string( T const& t ) ;


//
// Implementation
//

template< typename T >
std::string to_string( T const& t ) {
	std::ostringstream os ;
	os << t ;
	return os.str() ;
}

#endif
