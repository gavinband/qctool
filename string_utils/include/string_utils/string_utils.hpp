#ifndef STRING_UTILS_HPP__
#define STRING_UTILS_HPP__

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <exception>

namespace string_utils {
	std::vector< std::string > split( std::string string_to_split, std::string splitter ) ;
	std::string join( std::vector< std::string > const& strings, std::string joiner = ", " ) ;
	std::string strip( std::string string_to_strip, std::string chars = " " ) ;
	std::vector< std::string > split_and_strip_discarding_empty_entries( std::string string_to_split, std::string splitter, std::string strip_chars = " " ) ;
	std::vector< std::string > split_and_strip( std::string string_to_split, std::string splitter, std::string strip_chars = " " ) ;
	std::string to_lower( std::string ) ;
	std::string to_upper( std::string ) ;
	std::string wrap( std::string const& string_to_wrap, unsigned int wrap_column, unsigned int starting_column = 0u, std::size_t indent_amount = 0u ) ;

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

	struct StringConversionError: public std::exception { char const* what() const throw() { return "StringConversionError" ; } } ;

	template< typename T >
	T to_repr( std::string const& s ) {
		std::istringstream is( s ) ;
		T t ;
		is >> t;
		if( !is ) {
			throw StringConversionError() ;
		}
		is.peek() ;
		if( is ) {
			throw StringConversionError() ;
		}
		return t ;
	}
}

#endif
