#ifndef GENFILE_STRING_UTILS_HPP
#define GENFILE_STRING_UTILS_HPP

#include <vector>
#include <string>
#include <sstream>
#include <exception>

namespace genfile {
	namespace string_utils {

		std::string to_lower( std::string aString ) ;
		bool case_insensitive_equality( std::string const& left, std::string const& right ) ;
		std::string strip( std::string string_to_strip, std::string chars ) ;
		std::vector< std::string > split( std::string const& string_to_split, std::string const& split_chars ) ;
		std::vector< std::string > split_and_strip(
			std::string const& string_to_split,
			std::string const& split_chars = " \t\n",
			std::string const& strip_chars = " \t\n"
		) ;
		std::vector< std::string > split_and_strip_discarding_empty_entries( std::string const& string_to_split, std::string const& split_chars = " \t\n", std::string const& strip_chars = " \t\n" ) ;
		std::string strip_all_whitespace( std::string input ) ;
		
		std::vector< std::string > split_respecting_delimited_regions(
			std::string const& line,
			std::string const& split_chars,
			std::string const& delimiters
		) ;
		
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
		
		// We specialise to_repr<> for some types for performance reasons.
		template<> double to_repr< double >( std::string const& s ) ;
		template<> int to_repr< int >( std::string const& s ) ;
	}
}

#endif
