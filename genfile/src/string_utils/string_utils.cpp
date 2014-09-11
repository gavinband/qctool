
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <limits>
#include <vector>
#include <string>
#include <sstream>
#include <exception>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include "genfile/Error.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/string_utils/strtod.hpp"

namespace genfile {
	namespace string_utils {

		std::string to_lower( std::string aString ) {
			for( std::string::iterator i = aString.begin(); i != aString.end(); ++i ) {
				if( *i >= 'A' && *i <= 'Z' ) {
					*i += 32 ;
				}
			}
			return aString ;
		}

		void to_lower( std::string* aString ) {
			for( std::string::iterator i = aString->begin(); i != aString->end(); ++i ) {
				if( *i >= 'A' && *i <= 'Z' ) {
					*i += 32 ;
				}
			}
		}

		std::string to_upper( std::string aString ) {
			for( std::string::iterator i = aString.begin(); i != aString.end(); ++i ) {
				if( *i >= 'a' && *i <= 'z' ) {
					*i -= 32 ;
				}
			}
			return aString ;
		}

		void to_upper( std::string* aString ) {
			for( std::string::iterator i = aString->begin(); i != aString->end(); ++i ) {
				if( *i >= 'a' && *i <= 'z' ) {
					*i -= 32 ;
				}
			}
		}

		bool case_insensitive_equality( std::string const& left, std::string const& right ) {
			return to_lower( left ) == to_lower( right ) ;
		}

		std::string strip( std::string string_to_strip, std::string chars ) {
			if( chars.empty() ) {
				return string_to_strip ;
			}

			std::size_t lpos = string_to_strip.find_first_not_of( chars ) ;
			if( lpos == std::string::npos ) {
				return "" ;
			}

			std::size_t rpos = string_to_strip.find_last_not_of( chars ) ;
			return string_to_strip.substr( lpos, rpos - lpos + 1 ) ;
		}

		// Specialisations of to_repr for speed purposes.
		// (std::istringstream appears painfully slow in this case).
		template<> double to_repr( std::string const& s ) {
			return strtod( s ) ;
		}

		template<> int to_repr( std::string const& s ) {
			if( s.empty() || std::isspace( s[0] )) {
				throw StringConversionError() ;
			}
			for( std::size_t i = 0; i < s.size(); ++i ) {
				if( s[i] == char( 0 ) ) {
					throw StringConversionError() ;
				}
			}
			char* endptr ;
			long int result = strtol( s.c_str(), &endptr, 10 ) ;
			if( endptr == s.c_str() || std::size_t( endptr - s.c_str() ) != s.size() || result > std::numeric_limits< int >::max() ) {
				throw StringConversionError() ;
			}
			return int( result ) ;
		}

		template<> long to_repr( std::string const& s ) {
			if( s.empty() || std::isspace( s[0] )) {
				throw StringConversionError() ;
			}
			for( std::size_t i = 0; i < s.size(); ++i ) {
				if( s[i] == char( 0 ) ) {
					throw StringConversionError() ;
				}
			}
			char* endptr ;
			long result = strtol( s.c_str(), &endptr, 10 ) ;
			if( endptr == s.c_str() || std::size_t( endptr - s.c_str() ) != s.size() || result > std::numeric_limits< int >::max() ) {
				throw StringConversionError() ;
			}
			return result ;
		}

		template<> long long to_repr( std::string const& s ) {
			if( s.empty() || std::isspace( s[0] )) {
				throw StringConversionError() ;
			}
			for( std::size_t i = 0; i < s.size(); ++i ) {
				if( s[i] == char( 0 ) ) {
					throw StringConversionError() ;
				}
			}
			char* endptr ;
			long long result = strtol( s.c_str(), &endptr, 10 ) ;
			if( endptr == s.c_str() || std::size_t( endptr - s.c_str() ) != s.size() || result > std::numeric_limits< int >::max() ) {
				throw StringConversionError() ;
			}
			return result ;
		}

		namespace impl {
			std::string strip( std::string const& string_to_strip, std::size_t start, std::size_t end, std::string chars ) {
				assert( end >= start ) ;
				if( chars.empty() || start == end ) {
					return string_to_strip.substr( start, end - start ) ;
				}

				std::size_t lpos = string_to_strip.find_first_not_of( chars, start ) ;
				if( lpos == std::string::npos || lpos >= end ) {
					return "" ;
				}

				std::size_t rpos = string_to_strip.find_last_not_of( chars, end - 1 ) ;
				assert( rpos >= start ) ;
				return string_to_strip.substr( lpos, rpos - lpos + 1 ) ;
			}

			std::vector< std::string > split(
				std::string string_to_split,
				std::string const& split_chars,
				std::string const& strip_chars,
				bool preserve_empty_entries
			) {
				std::vector< std::string > substrings ;
				std::size_t last_pos = 0, pos = 0 ;
				for( ; pos != string_to_split.size(); last_pos = pos + 1 ) {
					pos = string_to_split.find_first_of( split_chars, last_pos ) ;
					if( pos == std::string::npos ) {
						pos = string_to_split.size() ;
					}
					std::string substr = strip( string_to_split, last_pos, pos, strip_chars ) ;
					if( preserve_empty_entries || substr.size() > 0 ) {
						substrings.push_back( substr ) ;
					}
				}

				return substrings ; 
			}
		}
		
		std::vector< std::string > split( std::string const& string_to_split, std::string const& split_chars ) {
			return impl::split( string_to_split, split_chars, "", true ) ;
		}
		
		std::vector< std::string > split_and_strip( std::string const& string_to_split, std::string const& split_chars, std::string const& strip_chars ) {
			return impl::split( string_to_split, split_chars, strip_chars, true ) ;
		}

		std::vector< std::string > split_and_strip_discarding_empty_entries( std::string const& string_to_split, std::string const& split_chars, std::string const& strip_chars ) {
			return impl::split( string_to_split, split_chars, strip_chars, false ) ;
		}
		
		std::string strip_all_whitespace( std::string input ) {
			std::string::iterator i ;
			i = std::remove( input.begin(), input.end(), ' ' ) ;
			i = std::remove( input.begin(), i, '\t' ) ;
			i = std::remove( input.begin(), i, '\n' ) ;
			return input.substr( 0, ( i - input.begin() )) ;
		}
		
		std::vector< std::string > split_respecting_delimited_regions(
			std::string const& line,
			std::string const& split_chars,
			std::string const& delimiters
		) {
			assert( delimiters.size() == 2 ) ;
			assert( split_chars.find( delimiters[0] ) == std::string::npos ) ;
			std::vector< std::string > result ;
			int in_quote = 0 ;
			std::size_t last_i = 0 ;
			for( std::size_t i = 0; i < line.size(); ++i ) {
				if( !in_quote && split_chars.find( line[i] ) != std::string::npos ) {
					result.push_back( line.substr( last_i, i - last_i )) ;
					last_i = i + 1 ;
				}
				else if( line[i] == delimiters[0] || line[i] == delimiters[1] ) {
					in_quote = ( in_quote + 1 ) % 2 ;
				}
			}
			if( in_quote ) {
				throw genfile::BadArgumentError( "genfile::string_utils::split_respecting_delimited_regions()", "line=\"" + line + "\"" ) ;
			}
			result.push_back( line.substr( last_i, line.size() - last_i )) ;

			return result ;
		}
		
		std::string join( std::vector< std::string > const& strings, std::string const& joiner ) {
			std::string result ;
			for( std::size_t i = 0; i < strings.size(); ++i ) {
					if( i > 0 ) {
							result += joiner ;
					}
					result += strings[i] ;
			}
			return result ;
		}
		
		namespace impl {
				bool is_alphabetic( char c) {
						return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ;
				}

				bool is_blank( char c) {
						return c == ' ' || c == '\n' ;
				}
		}
		
		std::string replace_all( std::string in, std::string const& pattern, std::string const& replacement ) {
			for( std::size_t pos = in.find( pattern ) ; pos != std::string::npos; pos = in.find( pattern ) ) {
				in.replace( pos, pattern.size(), replacement ) ;
			}
			return in ;
		}
		
		std::string wrap( std::string const& string_to_wrap, unsigned int wrap_column, unsigned int starting_column, std::size_t indent_amount ) { 
			assert( wrap_column > starting_column ) ; 
			assert( wrap_column > indent_amount ) ; 

			if( string_to_wrap.size() < (wrap_column - starting_column) ) { 
					return string_to_wrap ;
			}	

			std::string result ;
			unsigned int current_column = starting_column ;
			std::string indent( indent_amount, ' ' ) ; 

			std::string::const_iterator
					this_char = string_to_wrap.begin(),
					next_char = string_to_wrap.begin() ;
			for( ++next_char; this_char < string_to_wrap.end(); ++this_char, ++next_char ) { 
					if( current_column == 0 ) { 
							result.append( indent ) ; 
							current_column = indent_amount ;
							// skip whitespace at beginning of line
							for( ; this_char < string_to_wrap.end() && impl::is_blank( *this_char ); ++this_char, ++next_char ) ; 
					}	

					if( this_char < string_to_wrap.end() ) { 
							result.push_back( *this_char ) ; 
							if( *this_char == '\n' ) { 
									current_column = 0 ; 
							}	
							else if( ++current_column < wrap_column || next_char == string_to_wrap.end() ) { 
							}	
							else {
									if( impl::is_alphabetic( *this_char ) && impl::is_alphabetic( *next_char )) {
											// we are in the middle of a word.
											result.push_back( '-' ) ;
											result.push_back( '\n' ) ;
											current_column = 0 ;
									}
									else if( impl::is_blank( *this_char ) || impl::is_blank( *next_char )) {
											result.push_back( '\n' ) ;
											current_column = 0 ;
									}
									else {
											// do nothing in case this is something which can't be split.
									}
							}
					}
			}
			return result ;
		}
	}
	
}
