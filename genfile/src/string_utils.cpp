#include <vector>
#include <string>
#include <sstream>
#include <exception>
#include <algorithm>
#include "genfile/string_utils.hpp"

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

		namespace impl {
			std::vector< std::string > split( std::string string_to_split, std::string const& split_chars, std::string const& strip_chars, bool preserve_empty_entries ) {
				std::vector< std::string > substrings ;
				std::size_t pos = 0 ;
				do {
					pos = string_to_split.find_first_of( split_chars ) ;
					std::string substr = strip( string_to_split.substr(0, pos), strip_chars ) ;
					if( preserve_empty_entries || substr.size() > 0 ) {
						substrings.push_back( substr ) ;
					}
					if( pos != std::string::npos ) {
						string_to_split = string_to_split.substr( pos + 1, string_to_split.size() ) ;
					}
				}
				while( pos != std::string::npos ) ;

				return substrings ;	
			}
		}
		
		std::vector< std::string > split_and_strip( std::string string_to_split, std::string const& split_chars, std::string const& strip_chars ) {
			return impl::split( string_to_split, split_chars, strip_chars, true ) ;
		}

		std::vector< std::string > split_discarding_empty_entries( std::string string_to_split, std::string const& split_chars, std::string const& strip_chars ) {
			return impl::split( string_to_split, split_chars, strip_chars, false ) ;
		}
		
		std::string strip_all_whitespace( std::string input ) {
			std::string::iterator i ;
			i = std::remove( input.begin(), input.end(), ' ' ) ;
			i = std::remove( input.begin(), i, '\t' ) ;
			i = std::remove( input.begin(), i, '\n' ) ;
			return input.substr( 0, ( i - input.begin() )) ;
		}
	}
}
