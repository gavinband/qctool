#include <cassert>
#include <vector>
#include <string>
#include "string_utils/string_utils.hpp"

namespace string_utils {
	namespace impl {

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

		enum EmptyEntryTreatmentSelector {ePreserveEmptyEntries, eDiscardEmptyEntries}  ;

		std::vector< std::string > split_and_strip( std::string string_to_split, std::string splitter, EmptyEntryTreatmentSelector empty_entry_treatment_selector, std::string strip_chars ) {
				std::vector< std::string > substrings ;
				std::size_t pos = 0 ;
				do {
					pos = string_to_split.find( splitter ) ;
					std::string substr = strip( string_to_split.substr(0, pos), strip_chars ) ;
					if( empty_entry_treatment_selector == ePreserveEmptyEntries || substr.size() > 0 ) {
						substrings.push_back( substr ) ;
					}
					if( pos != std::string::npos ) {
						string_to_split = string_to_split.substr( pos + splitter.size(), string_to_split.size() ) ;
					}
				}
				while( pos != std::string::npos ) ;

				return substrings ;
		}

	}

	std::vector< std::string > split( std::string string_to_split, std::string splitter ) {
		return impl::split_and_strip( string_to_split, splitter, impl::ePreserveEmptyEntries, "" ) ;
	}

	std::vector< std::string > split_and_strip_discarding_empty_entries( std::string string_to_split, std::string splitter, std::string strip_chars ) {
		return impl::split_and_strip( string_to_split, splitter, impl::eDiscardEmptyEntries, strip_chars ) ;
	}

	std::string join( std::vector< std::string > const& strings, std::string joiner ) {
		std::string result ;
		for( std::size_t i = 0; i < strings.size(); ++i ) {
			if( i > 0 ) {
				result += joiner ;
			}
			result += strings[i] ;
		}
		return result ;
	}

	std::string strip( std::string string_to_strip, std::string chars ) {
		return impl::strip( string_to_strip, chars ) ;
	}

	std::vector< std::string > split_and_strip( std::string string_to_split, std::string splitter, std::string strip_chars ) {
		return impl::split_and_strip( string_to_split, splitter, impl::ePreserveEmptyEntries, strip_chars ) ;
	}

	std::string to_lower( std::string aString ) {
		for( std::string::iterator i = aString.begin(); i != aString.end(); ++i ) {
			if( *i >= 'A' && *i <= 'Z' ) {
				*i += 32 ;
			}
		}
		return aString ;
	}

	std::string to_upper( std::string aString ) {
		for( std::string::iterator i = aString.begin(); i != aString.end(); ++i ) {
			if( *i >= 'a' && *i <= 'z' ) {
				*i -= 32 ;
			}
		}
		return aString ;
	}

	namespace impl {
		bool is_alphabetic( char c) {
			return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ;
		}
	
		bool is_blank( char c) {
			return c == ' ' || c == '\n' ;
		}
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
