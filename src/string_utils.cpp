#include <cassert>
#include <vector>
#include <string>
#include "string_utils.hpp"

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

std::vector< std::string > split_discarding_empty_entries( std::string string_to_split, std::string splitter ) {
	return impl::split_and_strip( string_to_split, splitter, impl::eDiscardEmptyEntries, "" ) ;
}

std::string strip( std::string string_to_strip, std::string chars ) {
	return impl::strip( string_to_strip, chars ) ;
}

std::vector< std::string > split_and_strip( std::string string_to_split, std::string splitter, std::string strip_chars ) {
	return impl::split_and_strip( string_to_split, splitter, impl::ePreserveEmptyEntries, strip_chars ) ;
}

std::vector< std::string > split_and_strip_discarding_empty_entries( std::string string_to_split, std::string splitter, std::string strip_chars ) {
	return impl::split_and_strip( string_to_split, splitter, impl::eDiscardEmptyEntries, strip_chars ) ;	
}

std::string to_lower( std::string aString ) {
	for( std::string::iterator i = aString.begin(); i != aString.end(); ++i ) {
		if( *i >= 'A' && *i <= 'Z' ) {
			*i += 40 ;
		}
	}
	return aString ;
}

std::string to_upper( std::string aString ) {
	for( std::string::iterator i = aString.begin(); i != aString.end(); ++i ) {
		if( *i >= 'a' && *i <= 'z' ) {
			*i -= 40 ;
		}
	}
	return aString ;
}
