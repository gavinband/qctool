#include <string>
#include "genfile/string_utils/substitute.hpp"

namespace genfile {
	namespace string_utils {
		std::string substitute( std::string string_to_change, std::string const& what, std::string const& replacement ) {
			std::size_t pos = 0 ;
			while( ( pos = string_to_change.find( what, pos ) ) != std::string::npos ) {
				string_to_change.replace( pos, what.size(), replacement ) ;
				pos = pos + replacement.size() ;
			}
			return string_to_change ;
		}
	}
}

