#include <string>
#include <sstream>
#include <iomanip>
#include "genfile/string_utils/hex.hpp"

namespace genfile {
	namespace string_utils {
		std::string to_hex( std::string const& str ) {
			std::ostringstream o ;
			for( std::size_t i = 0; i < str.size(); ++i ) {
				if( i % 4 == 0 )
					o << "|" ;
				o << std::hex << std::setw(2) << std::setfill('0') << static_cast<int> (str[i]) ;
			}
			return o.str() ;
		}
	}
}