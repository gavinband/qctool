#ifndef GENFILE_STRING_UTILS_SUBSTITUTE_HPP
#define GENFILE_STRING_UTILS_SUBSTITUTE_HPP

#include <string>

namespace genfile {
	namespace string_utils {
		std::string substitute( std::string string_to_change, std::string const& what, std::string const& replacement ) ;
	}
}

#endif
