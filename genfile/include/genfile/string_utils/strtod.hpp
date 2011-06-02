#ifndef GENFILE_STRING_UTILS_STRTOD_HPP
#define GENFILE_STRING_UTILS_STRTOD_HPP

#include <string>
#include "genfile/string_utils/slice.hpp"

namespace genfile {
	namespace string_utils {
		// Convert string to a double.
		// Throw a StringConversionError if not possible.
		double strtod( slice const& string ) ;
	}
}

#endif
