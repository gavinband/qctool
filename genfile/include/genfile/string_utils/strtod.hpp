
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

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
