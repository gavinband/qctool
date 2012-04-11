
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_STRING_UTILS_SUBSTITUTE_HPP
#define GENFILE_STRING_UTILS_SUBSTITUTE_HPP

#include <string>

namespace genfile {
	namespace string_utils {
		std::string substitute( std::string string_to_change, std::string const& what, std::string const& replacement ) ;
	}
}

#endif
