//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <exception>
#include <string>
#include "integration/Error.hpp"

namespace integration {
	std::string NumericalError::format_message() const {
		return std::string( what() ) + ": " + m_source + ": " + m_message ;
	}
}

