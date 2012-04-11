
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "genfile/Error.hpp"

namespace genfile {
	std::string BadArgumentError::format_message() const {
		return "In argument(s) " + arguments() + " to function " + function() ;
	}

	std::string BadArgumentWithMessageError::format_message() const {
		return "In argument(s) " + arguments() + " to function " + function() + ": " + m_message ;
	}
}
