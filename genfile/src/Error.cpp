#include "genfile/Error.hpp"

namespace genfile {
	std::string BadArgumentError::format_message() const {
		return "In argument(s) " + arguments() + " to function " + function() ;
	}

	std::string BadArgumentWithMessageError::format_message() const {
		return "In argument(s) " + arguments() + " to function " + function() + ": " + m_message ;
	}
}
