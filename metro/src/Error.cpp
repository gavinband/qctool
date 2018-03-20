
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include "metro/Error.hpp"

namespace metro {
	ModelFitError::ModelFitError( std::string const& caller, std::string const& reason ):
			m_caller( caller ),
			m_reason( reason )
	{}

	ModelFitError::~ModelFitError() throw() {}

	char const* ModelFitError::what() const throw() {
		return "ModelFitError" ;
	}

	std::string const& ModelFitError::get_caller() const {
		return m_caller ;
	}
	std::string const& ModelFitError::get_reason() const {
		return m_reason ;
	}
	
	std::string NumericalError::format_message() const {
		return std::string( what() ) + ": " + m_source + ": " + m_message ;
	}
}

