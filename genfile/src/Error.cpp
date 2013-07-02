
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "genfile/Error.hpp"

namespace genfile {
	std::string BadArgumentError::format_message() const {
		return "In argument(s) " + arguments() + " to function " + function() + ( message().size() > 0 ? ( ": " + message() ) : "" ) ;
	}

	std::string InputError::format_message() const {
		std::ostringstream ostr ;
		ostr << "An input error occured in source \"" << source() << "\"" ;
		if( m_message.size() > 0 ) {
			ostr << ": " << m_message ;
		}
		ostr << "." ;
		return ostr.str() ;
	}

	std::string MalformedInputError::format_message() const {
		std::ostringstream ostr ;
		ostr << "Source \"" << source() << "\" is malformed on line " << ( line() + 1 ) ;
		if( has_column() ) {
			ostr << ", column " << ( column() + 1 ) ;
		}
		if( message().size() > 0 ) {
			ostr << ": " << message() ;
		}
		ostr << "." ;
		return ostr.str() ;
	}
	
	std::string FormatUnsupportedError::format_message() const {
		return "Source \"" + source() + "\" is in format \"" + format() + "\", which is unsupported." ;
	}
	
	std::string DuplicateEntryError::format_message() const {
		std::ostringstream ostr ;
		ostr << "Source \"" << source() << "\" has a duplicate entry for key \"" << m_key << "\", variable \"" << m_variable << "\"" ;
		if( message().size() > 0 ) {
			ostr << ": " << message() ;
		}
		ostr << "." ;
		return ostr.str() ;
	}

	std::string DuplicateIndividualError::format_message() const {
		std::ostringstream ostr ;
		ostr << "Source \"" << source() << "\" has a duplicate entry for individual \"" << m_id_1 << "\"" ;
		if( message().size() > 0 ) {
			ostr << ": " << message() ;
		}
		ostr << "." ;
		return ostr.str() ;
	}
	
}
