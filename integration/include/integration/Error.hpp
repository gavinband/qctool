//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef INTEGRATION_ERROR_HPP
#define INTEGRATION_ERROR_HPP

#include <exception>

namespace integration {
	struct NumericalError: public std::exception 
	{
		virtual ~NumericalError() throw() {}

		NumericalError( std::string const& source, std::string const& message = "" ):
			m_source( source ),
			m_message( message )
		{}
		
		NumericalError( NumericalError const& other ):
			m_source( other.m_source ),
			m_message( other.m_message )
		{}

		const char* what() const throw() { return "integration:NumericalError" ; }
		virtual std::string format_message() const ;
		std::string const& source() const { return m_source ; }
		std::string const& message() const { return m_message ; }
	private:
		std::string const m_source ;
		std::string const m_message ;
	} ;
}

#endif
