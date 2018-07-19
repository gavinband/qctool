
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef SNPTEST_ERROR_HPP
#define SNPTEST_ERROR_HPP

#include <exception>
#include <string>

namespace metro {
	struct LimitExceededError: public std::exception {
	public:
		LimitExceededError( std::string const& source, std::string const& limit, std::string const& reason ):
			m_source( source ),
			m_limit( limit ),
			m_reason( reason )
		{}
		~LimitExceededError() throw() {}
		char const* what() const throw() { return "core::LimitExceededError" ; }
		std::string const& source() const { return m_source ; }
		std::string const& limit_description() const { return m_limit ; }
		std::string const& reason() const { return m_reason ; }
	private:
		std::string const m_source ;
		std::string const m_limit ;
		std::string const m_reason ;
	} ;

	struct ModelFitError: public std::exception {
		ModelFitError( std::string const& caller, std::string const& reason ) ;
		~ModelFitError() throw() ;

		char const* what() const throw() ;
		std::string const& get_caller() const ;
		std::string const& get_reason() const ;
	private:
		std::string const m_caller ;
		std::string const m_reason ;
	} ;
	
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
