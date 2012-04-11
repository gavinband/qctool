
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef INTEGRATION_INTEGRATOR_HPP
#define INTEGRATION_INTEGRATOR_HPP

#include <exception>

namespace integration {
	class Integrator
	{
	protected:
		Integrator( double desired_error ) ;
		double desired_error() const { return m_desired_error ; }
	private:
		double m_desired_error ;
	} ;
	
	struct IntegrationError: public std::exception
	{
		char const* what() const throw() { return "IntegrationError" ; }
	} ;
	
	struct IntervalTooSmallToSubdivideError: public IntegrationError
	{
		IntervalTooSmallToSubdivideError( double xmin, double xmax ): m_xmin( xmin ), m_xmax( xmax ) {}
		char const* what() const throw() { return "IntervalTooSmallToSubdivideError" ; }
		double xmin() const { return m_xmin ; }
		double xmax() const { return m_xmax ; }
	private:
		double const m_xmin, m_xmax ;
	} ;

	struct IntegralIsNaNError: public IntegrationError
	{
		IntegralIsNaNError() {}
		char const* what() const throw() { return "IntegralIsNaNError" ; }
	} ;
}

#endif
