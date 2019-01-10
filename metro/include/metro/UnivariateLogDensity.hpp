
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_DISTRIBUTIONS_UNIVARIATELOGDENSITY_HPP
#define METRO_DISTRIBUTIONS_UNIVARIATELOGDENSITY_HPP

#include <memory>

namespace metro {
	struct UnivariateLogDensity {
		typedef std::unique_ptr< UnivariateLogDensity > UniquePtr ;
		
		virtual ~UnivariateLogDensity() ;
		virtual void evaluate_at( double x ) = 0 ;
		virtual double get_value_of_function() const = 0 ;
		virtual double get_value_of_first_derivative() const = 0 ;
		virtual double get_value_of_second_derivative() const = 0 ;
		virtual std::string get_summary() const = 0 ;
	} ;
}

#endif
