
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_LIKELIHOOD_SMOOTH_FUNCTION_HPP
#define METRO_LIKELIHOOD_SMOOTH_FUNCTION_HPP

#include <memory>
#include <boost/noncopyable.hpp>
#include <Eigen/Core>
#include "metro/DataSubset.hpp"

namespace metro {
	// This class represents a smooth function, i.e. a function that is
	// twice differentiable.
	// SmoothFunctions are intended to be stateful - they are properly thought
	// of as a function evaluated at a specific point.  This design choice
	// was made because it is typically the case that a large amount of computation can be
	// shared between derivatives.
	
	// a great deal of comp
	// The intentino 
	// This choice means that the accessor functions, named "get_value[...]",
	// are const and return references, while the 
	struct SmoothFunction: public boost::noncopyable {
	public:
		typedef std::unique_ptr< SmoothFunction > UniquePtr ;
		typedef double Scalar ;
		typedef Eigen::VectorXd Vector ;
		typedef Eigen::MatrixXd Matrix ;
	public:
		virtual ~SmoothFunction() {}
		virtual int number_of_parameters() const = 0 ;
		virtual Vector const& parameters() const = 0 ;

		// Evaluate the function at a given parameter value
		virtual void evaluate_at( Vector const& parameters, int const numberOfDerivatives = 2 ) = 0 ;
		// Re-evaluate the function at a previously set parameter value, previously
		// set by a call to evaluate_at().
		// This allows to evaluate a specified number of derivatives.
		virtual void evaluate( int const numberOfDerivatives = 2 ) = 0 ;
		
		virtual Scalar get_value_of_function() const = 0 ;
		virtual Vector get_value_of_first_derivative() const = 0 ;
		virtual Matrix get_value_of_second_derivative() const = 0 ;
		virtual std::string get_summary() const = 0 ;
	} ;
}

#endif
