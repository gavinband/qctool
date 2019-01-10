
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
	struct SmoothFunction: public boost::noncopyable {
	public:
		typedef std::unique_ptr< SmoothFunction > UniquePtr ;
		typedef double Scalar ;
		typedef Eigen::VectorXd Vector ;
		typedef Eigen::MatrixXd Matrix ;
	public:
		virtual ~SmoothFunction() {}
		virtual int number_of_parameters() const = 0 ;
		virtual void evaluate_at( Vector const& parameters, int const numberOfDerivatives = 2 ) = 0 ;
		virtual Scalar get_value_of_function() const = 0 ;
		virtual Vector get_value_of_first_derivative() const = 0 ;
		virtual Matrix get_value_of_second_derivative() const = 0 ;
		virtual std::string get_summary() const = 0 ;
	} ;
}

#endif
