
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef SNPTEST_REGRESSION_LOG_LIKELIHOOD_HPP
#define SNPTEST_REGRESSION_LOG_LIKELIHOOD_HPP

#include <vector>
#include <memory>
#include <boost/noncopyable.hpp>
#include "Eigen/Core"
#include "metro/regression/Design.hpp"

namespace metro {
	namespace regression {
		/*
		* Base class for classes implementing regression log-likelihoods
		*/
		struct LogLikelihood: public boost::noncopyable
		{
		public:
			typedef std::auto_ptr< LogLikelihood > UniquePtr ;
			typedef regression::Design::Point Point ;
			typedef regression::Design::Vector Vector ;
			typedef regression::Design::RowVector RowVector ;
			typedef regression::Design::Matrix Matrix ;
			typedef Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic > IntegerMatrix ;
		public:
			virtual ~LogLikelihood() ;
			virtual regression::Design& design() const = 0 ;
			virtual std::string get_parameter_name( std::size_t i ) const = 0 ;

			// Return a lx2 matrix identifying the l parameters.
			// For each parameter the row contains:
			// column 0: the index of the outcome column.
			// column 1: the index of the design matrix column.
			virtual IntegerMatrix identify_parameters() const = 0 ;
			virtual int number_of_outcomes() const = 0 ;

			virtual void evaluate_at( Point const& parameters, int const numberOfDerivatives = 2 ) = 0 ;
			virtual Point const& get_parameters() const = 0 ;
			virtual double get_value_of_function() const = 0 ;
			virtual Vector get_value_of_first_derivative() const = 0 ;
			virtual Matrix get_value_of_second_derivative() const = 0 ;
			
			virtual std::string get_summary() const = 0 ;
		} ;
	}
}

#endif
