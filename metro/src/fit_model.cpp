
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <utility>
#include <Eigen/Core>
#include <limits>
#include <cassert>
#include "metro/regression/LogLikelihood.hpp"
#include "metro/ModifiedCholesky.hpp"
#include "metro/ModifiedNewtonRaphson.hpp"
#include "metro/fit_model.hpp"

#define USE_MODIFIED_NR 1
#define DEBUG 1

namespace metro {
	Stepper::~Stepper() {}

	std::pair< bool, int > fit_model(
		metro::SmoothFunction& ll,
		std::string const& model_name,
		Eigen::VectorXd const& starting_point,
		Stepper& stopping_condition,
		std::vector< std::string >* comments
	) {
		typedef metro::regression::Design::Matrix Matrix ;
		// Fit null model
		assert( starting_point.size() == ll.number_of_parameters() ) ;
		ll.evaluate_at( starting_point ) ;

		bool success = true ;
		int iterations = 0 ;
		// Compute negative definiteness condition number.
		Eigen::SelfAdjointEigenSolver< Matrix > eigenSolver( ll.get_value_of_second_derivative() ) ;
		if( eigenSolver.eigenvalues().maxCoeff() > 0 ) {
			success = false ;
			if( comments ) {
				comments->push_back( model_name + ":not-negative-definite" ) ;
			}
		} else {
#if 1
		Eigen::VectorXd parameters = metro::find_maximum_by_modified_newton_raphson_with_line_search(
			ll, starting_point, stopping_condition
		) ;
#else
			metro::Snptest25StoppingCondition< metro::regression::LogLikelihood > stopping_condition(
				ll,
				0.01,
				100,
				0
				//&get_ui_context().logger()
			) ;

			Eigen::ColPivHouseholderQR< Matrix > solver ;
			Eigen::VectorXd parameters = maximise_by_newton_raphson( ll, starting_point, stopping_condition, solver ) ;
#endif

			if( stopping_condition.diverged() ) {
				if( comments ) {
					comments->push_back( model_name + ":model_fit_error:failed_to_converge_to_mle" ) ;
				}
			}

			success = !stopping_condition.diverged() ;
			iterations = stopping_condition.number_of_iterations() ;
		}

		return std::make_pair( success, iterations ) ;
	}
}