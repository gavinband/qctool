
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef SNPTEST_CASE_CONTROL_ASSOCIATION_LOG_LIKELIHOOD_HPP
#define SNPTEST_CASE_CONTROL_ASSOCIATION_LOG_LIKELIHOOD_HPP

#include <vector>
#include <memory>
#include <boost/noncopyable.hpp>
#include "Eigen/Core"
#include "metro/RegressionDesign.hpp"
#include "metro/case_control/LogLikelihood.hpp"

namespace metro {
	namespace case_control {
		/*
		* This class implements a log-likelihood for binomial logistic regression, allowing predictors
		* to take one of a finite set of values with associated probabilities. (A "missing data log-likelihood" ).
		* The outcome variables must be 0 or 1.
		*/
		struct LogisticRegressionLogLikelihood: public LogLikelihood
		{
		public:
			typedef RegressionDesign::Point Point ;
			typedef RegressionDesign::Vector Vector ;
			typedef RegressionDesign::RowVector RowVector ;
			typedef RegressionDesign::Matrix Matrix ;
			typedef Eigen::Block< Matrix > MatrixBlock ;
			typedef Eigen::Block< Matrix const > ConstMatrixBlock ;
			typedef boost::function< std::string( std::string const& predictor_name, int outcome_level ) > GetParameterName ;
		public:
			typedef std::auto_ptr< LogisticRegressionLogLikelihood > UniquePtr ;
			static UniquePtr create( RegressionDesign::UniquePtr ) ;
			
		public:
			LogisticRegressionLogLikelihood( RegressionDesign::UniquePtr ) ;

			RegressionDesign const& get_design() const { return *m_design ; }
			void set_predictor_levels(
				Matrix const& levels,
				Matrix const& probabilities,
				std::vector< metro::SampleRange > const& included_samples
			) ;
		
			void set_parameter_naming_scheme( GetParameterName ) ;
			std::string get_parameter_name( std::size_t i ) const ;
	
			LogLikelihood::IntegerMatrix identify_parameters() const ;
			int number_of_outcomes() const ;

			void evaluate_at( Point const& parameters, int const numberOfDerivatives = 2 ) ;

			Point const& get_parameters() const ;
			double get_value_of_function() const ;
			Vector get_value_of_first_derivative() const ;
			Matrix get_value_of_second_derivative() const ;
			
			std::string get_summary() const ;

		private:
			RegressionDesign::UniquePtr m_design ;
			GetParameterName m_get_parameter_name ;
			Point m_parameters ;
			enum State {
				e_Uncomputed = 0,
				e_ComputedFunction = 1,
				e_Computed1stDerivative = 2,
				e_Computed2ndDerivative = 4
			} ;
			uint32_t m_state ;
			Matrix m_outcome_probabilities ;
			Matrix m_A ;
			Matrix m_B ;
			Vector m_temp ;
			Matrix m_first_derivative_terms ;

			double m_value_of_function ;
			Vector m_value_of_first_derivative ;
			Matrix m_value_of_second_derivative ;

		private:
			// Evaluate
			void evaluate_at_impl( Point const& parameters, std::vector< metro::SampleRange > const& included_samples, int const numberOfDerivatives ) ;
			// Calculate the probability of outcome given the genotype, parameters, and covariates.
			Vector evaluate_mean_function( Vector const& linear_combinations, Vector const& outcomes ) const ;
			// Calculate matrix of probabilities of outcome per genotype, given the parameters.
			void calculate_outcome_probabilities( Vector const& parameters, Vector const& phenotypes, Matrix* result ) const ;
			void compute_value_of_function( Matrix const& V, std::vector< metro::SampleRange > const& included_samples ) ;
			void compute_value_of_first_derivative( Matrix const& A, std::vector< metro::SampleRange > const& included_samples, Matrix* B ) ;
			void compute_value_of_second_derivative( Matrix const& B, std::vector< metro::SampleRange > const& included_samples ) ;
		} ;
	}
}

#endif
