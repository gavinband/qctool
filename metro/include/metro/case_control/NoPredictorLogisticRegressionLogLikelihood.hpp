
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef SNPTEST_CASE_CONTROL_ASSOCIATION_NO_PREDICTOR_LOG_LIKELIHOOD_HPP
#define SNPTEST_CASE_CONTROL_ASSOCIATION_NO_PREDICTOR_LOG_LIKELIHOOD_HPP

#include <vector>
#include <memory>
#include <boost/noncopyable.hpp>
#include <boost/function.hpp>
#include "Eigen/Core"
#include "metro/case_control/RegressionDesign.hpp"
#include "metro/case_control/LogLikelihood.hpp"

namespace metro {
	namespace case_control {
		struct NoPredictorLogisticRegressionLogLikelihood: public LogLikelihood
		{
		public:
			typedef Eigen::VectorXd Point ;
			typedef Eigen::VectorXd Vector ;
			typedef Eigen::RowVectorXd RowVector ;
			typedef Eigen::MatrixXd Matrix ;
			typedef boost::function< std::string( std::size_t ) > GetPredictorNames ;

		public:
			typedef std::auto_ptr< NoPredictorLogisticRegressionLogLikelihood > UniquePtr ;
			static UniquePtr create( RegressionDesign::UniquePtr ) ;
			
		public:
			NoPredictorLogisticRegressionLogLikelihood( RegressionDesign::UniquePtr ) ;

			RegressionDesign const& get_design() const { return *m_design ; }
			void set_predictor_levels( Matrix const& levels, Matrix const& probabilities, std::vector< metro::SampleRange > const& included_samples ) ;
			int number_of_outcomes() const ;
		
			std::string get_parameter_name( std::size_t i ) const ;	
			LogLikelihood::IntegerMatrix identify_parameters() const ;
			
			void evaluate_at( Point const& parameters, int const numberOfDerivatives = 2 ) ;
			Point const& get_parameters() const ;
			double get_value_of_function() const ;
			Vector get_value_of_first_derivative() const ;
			Matrix get_value_of_second_derivative() const ;
			
			std::string get_summary() const ;

		private:
			double const m_threshhold_weight ;
			RegressionDesign::UniquePtr m_design ;
			std::vector< metro::SampleRange > m_included_samples ;
			
			Point m_parameters ;
			Matrix m_design_matrix ;
			Matrix m_design_matrix_row_tensor_squares ;
			Vector m_outcome_probabilities ;
			Vector m_V ;

			double m_value_of_function ;
			Vector m_value_of_first_derivative ;
			Matrix m_value_of_second_derivative ;

		private:
			std::vector< int > compute_exclusions() const ;
			void calculate_design_matrix( int const, Matrix const& covariates, int const number_of_predictors, Matrix* result ) const ;
			void compute_tensor_squares_of_rows( Matrix* result ) const ;

			void evaluate_at_impl( Point const& parameters, std::vector< metro::SampleRange > const& included_samples, int const numberOfDerivatives ) ;
			// Calculate the probability of outcome given the genotype, parameters, and covariates.
			void evaluate_mean_function( Vector const& linear_combinations, Vector const& outcomes, Vector* result ) const ;
			// Calculate matrix of probabilities of outcome per genotype, given the parameters.
			void calculate_outcome_probabilities( Vector const& parameters, Vector const& phenotypes, Vector* result ) const ;
			void compute_value_of_function( std::vector< metro::SampleRange > const& included_samples ) ;
			void compute_value_of_first_derivative( std::vector< metro::SampleRange > const& included_samples ) ;
			void compute_value_of_second_derivative( std::vector< metro::SampleRange > const& included_samples ) ;
		} ;
	}
}

#endif
