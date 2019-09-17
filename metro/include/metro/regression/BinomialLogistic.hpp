
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
#include "metro/regression/Design.hpp"
#include "metro/regression/LogLikelihood.hpp"

namespace metro {
	namespace regression {
		/*
		* This class implements a log-likelihood for binomial logistic regression, allowing predictors
		* to take one of a finite set of values with associated probabilities. (A "missing data log-likelihood" ).
		* The outcome must be a two-column matrix containing positive integers; in the binomial formulation
		* in terms of 'successes/failures', the 1st column is the number of failure (baseline) outcomes
		* and the 2nd column the number of successful (non-baseline) outcomes.
		*/
		struct BinomialLogistic: public LogLikelihood
		{
		public:
			typedef regression::Design::Point Point ;
			typedef regression::Design::Vector Vector ;
			typedef regression::Design::RowVector RowVector ;
			typedef regression::Design::Matrix Matrix ;
			typedef Eigen::Block< Matrix > MatrixBlock ;
			typedef Eigen::Block< Matrix const > ConstMatrixBlock ;
			typedef boost::function< std::string( std::string const& predictor_name, std::string const& outcome_name ) > GetParameterName ;
		public:
			typedef std::auto_ptr< BinomialLogistic > UniquePtr ;
			static UniquePtr create( Design& ) ;
			static UniquePtr create( Design::UniquePtr ) ;
			static UniquePtr create( Design&, std::vector< metro::SampleRange > ) ;
			static UniquePtr create( Design::UniquePtr, std::vector< metro::SampleRange > ) ;
			BinomialLogistic( Design& ) ;
			BinomialLogistic( Design&, std::vector< metro::SampleRange > included_samples ) ;
			BinomialLogistic( Design::UniquePtr ) ;
			BinomialLogistic( Design::UniquePtr, std::vector< metro::SampleRange > included_samples ) ;
			~BinomialLogistic() ;
			
			regression::Design& design() const { return *m_design ; }
			void set_parameter_naming_scheme( GetParameterName ) ;
			std::string get_parameter_name( std::size_t i ) const ;
	
			int number_of_outcomes() const ;
			int number_of_parameters() const { return m_parameters.size() ; }
			IntegerMatrix identify_parameters() const ;

			void evaluate_at( Point const& parameters, int const numberOfDerivatives = 2 ) ;
			void evaluate( int const numberOfDerivatives = 2 ) ;

			Point const& parameters() const ;
			double get_value_of_function() const ;
			Vector get_value_of_first_derivative() const ;
			Matrix get_value_of_second_derivative() const ;
			
			std::string get_summary() const ;

		private:
			regression::Design* m_design ;
			bool m_design_owned ;
			GetParameterName m_get_parameter_name ;
			Point m_parameters ;
			std::vector< metro::SampleRange > m_included_samples ;
			std::vector< metro::SampleRange > m_evaluated_samples ;
			std::vector< metro::SampleRange > m_storage_samples ;
			int m_number_of_stored_samples ;
			int m_numberOfDerivativesComputed ;
			int m_numberOfCDLDerivativesComputed ;

			Matrix m_outcome_probabilities ;
			
			// Given outcome, likelihood function is of the form f(x^t theta) where
			// theta are the parameters (currently we don't handle the case e.g. normal
			// linear regression where there are additional parameters).  
			// Therefore derivatives are of the form (Df)|(x^t theta).x^t
			// And D^2 f(x^t theta) . x^t
			
			Matrix m_hx ; // coefficient for each level x in complete data likelihood for each sample.
			Matrix m_normalisedDhx ; // coefficient of x in 1st derivative of mean function.
			Matrix m_normalisedDdhx ; // coefficient of x^t ⊗ x in 2nd derivative of mean function.
			Matrix m_f1 ; // temp storage used in likelihood and derivative computations

			Matrix m_first_derivative_terms ;
			double m_value_of_function ;
			Vector m_value_of_first_derivative ;
			Matrix m_value_of_second_derivative ;

		private:
			// Evaluate
			void evaluate_at_impl( Point const& parameters, std::vector< metro::SampleRange > const& included_samples, int const numberOfDerivatives ) ;

			void setup_storage() ;

			// Compute complete data likelihood, and its derivatives, for each predictor level
			// and each sample
			virtual void compute_complete_data_likelihood_and_derivatives(
				Point const& parameters,
				Matrix* fx,
				Matrix* dfx,
				Matrix* ddfx,
				int const numberOfDerivatives
			) ;

			// Calculate matrix of probabilities of outcome per genotype, given the parameters.
			void calculate_outcome_probabilities( Vector const& parameters, Matrix const& phenotypes, Matrix* result ) const ;
			void compute_value_of_loglikelihood(
				Matrix const& hx,
				double* result
			) ;
			void compute_value_of_first_derivative(
				Matrix const& normalisedDhx,
				Matrix* result_terms,
				Vector* result
			) ;
			void compute_value_of_second_derivative(
				Matrix const& first_derivative_terms,
				Matrix const& normalisedDdhx,
				Matrix* result
			) ;
		} ;
	}
}

#endif
