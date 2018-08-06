
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
		* This class implements a log-likelihood for Bernoulli logistic regression, allowing predictors
		* to take one of a finite set of values with associated probabilities. (A "missing data log-likelihood" ).
		* The outcome must be a two-column matrix containing 0s and 1s; a 1 in the 2nd column indicates
		* a non-baseline outcome and a 1 in the 1st column indicates a baseline outcome.
		*/
		struct Logistic: public LogLikelihood
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
			typedef std::auto_ptr< Logistic > UniquePtr ;
			static UniquePtr create( Design& ) ;
			static UniquePtr create( Design::UniquePtr ) ;
			Logistic( Design& ) ;
			Logistic( Design::UniquePtr ) ;
			~Logistic() ;
			
			regression::Design& design() const { return *m_design ; }
			void set_parameter_naming_scheme( GetParameterName ) ;
			std::string get_parameter_name( std::size_t i ) const ;
	
			IntegerMatrix identify_parameters() const ;
			int number_of_outcomes() const ;

			void evaluate_at( Point const& parameters, int const numberOfDerivatives = 2 ) ;

			Point const& get_parameters() const ;
			double get_value_of_function() const ;
			Vector get_value_of_first_derivative() const ;
			Matrix get_value_of_second_derivative() const ;
			
			std::string get_summary() const ;

		private:
			regression::Design* m_design ;
			bool m_design_owned ;
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
			Vector evaluate_mean_function( Vector const& linear_combinations, Matrix const& outcomes ) const ;
			// Calculate matrix of probabilities of outcome per genotype, given the parameters.
			void calculate_outcome_probabilities( Vector const& parameters, Matrix const& phenotypes, Matrix* result ) const ;
			void compute_value_of_function( Matrix const& V, std::vector< metro::SampleRange > const& included_samples ) ;
			void compute_value_of_first_derivative( Matrix const& A, std::vector< metro::SampleRange > const& included_samples, Matrix* B ) ;
			void compute_value_of_second_derivative( Matrix const& B, std::vector< metro::SampleRange > const& included_samples ) ;
		} ;
	}
}

#endif
