
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef SNPTEST_CASE_CONTROL_MULTINOMIAL_LOG_LIKELIHOOD_HPP
#define SNPTEST_CASE_CONTROL_MULTINOMIAL_LOG_LIKELIHOOD_HPP

#include <vector>
#include <memory>
#include <boost/noncopyable.hpp>
#include "Eigen/Core"
#include "metro/regression/Design.hpp"
#include "metro/regression/LogLikelihood.hpp"

namespace metro {
	namespace regression {
		/*
		* This class implements a log-likelihood for multinomial logistic regression, allowing predictors
		* to take one of a finite set of values with associated probabilities. (A "missing data log-likelihood" ).
		* The outcome variables must form a contiguous set of numbers in the range 0,1,...,M
		* for some positive integer M > 0.
		*
		* Notation:
		* N - number of samples
		* L - number of predictor levels
		* M - number of outcome levels not counting the baseline.
		* 
		*/
		struct MultinomialRegressionLogLikelihood: public LogLikelihood
		{
		public:
			typedef regression::Design::Vector Vector ;
			typedef regression::Design::RowVector RowVector ;
			typedef regression::Design::Matrix Matrix ;
			typedef Eigen::Block< Matrix > MatrixBlock ;
			typedef Eigen::Block< Matrix const > ConstMatrixBlock ;
			typedef Eigen::PermutationMatrix< Eigen::Dynamic, Eigen::Dynamic > PermutationMatrix ;
			typedef boost::function< std::string( std::string const& predictor_name, int outcome_level ) > GetParameterName ;
		public:
			typedef std::auto_ptr< MultinomialRegressionLogLikelihood > UniquePtr ;
			static UniquePtr create( regression::Design::UniquePtr ) ;
			
		public:
			MultinomialRegressionLogLikelihood( regression::Design::UniquePtr ) ;

			regression::Design& design() const { return *m_design ; }
			void set_predictor_levels( Matrix const& levels, Matrix const& probabilities, std::vector< metro::SampleRange > const& included_samples ) ;
		
			void set_parameter_naming_scheme( GetParameterName ) ;
			int number_of_parameters() const ;
			int number_of_outcomes() const ;
			std::string get_parameter_name( std::size_t i ) const ;		
			IntegerMatrix identify_parameters() const ;
			
			void evaluate_at( Vector const& parameters, int const numberOfDerivatives = 2 ) ;
			Vector const& parameters() const ;
			double get_value_of_function() const ;
			Vector get_value_of_first_derivative() const ;
			Matrix get_value_of_second_derivative() const ;

			std::string get_summary() const ;

		private:
			regression::Design::UniquePtr m_design ;
			int const m_number_of_samples ;
			GetParameterName m_get_parameter_name ;
			Vector m_parameter_vector ;
			Matrix m_parameter_matrix ;
			IntegerMatrix m_parameter_identity ;
			std::size_t m_numberOfDerivativesComputed ;

			// rearranger.  This matrix rearranges columns of an NxLM matrix from L blocks of NxM to M blocks of NxL.
			PermutationMatrix m_outcome_wise_to_predictor_wise_rearranger ;
			// Psi matrix, (NxLM).  This represents phenotype data as indicator 1's in each row.
			Matrix m_psi ;
			// F matrix (Nx(L(M+1))).  Stores outcome probabilities.
			Matrix m_F ;
			// A matrix (NxL).  Stores outcome probabilities times predictor probability, renormalised.
			Matrix m_A ;
			// B matrix (Nx(LM)).  Used in the computation of 1st and 2nd derivatives.
			Matrix m_B ;
			// C matrix (Nx(LM^2)).  Used in the computation of 2nd derivative.
			Matrix m_C ;
			// temp matrix.  Stores the terms that are summed to make the 1st derivative
			Matrix m_first_derivative_terms ;
			// temp matrix.  Stores rows of the design matrix tensor squares.
			mutable Matrix m_design_matrix_tensor_square_rows ;

			double m_value_of_function ;
			Vector m_value_of_first_derivative ;
			Matrix m_value_of_second_derivative ;

		private:
			
			void compute_psi(
				Matrix const& outcome,
				int const number_of_levels,
				std::vector< metro::SampleRange > const& included_samples,
				Matrix* result
			) const ;
			void compute_rearranger( int const number_of_predictor_levels, PermutationMatrix* result ) const ;
			IntegerMatrix compute_parameter_identity() const ;

			void evaluate_at_impl(
				Vector const& parameters,
				int const numberOfDerivatives
			) ;
			void compute_F( Matrix const& parameters, Matrix* result ) const ;
			void rearrange_F( Matrix* result ) const ;
			void compute_A_and_function_value(
				Matrix const& F,
				Matrix const& predictor_probabilities,
				std::vector< metro::SampleRange > const& included_samples,
				Matrix* result,
				double* value_of_function
			) const ;
			void compute_B(
				Matrix const& A,
				Matrix const& F,
				Matrix const& psi,
				Matrix* result
			) const ;
			void compute_value_of_first_derivative(
				Matrix const& B,
				std::vector< metro::SampleRange > const& included_samples,
				Vector* result,
				Matrix* terms
			) const ;
			// void compute_C( Matrix const& B, Matrix const& F, Matrix const& Gamma, Matrix* result ) const ;
			void compute_C(
				Matrix const& B,
				Matrix const& F,
				Matrix const& Gamma,
				std::vector< metro::SampleRange > const& included_samples,
				Matrix* result
			) const ;
			void compute_value_of_second_derivative(
				Matrix const& B,
				Matrix const& C,
				std::vector< metro::SampleRange > const& included_samples,
				Matrix* result
			) const ;
		} ;
	}
}

#endif
