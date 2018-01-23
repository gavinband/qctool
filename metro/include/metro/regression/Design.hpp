
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef SNPTEST_CASE_CONTROL_REGRESSION_DESIGN_HPP
#define SNPTEST_CASE_CONTROL_REGRESSION_DESIGN_HPP

#include <vector>
#include <iostream>
#include <memory>
#include <boost/noncopyable.hpp>
#include <boost/function.hpp>
#include "Eigen/Core"
#include "metro/SampleRange.hpp"

namespace metro {
	namespace regression {
		struct Design: public boost::noncopyable
		{
			typedef std::auto_ptr< Design > UniquePtr ;
		public:
			typedef Eigen::VectorXd Point ;
			typedef Eigen::VectorXd Vector ;
			typedef Eigen::RowVectorXd RowVector ;
			typedef Eigen::MatrixXd Matrix ;
			//typedef Eigen::Matrix< char, Eigen::Dynamic, Eigen::Dynamic > NonmissingnessMatrix ;
			typedef Matrix NonmissingnessMatrix ;
			typedef Eigen::Block< Matrix > MatrixBlock ;
			typedef Eigen::Block< Matrix const > ConstMatrixBlock ;
			typedef std::vector< metro::SampleRange > SampleRanges ;
			typedef boost::function< std::string const& ( int parameter_i ) > GetPredictorNames ;
			enum Transform { eIdentity = 0, eMeanCentre = 1 } ;

		public:
		
			static UniquePtr create(
				Matrix const& outcome, SampleRanges const& nonmissing_outcome, std::vector< std::string > const& outcome_names,
				Matrix const& covariates, SampleRanges const& nonmissing_covariates, std::vector< std::string > const& covariate_names,
				std::vector< std::string > const& predictor_names,
				Transform transform = eIdentity,
				std::vector< int > const& interacting_covariates = std::vector< int >()
			) ;

			Design(
				Matrix const& outcome, SampleRanges const& nonmissing_outcome, std::vector< std::string > const& outcome_names,
				Matrix const& covariates, SampleRanges const& nonmissing_covariates, std::vector< std::string > const& covariate_names,
				std::vector< std::string > const& predictor_names,
				Transform transform = eIdentity,
				std::vector< int > const& interacting_covariates = std::vector< int >()
			) ;
		
			Design( Design const& other ) ;

			int const number_of_predictors() const { return m_predictor_names.size() ; }
			int const number_of_interaction_terms() const { return number_of_predictors() * m_design_matrix_interaction_columns.size() ; }

			std::string const& get_predictor_name( std::size_t i ) const ;
			std::vector< std::string > const& design_matrix_column_names() const ;
			std::string const& get_outcome_name( std::size_t i ) const ;

			Design& set_outcome(
				Matrix const& outcome,
				SampleRanges const& nonmissingness,
				std::vector< std::string > const& names
			) ;

			// levels is a K x d matrix
			// and probabilities is N x K
			// where
			// N = number of samples
			// K = number of different possible levels of the predictors
			// and d = number of predictors.
			// The interpretation is that individual n has predictor levels given by the kth row of levels with probability probabilities(n,k)
			Design& set_predictors(
				Matrix const& levels,
				Matrix const& probabilities,
				SampleRanges const& included_samples
			) ;

			Matrix const& matrix() const { return m_design_matrix ; }
			Design& set_predictor_level( int level ) ;
		
			Matrix const& get_predictor_level_probabilities() const { return m_predictor_level_probabilities ; }
			int const get_number_of_predictor_levels() const { return m_predictor_level_probabilities.cols() ; }
			// Matrix const& get_predictor_levels() const { return m_predictor_levels ; }
		
			Matrix const& outcome() const { return m_outcome ; }

			SampleRanges const& nonmissing_samples() const { return m_nonmissing_samples ; }

			std::string get_summary() const ;
		private:
			// Outcome variables
			Matrix m_outcome ;
			std::vector< std::string > m_outcome_names ;
			SampleRanges m_nonmissing_outcome ;

			// Predictor variables
			Matrix m_predictor_level_probabilities ;
			Matrix m_predictor_levels ;
			std::vector< std::string > const m_predictor_names ;
			SampleRanges m_nonmissing_predictors ;

			// Covariate variables
			std::vector< std::string > const m_covariate_names ;
			SampleRanges m_nonmissing_covariates ;

			// Design matrix column transform - e.g. mean centre.
			Transform const m_transform ;

			// The current design matrix that holds current predictor values and covariates.
			Matrix m_design_matrix ;

			// SampleRange that tracks samples that have no nonmissing data
			// in outcome, covariates, or predictors
			// Invariant: this the intersection of nonmissing samples implied
			// by m_nonmissing_outcome, m_nonmissing_predictors, and m_nonmissing_covariates.
			SampleRanges m_nonmissing_samples ;
			
			std::vector< int > m_design_matrix_interaction_columns ;
			std::vector< std::string > m_design_matrix_column_names ;

		private:
			SampleRanges compute_nonmissing_samples(
				SampleRanges const& nonmissing_outcome,
				SampleRanges const& nonmissing_predictors,
				SampleRanges const& nonmissing_covariates
			) const ;
			void calculate_design_matrix(
				int const,
				Matrix const& covariates,
				std::vector< std::string > const& predictor_names,
				std::vector< int > const& interaction_cols,
				Matrix* result,
				std::vector< int >* design_matrix_interaction_cols,
				std::vector< std::string > * design_matrix_column_names
			) const ;
			
			void mean_centre_predictor_levels(
				Matrix const& probs,
				Matrix* predictor_levels,
				SampleRanges const& included_samples
			) const ;
		} ;
	}
}

#endif
