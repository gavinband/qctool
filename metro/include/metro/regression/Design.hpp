
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
				Matrix const& outcome, NonmissingnessMatrix const& phenotype_nonmissingness, std::vector< std::string > const& outcome_names,
				Matrix const& covariates, NonmissingnessMatrix const& covariate_nonmissingness, std::vector< std::string > const& covariate_names,
				std::vector< std::string > const& predictor_names,
				Transform transform = eIdentity,
				std::vector< int > const& interacting_covariates = std::vector< int >()
			) ;

			Design(
				Matrix const& outcome, NonmissingnessMatrix const& phenotype_nonmissingness, std::vector< std::string > const& outcome_names,
				Matrix const& covariates, NonmissingnessMatrix const& covariate_nonmissingness, std::vector< std::string > const& covariate_names,
				std::vector< std::string > const& predictor_names,
				Transform transform = eIdentity,
				std::vector< int > const& interacting_covariates = std::vector< int >()
			) ;
		
			Design( Design const& other ) ;

			int const number_of_uncertain_predictors() const { return m_number_of_predictors * ( 1 + m_design_matrix_interaction_columns.size() ); }

			std::string const& get_predictor_name( std::size_t i ) const ;
			std::vector< std::string > const& design_matrix_column_names() const ;
			std::string const& get_outcome_name( std::size_t i ) const ;

			Design& set_outcome(
				Matrix const& outcome,
				NonmissingnessMatrix const& nonmissingness,
				std::vector< std::string > const& names
			) ;

			// levels is a K x d matrix
			// and probabilities is N x K
			// where
			// N = number of samples
			// K = number of different possible levels of the predictors
			// and d = number of predictors.
			// The interpretation is that individual n has predictor levels given by the kth row of levels with probability probabilities(n,k)
			Design& set_predictor_levels(
				Matrix const& levels,
				Matrix const& probabilities,
				std::vector< metro::SampleRange > const& included_samples
			) ;

			Design& set_covariates(
				Matrix const& covariates,
				Matrix const& covariate_nonmissingness,
				int start_column = 0
			) ;

			Matrix const& matrix() const { return m_design_matrix ; }
			Design& set_predictor_level( int level ) ;
		
			Matrix const& get_predictor_level_probabilities() const { return m_predictor_level_probabilities ; }
			int const get_number_of_predictor_levels() const { return m_number_of_predictor_levels ; }
			// Matrix const& get_predictor_levels() const { return m_predictor_levels ; }
		
			Matrix const& outcome() const { return m_outcome ; }
			std::vector< metro::SampleRange > const& globally_included_samples() const { return m_globally_included_samples ; }
			std::vector< metro::SampleRange > const& per_predictor_included_samples() const { return m_predictor_included_samples ; }
			Vector const& sample_weights() const { return m_sample_weights ; }

			std::string get_summary() const ;
		private:
			Matrix m_outcome ;
			NonmissingnessMatrix m_nonmissing_outcome ;
			Matrix m_nonmissing_covariates ;
			int const m_number_of_predictors ;
			std::vector< std::string > m_outcome_names ;
			std::vector< std::string > const m_covariate_names ;
			std::vector< std::string > const m_predictor_names ;
			Transform const m_transform ;

			Vector m_sample_weights ;
			SampleRanges m_globally_included_samples ;
			SampleRanges m_predictor_included_samples ;
			Matrix m_design_matrix ;
			NonmissingnessMatrix m_design_matrix_nonmissingness ;
			std::vector< int > m_design_matrix_interaction_columns ;
			std::vector< std::string > m_design_matrix_column_names ;

			Matrix m_predictor_level_probabilities ;
			Matrix m_predictor_levels ;
			int m_number_of_predictor_levels ;
		
		private:
			SampleRanges compute_included_samples(
				NonmissingnessMatrix const& nonmissing_outcome,
				NonmissingnessMatrix const& nonmissing_covariates
			) const ;
			void calculate_design_matrix(
				int const,
				Matrix const& covariates, NonmissingnessMatrix const& covariate_nonmissingness,
				int const number_of_predictors,
				std::vector< int > const& interaction_cols,
				Matrix* result,
				NonmissingnessMatrix* result_nonmissingness,
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
