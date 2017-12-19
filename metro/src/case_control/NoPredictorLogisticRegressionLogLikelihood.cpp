
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include <vector>
#include <utility>
#include <algorithm>
#include <boost/iterator/counting_iterator.hpp>
#include "metro/case_control/NoPredictorLogisticRegressionLogLikelihood.hpp"
#include "metro/intersect_ranges.hpp"
#include "genfile/Error.hpp"
#include "genfile/string_utils/string_utils.hpp"

// #define DEBUG_LOGLIKELIHOOD 1

namespace metro {
	namespace case_control {
		NoPredictorLogisticRegressionLogLikelihood::UniquePtr NoPredictorLogisticRegressionLogLikelihood::create(
			RegressionDesign::UniquePtr design
		) {
			return NoPredictorLogisticRegressionLogLikelihood::UniquePtr(
				new NoPredictorLogisticRegressionLogLikelihood( design )
			) ;
		}
		
		
		NoPredictorLogisticRegressionLogLikelihood::NoPredictorLogisticRegressionLogLikelihood(
			RegressionDesign::UniquePtr design
		):
			m_threshhold_weight( 0.1 ),
			m_design( design ),
			m_parameters( Vector::Zero( m_design->matrix().cols() ) )
		{
			assert( m_design.get() ) ;
			assert( m_design->number_of_uncertain_predictors() == 0 ) ;
			compute_tensor_squares_of_rows( &m_design_matrix_row_tensor_squares ) ;
		}

		std::string NoPredictorLogisticRegressionLogLikelihood::get_summary() const {
			using genfile::string_utils::to_string ;
			return "NoPredictorLogisticRegressionLogLikelihood( " + to_string( m_design->outcome().size() ) + " samples )" ;
		}

		void NoPredictorLogisticRegressionLogLikelihood::set_predictor_levels( Matrix const& levels, Matrix const& probabilities, std::vector< metro::SampleRange > const& included_samples ) {
			// No varying predictor, so there should be a single level with no columns.
			assert( levels.rows() == 1 ) ;
			assert( levels.cols() == 0 ) ;
			assert( probabilities.cols() == 1 ) ;

			m_design->set_sample_weights( probabilities ) ;
			m_included_samples = included_samples ;
		}
		
		// Construct a big matrix of size ( D * #of samples) x ( D * #of genotype levels ).
		// The block at position (D*sample, D*g) is the tensor square of the design matrix row
		// for the given sample and genotype level.
		void NoPredictorLogisticRegressionLogLikelihood::compute_tensor_squares_of_rows( Matrix* result ) const {
			Matrix const& design_matrix = m_design->matrix() ;
			int const D = design_matrix.cols() ;
			result->resize( D * design_matrix.rows(), D ) ;
			for( int sample = 0; sample < design_matrix.rows(); ++sample ) {
				result->block( D * sample, 0, D, D ) = design_matrix.row( sample ).transpose() * design_matrix.row( sample ) ;
			}
		}

		std::string NoPredictorLogisticRegressionLogLikelihood::get_parameter_name( std::size_t i ) const {
			return "PARAMETER_" + genfile::string_utils::to_string( i ) ;
		}

		LogLikelihood::IntegerMatrix NoPredictorLogisticRegressionLogLikelihood::identify_parameters() const {
			LogLikelihood::IntegerMatrix result = LogLikelihood::IntegerMatrix::Zero( m_design->matrix().cols(), 2 ) ;
			result.col(0).setConstant( 1.0 ) ;
			for( int i = 0; i < m_design->matrix().cols(); ++i ) {
				result( i, 1 ) = i ;
			}
			return result ;
		}

		int NoPredictorLogisticRegressionLogLikelihood::number_of_outcomes() const {
			return 2 ;
		}

		void NoPredictorLogisticRegressionLogLikelihood::evaluate_at( Point const& parameters, int const numberOfDerivatives ) {
			evaluate_at_impl( parameters, m_included_samples, numberOfDerivatives ) ;
		}

		void NoPredictorLogisticRegressionLogLikelihood::evaluate_at_impl(
			Point const& parameters,
			std::vector< metro::SampleRange > const& included_samples,
			int const numberOfDerivatives
		) {
			m_parameters = parameters ;
			calculate_outcome_probabilities( parameters, m_design->outcome(), &m_outcome_probabilities ) ;
			assert( numberOfDerivatives < 3 ) ;

#if DEBUG_LOGLIKELIHOOD
			std::cerr << std::fixed << std::setprecision(4) ;
			int const N = m_design->outcome().size() ;
			std::cerr << "==== NoPredictorLogisticRegressionLogLikelihood::evaluate_at() ====\n" ;
			std::cerr << "Params: " << parameters.transpose() << "\n" ;
			std::cerr << "Number of samples: " << N << ".\n" ;
			std::cerr << "Included samples:" ;
			for( std::size_t i = 0; i < included_samples.size(); ++i ) {
				std::cerr << " " << included_samples[i].begin() << "-" << included_samples[i].end() ;
			}
			std::cerr << "\n" ;
			std::cerr << "Parameters: " << parameters << "\n" ;
			std::cerr << "Phenotypes:\n" << m_design->outcome().head( std::min( N, 5 ) ) << "...\n";
			std::cerr << "Outcome probabilities:\n"
				<< m_outcome_probabilities.topRows( std::min( N, 5 ) )
				<< "...\n" ;
#endif
			// Calculate log-likelihood.
			// We sum over samples ;
			compute_value_of_function( included_samples ) ;

#if DEBUG_LOGLIKELIHOOD
			std::cerr << "function value = \n"
				<< m_value_of_function << "\n" ;
#endif
			if( numberOfDerivatives > 0 ) {
				compute_value_of_first_derivative( included_samples ) ;
#if DEBUG_LOGLIKELIHOOD
				std::cerr << "derivative = \n"
					<< m_value_of_first_derivative << "\n" ;
#endif
				if( numberOfDerivatives > 1 ) {
					compute_value_of_second_derivative( included_samples ) ;
#if DEBUG_LOGLIKELIHOOD
					std::cerr << "2nd derivative = \n"
						<< m_value_of_second_derivative << ".\n" ;
#endif
				}
			}
		}
		
		NoPredictorLogisticRegressionLogLikelihood::Point const& NoPredictorLogisticRegressionLogLikelihood::get_parameters() const {
			return m_parameters ;
		}

		// Calculate P( outcome | genotype, covariates, parameters ).
		void NoPredictorLogisticRegressionLogLikelihood::calculate_outcome_probabilities(
			Vector const& parameters,
			Vector const& outcomes,
			NoPredictorLogisticRegressionLogLikelihood::Vector* result
		) const {
			assert( result != 0 ) ;
			result->resize( outcomes.size() ) ;
			Matrix const& design_matrix = m_design->matrix() ;
			evaluate_mean_function( design_matrix * parameters, outcomes, result ) ;
		}

		void NoPredictorLogisticRegressionLogLikelihood::evaluate_mean_function(
			Vector const& linear_combinations,
			Vector const& outcomes,
			Vector* result
		) const {
			assert( result != 0 ) ;
			assert( linear_combinations.size() == outcomes.size() ) ;
			Vector const ones = Vector::Ones( linear_combinations.size() ) ;
			Vector const exps = linear_combinations.array().exp() ;
			(*result) = ( ones.array() + outcomes.array() * ( exps.array() - ones.array() ) )  / ( ones.array() + exps.array() ) ;
		}

		void NoPredictorLogisticRegressionLogLikelihood::compute_value_of_function( std::vector< metro::SampleRange > const& included_samples ) {
			m_value_of_function = 0.0 ;
			m_V = m_design->sample_weights().array() * m_outcome_probabilities.array() ;
			for( std::size_t i = 0; i < included_samples.size(); ++i ) {
				int const start_row = included_samples[i].begin() ;
				int const end_row = included_samples[i].end() ;
				m_value_of_function += (
					m_V.segment( start_row, end_row - start_row ).array().log().sum()
				) ;
			}
		}

		void NoPredictorLogisticRegressionLogLikelihood::compute_value_of_first_derivative( std::vector< metro::SampleRange > const& included_samples ) {
			int const N = m_design->outcome().size() ;
			int const D = m_design->matrix().cols() ;
			Matrix const& design_matrix = m_design->matrix() ;
			
			m_V = Vector::Ones( N ) - m_outcome_probabilities ;

			m_value_of_first_derivative = Vector::Zero( D ) ;
			Vector signs = (( m_design->outcome() * 2.0 ) - Vector::Ones( N ) ) ; // == (-1)^{phenotype + 1}
			// multiply ith row of design matrix by sign times ith entry of column of V.
			for( std::size_t i = 0; i < included_samples.size(); ++i ) {
				int const start_row = included_samples[i].begin() ;
				int const end_row = included_samples[i].end() ;
				//std::cerr << "i = " << i << ", V.rows() = " << V.rows() << ", start_row = " << start_row << ", end_row = " << end_row << ".\n" ;
				m_value_of_first_derivative += (
					(
						signs.segment( start_row, end_row - start_row ).array()
						* m_V.segment( start_row, end_row - start_row ).array()
					).matrix()
					.asDiagonal()
					* design_matrix.block( start_row, 0, end_row - start_row, design_matrix.cols () )
				).colwise().sum() ;
			}
		}

		void NoPredictorLogisticRegressionLogLikelihood::compute_value_of_second_derivative( std::vector< metro::SampleRange > const& included_samples ) {
			int const D = m_design->matrix().cols() ;
			int const N = m_design->matrix().rows() ;
			
			m_value_of_second_derivative = Matrix::Zero( D, D ) ;
			m_V = m_V.array() * ( ( Vector::Ones( N ) - 2.0 * m_outcome_probabilities ).array() - m_V.array() );

			{
				for( std::size_t i = 0; i < included_samples.size(); ++i ) {
					int const start_row = included_samples[i].begin() ;
					int const end_row = included_samples[i].end() ;
					for( int sample = start_row; sample < end_row; ++sample ) {
						m_value_of_second_derivative += m_V( sample ) * m_design_matrix_row_tensor_squares.block( D * sample, 0, D, D ) ;
					}
				}
			}
		}
	

		double NoPredictorLogisticRegressionLogLikelihood::get_value_of_function() const {
			return m_value_of_function ;
		}

		NoPredictorLogisticRegressionLogLikelihood::Vector NoPredictorLogisticRegressionLogLikelihood::get_value_of_first_derivative() const {
			return m_value_of_first_derivative ;
		}

		NoPredictorLogisticRegressionLogLikelihood::Matrix NoPredictorLogisticRegressionLogLikelihood::get_value_of_second_derivative() const {
			return m_value_of_second_derivative ;
		}
	}
}
