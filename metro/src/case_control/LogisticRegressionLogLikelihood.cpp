
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include <vector>
#include <utility>
#include <algorithm>
#include <boost/format.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include "metro/case_control/LogisticRegressionLogLikelihood.hpp"
#include "metro/intersect_ranges.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/Error.hpp"

// #define DEBUG_LOGLIKELIHOOD 2

namespace metro {
	namespace case_control {
		LogisticRegressionLogLikelihood::UniquePtr LogisticRegressionLogLikelihood::create(
			RegressionDesign& design
		) {
			return LogisticRegressionLogLikelihood::UniquePtr(
				new LogisticRegressionLogLikelihood( design )
			) ;
		}
		
		LogisticRegressionLogLikelihood::LogisticRegressionLogLikelihood(
			RegressionDesign& design
		):
			m_design( &design ),
			m_parameters( Vector::Zero( m_design->matrix().cols() ) ),
			m_state( e_Uncomputed )
		{
			assert( m_design != 0 ) ;
			// Verify design looks sane.
			Eigen::VectorXd const& outcome = m_design->outcome() ;
			for( int i = 0; i < outcome.size(); ++i ) {
				if( outcome(i) == outcome(i) && outcome(i) != 0 && outcome(i) != 1 ) {
					throw genfile::BadArgumentError(
						"metro::case_control::LogisticRegressionLogLikelihood::LogisticRegressionLogLikelihood()",
						"design",
						( boost::format( "Outcome (%f for individual %d) has non-missing level other than zero or one." ) % outcome(i) % (i+1) ).str()
					) ;
				}
			}
		}

		std::string LogisticRegressionLogLikelihood::get_summary() const {
			using genfile::string_utils::to_string ;
			return "LogisticRegressionLogLikelihood( "
				+ to_string( m_design->outcome().size() )
				+ " samples )" ;
		}

		void LogisticRegressionLogLikelihood::set_parameter_naming_scheme( LogisticRegressionLogLikelihood::GetParameterName get_parameter_name ) {
			assert( get_parameter_name ) ;
			m_get_parameter_name = get_parameter_name ;		
		}
	
		std::string LogisticRegressionLogLikelihood::get_parameter_name( std::size_t i ) const {
			LogLikelihood::IntegerMatrix const identity = identify_parameters() ;
			return m_get_parameter_name( m_design->get_predictor_name( identity(i,1)), identity(i,0) ) ;
		}
	
		LogLikelihood::IntegerMatrix LogisticRegressionLogLikelihood::identify_parameters() const {
			LogLikelihood::IntegerMatrix result = LogLikelihood::IntegerMatrix::Zero( m_design->matrix().cols(), 2 ) ;
			// Each parameter corresponds to an effect for the non-baseline outcome level
			// and a column of the design matrix.
			result.col(0).setConstant( 1.0 ) ;
			// There is one parameter per design matrix column.
			for( int i = 0; i < m_design->matrix().cols(); ++i ) {
				result( i, 1 ) = i ;
			}
			return result ;
		}

		int LogisticRegressionLogLikelihood::number_of_outcomes() const {
			return 2 ;
		}

		void LogisticRegressionLogLikelihood::set_predictor_levels(
			Matrix const& levels,
			Matrix const& probabilities,
			std::vector< metro::SampleRange > const& included_samples
		) {
			m_design->set_predictor_levels( levels, probabilities, included_samples ) ;
			m_state = e_Uncomputed ;
#if DEBUG_LOGLIKELIHOOD
			std::cerr << "LogisticRegressionLogLikelihood::set_predictor_probs(): Included samples:" ;
			for( std::size_t i = 0; i < included_samples.size(); ++i ) {
				std::cerr << " " << included_samples[i].begin() << "-" << included_samples[i].end() ;
			}
			std::cerr << "\n" ;
#endif
		}

		void LogisticRegressionLogLikelihood::evaluate_at( Point const& parameters, int const numberOfDerivatives ) {
			evaluate_at_impl( parameters, m_design->per_predictor_included_samples(), numberOfDerivatives ) ;
		}

		void LogisticRegressionLogLikelihood::evaluate_at_impl(
			Point const& parameters,
			std::vector< metro::SampleRange > const& included_samples,
			int const numberOfDerivatives
		) {
//			m_state = eUncomputed ;
			assert( parameters.size() == m_design->matrix().cols() ) ;
			assert( numberOfDerivatives < 3 ) ;

			if( m_parameters != parameters ) {
				m_state = e_Uncomputed ;
				m_parameters = parameters ;
			}

#if DEBUG_LOGLIKELIHOOD
			std::cerr << std::fixed << std::setprecision(4) ;
			int const N = m_design->outcome().size() ;
			std::cerr << "==== LogisticRegressionLogLikelihood::evaluate_at() ====\n" ;
			std::cerr << "  m_state: " << ( boost::format( "0x%x" ) % m_state ).str() << ", "
				<< "Number of samples: " << N << ", "
				<< "Included samples:" ;
			for( std::size_t i = 0; i < included_samples.size(); ++i ) {
				std::cerr << " " << included_samples[i].begin() << "-" << included_samples[i].end() ;
			}
			std::cerr << "\n" ;
			std::cerr << "  Params: " << parameters.transpose() << ".\n" ;
			std::cerr << "  Phenotypes: " << m_design->outcome().head( std::min( N, 5 ) ).transpose() << "...\n";
			//std::cerr << "  Predictor levels: " << m_design->get_predictor_levels().transpose() << ".\n" ;
#endif

			bool computeFunction = ( m_state == e_Uncomputed ) ;
			bool compute1stDerivative = ( (numberOfDerivatives > 0) & !(m_state & e_Computed1stDerivative )) ;
			bool compute2ndDerivative = ( (numberOfDerivatives > 1) & !(m_state & e_Computed2ndDerivative )) ;

			// Calculate log-likelihood if needed.
			if( computeFunction ) {
	 			calculate_outcome_probabilities( parameters, m_design->outcome(), &m_outcome_probabilities ) ;
				m_A = m_design->get_predictor_level_probabilities().array() * m_outcome_probabilities.array() ;
				compute_value_of_function( m_A, included_samples ) ;
				m_state |= e_ComputedFunction ;
				// (m_state | e_ComputedFunction) = 1 implies m_A and m_outcome_probabilities, as well as function value, are computed
			}

#if DEBUG_LOGLIKELIHOOD
			std::cerr << "  " << ( computeFunction ? "(computed) " : "(memoized) " ) << "function value = " << m_value_of_function << "\n" ;
#endif

			if( compute1stDerivative ) {
				compute_value_of_first_derivative( m_A, included_samples, &m_B ) ;
				m_state |= e_Computed1stDerivative ;
				// m_state | e_Computed1stDerivative implies m_B, as well as 1st derivative, are computed
			}
#if DEBUG_LOGLIKELIHOOD
			if( m_state & e_Computed1stDerivative ) {
				std::cerr << "  " << ( compute1stDerivative ? "(computed) " : "(memoized) " ) << "derivative = " << m_value_of_first_derivative.transpose() << "\n" ;
			} else {
				std::cerr << "  (1st derivative not computed).\n" ;
			}
#endif

			if( compute2ndDerivative ) {
				compute_value_of_second_derivative( m_B, included_samples ) ;
				m_state |= e_Computed2ndDerivative ;
			}
#if DEBUG_LOGLIKELIHOOD
			if( m_state & e_Computed2ndDerivative ) {
				std::cerr << "  " << ( compute2ndDerivative ? "(computed) " : "(memoized) " ) << "2nd derivative = \n"
					<< m_value_of_second_derivative << ".\n" ;
			} else {
				std::cerr << "  (2nd derivative not computed).\n" ;
			}
			std::cerr << "====\n" ;
#endif
		}
		
		LogisticRegressionLogLikelihood::Point const& LogisticRegressionLogLikelihood::get_parameters() const {
			return m_parameters ;
		}

		// Calculate P( outcome | genotype, covariates, parameters ).
		void LogisticRegressionLogLikelihood::calculate_outcome_probabilities(
			Vector const& parameters,
			Vector const& outcome,
			LogisticRegressionLogLikelihood::Matrix* result
		) const {
			result->resize( outcome.size(), m_design->get_number_of_predictor_levels() ) ;
			Vector linear_combination ;
			int const number_of_predictor_levels = m_design->get_number_of_predictor_levels() ;
			for( int g = 0; g < number_of_predictor_levels; ++g ) {
				Matrix const& design_matrix = m_design->set_predictor_level( g ).matrix() ;
				result->col( g ) = evaluate_mean_function( design_matrix * parameters, outcome ) ;
			}
		}
		
		LogisticRegressionLogLikelihood::Vector LogisticRegressionLogLikelihood::evaluate_mean_function( Vector const& linear_combinations, Vector const& outcomes ) const {
			assert( linear_combinations.size() == outcomes.size() ) ;
			Vector exps = linear_combinations.array().exp() ;
			Vector ones = Vector::Ones( linear_combinations.size() ) ;
			return ( ones.array() + outcomes.array() * ( exps.array() - ones.array() ) )  / ( ones + exps ).array() ;
		}

		void LogisticRegressionLogLikelihood::compute_value_of_function( Matrix const& V, std::vector< metro::SampleRange > const& included_samples ) {
			m_value_of_function = 0.0 ;
			for( std::size_t i = 0; i < included_samples.size(); ++i ) {
				int const start_row = included_samples[i].begin() ;
				int const end_row = included_samples[i].end() ;
				m_value_of_function += (
					V.block( start_row, 0, end_row - start_row, V.cols() ).rowwise().sum().array().log().sum()
				) ;
			}
		}

		void LogisticRegressionLogLikelihood::compute_value_of_first_derivative( Matrix const& A, std::vector< metro::SampleRange > const& included_samples, Matrix* B ) {
			int const N = m_design->outcome().size() ;
			int const D = m_design->matrix().cols() ;
			int const N_predictor_levels = m_design->get_number_of_predictor_levels() ;
			// ...and its first derivative...
			m_value_of_first_derivative = Vector::Zero( D ) ;
			(*B) = A ;
			{
				RowVector ones( RowVector::Ones( N_predictor_levels )) ;
				for( int i = 0; i < B->rows(); ++i ) {
					B->row(i) /= B->row(i).sum() ;
					B->row(i).array() *= ( ones - m_outcome_probabilities.row(i) ).array() ;
				}
			}

			{
				m_first_derivative_terms.setZero( m_design->matrix().rows(), m_design->matrix().cols() ) ;
				Vector signs = (( m_design->outcome() * 2.0 ) - Vector::Ones( N ) ) ; // == (-1)^{phenotype + 1}
				for( int g = 0; g < N_predictor_levels; ++g ) {
					Matrix const& design_matrix = m_design->set_predictor_level( g ).matrix() ;
#if DEBUG_LOGLIKELIHOOD
					std::cerr << "LogisticRegressionLogLikelihood::compute_value_of_first_derivative(): design_matrix for predictor level " << g << " =\n"
						<< design_matrix.topRows( std::min( N, 5 ) ) << "\n" ;
#endif
					//m_design->matrix().block( 0, 1, m_design->matrix().rows(), predictor_levels.cols() ).rowwise() = predictor_levels.row( g ) ;
					for( std::size_t i = 0; i < included_samples.size(); ++i ) {
						int const start_row = included_samples[i].begin() ;
						int const end_row = included_samples[i].end() ;
						MatrixBlock block = m_first_derivative_terms.block( start_row, 0, end_row - start_row, m_first_derivative_terms.cols() ) ;
						// multiply ith row of design matrix by sign times ith entry of column of B.
						block += (
								signs.segment( start_row, end_row - start_row ).array()
								* B->col( g ).segment( start_row, end_row - start_row ).array()
							).matrix().asDiagonal()
							* design_matrix.block( start_row, 0, end_row - start_row, design_matrix.cols() )
						;
					}
				}
				m_value_of_first_derivative = m_first_derivative_terms.colwise().sum() ;
			}
		}

		void LogisticRegressionLogLikelihood::compute_value_of_second_derivative( Matrix const& B, std::vector< metro::SampleRange > const& included_samples ) {
			int const D = m_design->matrix().cols() ;
			int const N_predictor_levels = m_design->get_number_of_predictor_levels() ;
			m_value_of_second_derivative = Matrix::Zero( D, D ) ;
			// Compute second derivative...
			for( int row_i = 0; row_i < D; ++row_i ) {
				MatrixBlock resultRow = m_value_of_second_derivative.block(row_i, 0, 1, D ) ;
				for( int g1 = 0; g1 < N_predictor_levels; ++g1 ) {
					Matrix const& design_matrix = m_design->set_predictor_level( g1 ).matrix() ;
					for( std::size_t sample_i = 0; sample_i < included_samples.size(); ++sample_i ) {
						int const start_row = included_samples[sample_i].begin() ;
						int const end_row = included_samples[sample_i].end() ;
						ConstMatrixBlock design_matrix_block = design_matrix.block( start_row, 0, end_row - start_row, design_matrix.cols() ) ;
						resultRow += (
								B.col(g1).segment( start_row, end_row - start_row ).array()
								* (
									Vector::Constant( end_row - start_row, 1.0 )
									- 2.0 * m_outcome_probabilities.col( g1 ).segment( start_row, end_row - start_row )
								).array()
								* design_matrix_block.col( row_i ).array()
							).matrix().transpose()
							* design_matrix_block
						;
					}
				}
				resultRow -= ( m_first_derivative_terms.col(row_i).asDiagonal() * m_first_derivative_terms ).colwise().sum() ;
			}
		}
		
		double LogisticRegressionLogLikelihood::get_value_of_function() const {
			return m_value_of_function ;
		}

		LogisticRegressionLogLikelihood::Vector LogisticRegressionLogLikelihood::get_value_of_first_derivative() const {
			return m_value_of_first_derivative ;
		}

		LogisticRegressionLogLikelihood::Matrix LogisticRegressionLogLikelihood::get_value_of_second_derivative() const {
			return m_value_of_second_derivative ;
		}
	}
}
