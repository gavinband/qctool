#include <iostream>
#include <vector>
#include <utility>
#include <boost/iterator/counting_iterator.hpp>
#include "snptest/case_control/LogLikelihood.hpp"

// #define DEBUG_LOGLIKELIHOOD 1

namespace snptest {
	namespace case_control {
		LogLikelihood::LogLikelihood(
			Vector const& phenotypes,
			FinitelySupportedFunctionSet const& genotypes,
			Matrix const& covariates
		):
			m_phenotypes( phenotypes ),
			m_genotype_call_probabilities( genotypes.get_values() ),
			m_genotype_levels( genotypes.get_support() ),
			m_covariates( covariates ),
			m_design_matrix( calculate_design_matrix( covariates ) )
		{
			assert( m_phenotypes.size() == m_genotype_call_probabilities.rows() ) ;
			assert( m_covariates.cols() == 0 || m_covariates.rows() == m_phenotypes.size() ) ;
			compute_tensor_squares( m_design_matrix ) ;
		}

		LogLikelihood::LogLikelihood(
			Vector const& phenotypes,
			FinitelySupportedFunctionSet const& genotypes,
			Matrix const& covariates,
			std::vector< std::size_t > const& excluded_samples
		):
			m_phenotypes( phenotypes ),
			m_genotype_call_probabilities( genotypes.get_values() ),
			m_genotype_levels( genotypes.get_support() ),
			m_covariates( covariates ),
			m_design_matrix( calculate_design_matrix( covariates ) )
		{
			assert( m_phenotypes.size() == m_genotype_call_probabilities.rows() ) ;
			assert( excluded_samples.size() <= std::size_t( m_phenotypes.size() ) ) ;
			if( excluded_samples.size() > 0 ) {
				deal_with_exclusions( excluded_samples ) ;
			}
			compute_tensor_squares( m_design_matrix ) ;
		}

		LogLikelihood::Matrix LogLikelihood::calculate_design_matrix( Matrix const& covariates ) const {
			Matrix result( m_phenotypes.rows(), covariates.cols() + 2 ) ;
			result.leftCols( 1 ).setOnes() ;
			if( covariates.cols() > 0 ) {
				result.rightCols( covariates.cols() ) = covariates ;
			}
			return result ;
		}
		
		// Construct a big matrix of size (D*N) x (D*3).
		// The block at position (D*sample,D*g) is the tensor square of the design matrix
		// row for the given sample and given genotype.
		void LogLikelihood::compute_tensor_squares( Matrix& design_matrix ) {
			int const D = design_matrix.cols() ;
			m_design_matrix_row_tensor_squares.resize( D * m_phenotypes.size(), D * 3 ) ;
			for( std::size_t g = 0; g < 3; ++g ) {
				design_matrix.col(1).setConstant( m_genotype_levels( g ) ) ;
				for( int sample = 0; sample < design_matrix.rows(); ++sample ) {
					if( sample > 0 && design_matrix.row( sample ) == design_matrix.row( sample - 1 )) {
						m_design_matrix_row_tensor_squares.block( D * sample, D * g, D, D ) = m_design_matrix_row_tensor_squares.block( D * ( sample - 1 ), D * g, D, D ) ;
					}
					else {
						m_design_matrix_row_tensor_squares.block( D * sample, D * g, D, D ) = design_matrix.row( sample ).transpose() * design_matrix.row( sample ) ;
					}
				}
			}
		}
		
		void LogLikelihood::deal_with_exclusions( std::vector< std::size_t > exclusions ) {
			assert( exclusions.size() > 0 ) ;
			std::sort( exclusions.begin(), exclusions.end() ) ;
			assert( exclusions.back() < std::size_t( m_phenotypes.size() ) ) ;
			Matrix remove_exclusions = Matrix::Zero( m_phenotypes.rows() - exclusions.size(), m_phenotypes.rows() ) ;
			exclusions.push_back( m_phenotypes.size() ) ;
			int	block_start_column = 0 ;
			for( std::size_t i = 0; i < exclusions.size(); ++i ) {
				int block_start_row = block_start_column - i ;
				int block_end_column = exclusions[i] ;
				int block_size = block_end_column - block_start_column ;
				remove_exclusions.block( block_start_row, block_start_column, block_size, block_size )
					= Matrix::Identity( block_size, block_size ) ;
				block_start_column = block_end_column + 1 ;
			}
			
			m_design_matrix = remove_exclusions * m_design_matrix ;
			m_phenotypes = remove_exclusions * m_phenotypes ;
			m_genotype_call_probabilities = remove_exclusions * m_genotype_call_probabilities ;
		}

		void LogLikelihood::evaluate_at( Point const& parameters ) {
			calculate_outcome_probabilities( parameters, m_phenotypes, m_design_matrix, &m_outcome_probabilities ) ;
			assert( m_outcome_probabilities.rows() == m_genotype_call_probabilities.rows() && m_outcome_probabilities.cols() == m_genotype_call_probabilities.cols() ) ;

#if DEBUG_LOGLIKELIHOOD
			std::cerr << "==== LogLikelihood::evaluate_at() ====\n" ;
			std::cerr << "Number of samples: " << N << ".\n" ;
			std::cerr << "Parameters: " << parameters << "\n" ;
			std::cerr << "Phenotypes:\n" << m_phenotypes.head( std::min( N, 5 ) ) << "...\n";
			std::cerr << "Outcome probabilities:\n"
				<< m_outcome_probabilities.topRows( std::min( N, 5 ) )
				<< "...\n" ;
#endif
			// Calculate log-likelihood.
			// We sum over samples ignoring missing samples.
			m_V = ( m_genotype_call_probabilities.array() * m_outcome_probabilities.array() ) ;
#if DEBUG_LOGLIKELIHOOD
			std::cerr << "call probs:\n"
				<< m_genotype_call_probabilities.topRows( std::min( N, 5 ) ) << "...\n" ;
			std::cerr << "call values:\n"
				<< m_genotype_levels << ".\n" ;
			std::cerr << "V = call probs * outcome probs:\n"
				<< m_V.topRows( std::min( N, 5 ) )
				<< "...\n" ;
#endif
			compute_value_of_function( m_V ) ;
			compute_value_of_first_derivative( m_V ) ;
			compute_value_of_second_derivative( m_V ) ;
		}

		void LogLikelihood::compute_value_of_function( Matrix const& V ) {
			m_value_of_function = V.rowwise().sum().array().log().sum() ;
		#if DEBUG_LOGLIKELIHOOD
			std::cerr << "Function value = " << m_value_of_function << ".\n" ;
		#endif
		}

		void LogLikelihood::compute_value_of_first_derivative( Matrix& V ) {
			int const N = m_phenotypes.size() ;
			int const D = m_design_matrix.cols() ;
			
			// ...and its first derivative...
			{
				RowVector ones( RowVector::Ones( 3 )) ;
				for( int i = 0; i < V.rows(); ++i ) {
					V.row(i) /= V.row(i).sum() ;
					V.row(i).array() *= ( ones - m_outcome_probabilities.row(i) ).array() ;
				}
			}
		#if DEBUG_LOGLIKELIHOOD
			std::cerr << "V = ( call probs * outcome probs * ( 1 - outcome_probs )) / sum( call probs * outcome probs ):\n"
				<< V.topRows( std::min( N, 5 ) )
				<< "...\n" ;
		#endif

			m_value_of_first_derivative = Vector::Zero( D ) ;
			{
				Vector signs = (( m_phenotypes * 2.0 ) - Vector::Ones( N ) ) ; // == (-1)^{phenotype + 1}
				for( std::size_t g = 0; g < 3; ++g ) {
					m_design_matrix.col(1).setConstant( m_genotype_levels( g ) ) ;
					// multiply ith row of design matrix by sign times ith entry of column of V.
					m_value_of_first_derivative += (
						( signs.array() * V.col( g ).array() ).matrix().asDiagonal() * m_design_matrix
					).colwise().sum() ;
				}
			}
		}

		void LogLikelihood::compute_value_of_second_derivative( Matrix const& V ) {
			int const N = m_phenotypes.size() ;
			int const D = m_design_matrix.cols() ;

			// ...and its second derivative...
			m_value_of_second_derivative = Matrix::Zero( m_design_matrix.cols(), m_design_matrix.cols() ) ;
			Vector coeffs( N ) ;
			for( std::size_t g = 0; g < 3; ++g ) {
				coeffs = V.col( g ).array() * ( Vector::Ones( N ) - 2.0 * m_outcome_probabilities.col( g ) ).array() ;
				m_design_matrix.col(1).setConstant( m_genotype_levels( g ) ) ;
				for( int sample = 0; sample < N; ++sample ) {
					double const coeff = coeffs( sample ) - V( sample, g ) * V( sample, g ) ;
					m_value_of_second_derivative += coeff * m_design_matrix_row_tensor_squares.block( D * sample, D * g, D, D ) ;
				}
			}
			#if DEBUG_LOGLIKELIHOOD
				std::cerr << "second derivative:\n" << m_value_of_second_derivative << ".\n" ;
			#endif
		}
	
		// Calculate P( outcome | genotype, covariates, parameters ).
		void LogLikelihood::calculate_outcome_probabilities(
			Vector const& parameters,
			Vector const& phenotypes,
			Matrix& design_matrix,
			LogLikelihood::Matrix* result
		) const {
			result->resize( phenotypes.size(), 3 ) ;
			Vector linear_combination ;
			for( std::size_t g = 0; g < 3; ++g ) {
				design_matrix.col( 1 ).setConstant( m_genotype_levels( g ) ) ;
				result->col( g ) = evaluate_mean_function( design_matrix * parameters, phenotypes ) ;
			}
		}

		LogLikelihood::Vector LogLikelihood::evaluate_mean_function( Vector const& linear_combinations, Vector const& outcomes ) const {
			assert( linear_combinations.size() == outcomes.size() ) ;
			Vector exps = linear_combinations.array().exp() ;
			Vector ones = Vector::Ones( linear_combinations.size() ) ;
			return ( ones.array() + outcomes.array() * ( exps.array() - ones.array() ) )  / ( ones + exps ).array() ;
		}

		double LogLikelihood::get_value_of_function() const {
			return m_value_of_function ;
		}

		LogLikelihood::Vector LogLikelihood::get_value_of_first_derivative() const {
			return m_value_of_first_derivative ;
		}

		LogLikelihood::Matrix LogLikelihood::get_value_of_second_derivative() const {
			return m_value_of_second_derivative ;
		}
	}
}
