#include <iostream>
#include <vector>
#include <boost/iterator/counting_iterator.hpp>
#include "snptest/case_control/NullModelLogLikelihood.hpp"

namespace snptest {
	namespace case_control {
		NullModelLogLikelihood::NullModelLogLikelihood(
			Vector const& phenotypes,
			FinitelySupportedFunctionSet const& genotypes,
			Matrix const& covariates
		):
			m_phenotypes( phenotypes ),
			m_genotype_call_probabilities( genotypes.get_values() ),
			m_covariates( covariates ),
			m_design_matrix( calculate_design_matrix( covariates ) )
		{
			assert( m_covariates.cols() == 0 || m_covariates.rows() == m_phenotypes.size() ) ;
		}

		NullModelLogLikelihood::NullModelLogLikelihood(
			Vector const& phenotypes,
			FinitelySupportedFunctionSet const& genotypes,
			Matrix const& covariates,
			std::vector< std::size_t > const& excluded_samples
		):
			m_phenotypes( phenotypes ),
			m_genotype_call_probabilities( genotypes.get_values() ),
			m_covariates( covariates ),
			m_design_matrix( calculate_design_matrix( covariates ) )
		{
			assert( m_phenotypes.size() == m_genotype_call_probabilities.rows() ) ;
			assert( excluded_samples.size() <= std::size_t( m_phenotypes.rows() ) ) ;
			if( excluded_samples.size() > 0 ) {
				deal_with_exclusions( excluded_samples ) ;
			}
		}

		NullModelLogLikelihood::Matrix NullModelLogLikelihood::calculate_design_matrix( Matrix const& covariates ) const {
			Matrix result( m_phenotypes.rows(), covariates.cols() + 1 ) ;
			result.leftCols( 1 ).setOnes() ;
			if( covariates.cols() > 0 ) {
				result.rightCols( covariates.cols() ) = covariates ;
			}
			return result ;
		}

		void NullModelLogLikelihood::deal_with_exclusions( std::vector< std::size_t > exclusions ) {
			assert( exclusions.size() > 0 ) ;
			std::sort( exclusions.begin(), exclusions.end() ) ;
			assert( exclusions.back() < std::size_t( m_phenotypes.size() ) ) ;
			exclusions.push_back( m_phenotypes.size() ) ;
			Matrix remove_exclusions = Matrix::Zero( m_phenotypes.rows() - exclusions.size(), m_phenotypes.rows() ) ;
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

		void NullModelLogLikelihood::evaluate_at( Point const& parameters ) {
			int const N = m_phenotypes.size() ;
			int const D = m_design_matrix.cols() ;
			assert( parameters.size() == D ) ;

			m_outcome_probabilities = calculate_outcome_probabilities( parameters, m_phenotypes, m_design_matrix ) ;
			assert( m_outcome_probabilities.rows() == m_phenotypes.size() && m_outcome_probabilities.cols() == 1 ) ;

#if DEBUG_LOGLIKELIHOOD			
			std::cerr << "==== NullModelLogLikelihood::evaluate_at() ====\n" ;
			std::cerr << "Parameters: " << parameters << "\n" ;
			std::cerr << "Phenotypes:\n" << m_phenotypes.head( std::min( N, 5 ) ) << "...\n";
			std::cerr << "Outcome probabilities:\n"
				<< m_outcome_probabilities.topRows( std::min( N, 5 ) )
				<< "...\n" ;
#endif
			
			// Calculate log-likelihood...
			Vector V = ( m_genotype_call_probabilities.rowwise().sum().array() * m_outcome_probabilities.array() ) ;
			m_value_of_function = V.rowwise().sum().array().log().sum() ;

#if DEBUG_LOGLIKELIHOOD			
			std::cerr << "call prob sums:\n"
				<< m_genotype_call_probabilities.rowwise().sum().head( std::min( N, 5 ) ) << "...\n" ;

			std::cerr << "V = call prob sums * outcome probs:\n"
				<< V.head( std::min( N, 5 ) )
				<< "...\n"
				<< "Function value = " << m_value_of_function << ".\n" ;
#endif
			// ...and its first and second derivatives...
			V = Vector::Ones( N ) - m_outcome_probabilities ;
			
#if DEBUG_LOGLIKELIHOOD
			std::cerr << "V = ( 1 - outcome_probs ):\n"
				<< V.topRows( std::min( N, 5 ) )
				<< "...\n" ;
#endif

			{
				Vector signs = (( m_phenotypes * 2.0 ) - Vector::Ones( N ) ) ; // == (-1)^{phenotype + 1}
				// multiply ith row of design matrix by sign times ith entry of column of V.
				RowVector second_term = ( V.asDiagonal() * m_design_matrix ).colwise().sum() ;
				m_value_of_first_derivative = (
					( signs.array() * V.array() ).matrix().asDiagonal() * m_design_matrix
				).colwise().sum() ;
			}
			
			// ...and its second derivative...

			m_value_of_second_derivative = Matrix::Zero( m_design_matrix.cols(), m_design_matrix.cols() ) ;
			Vector coeffs = V.array() * ( Vector::Ones( N ) - 2.0 * m_outcome_probabilities ).array() ;
			for( int sample = 0; sample < N; ++sample ) {
				m_value_of_second_derivative +=
					coeffs( sample ) * ( m_design_matrix.row( sample ).transpose() * m_design_matrix.row( sample ) ) ;
				RowVector second_term = V( sample ) * m_design_matrix.row( sample ) ;
				m_value_of_second_derivative -= second_term.transpose() * second_term ;
			}
#if DEBUG_LOGLIKELIHOOD
			std::cerr << "first derivative:\n" << m_value_of_first_derivative << ".\n" ;
			std::cerr << "second derivative:\n" << m_value_of_second_derivative << ".\n" ;
#endif			
		}
		
		// Calculate P( outcome | genotype, covariates, parameters ).
		NullModelLogLikelihood::Matrix NullModelLogLikelihood::calculate_outcome_probabilities(
			Vector const& parameters,
			Vector const& phenotypes,
			Matrix& design_matrix
		) const {
			return evaluate_mean_function( design_matrix * parameters, phenotypes ) ;
		}

		NullModelLogLikelihood::Vector NullModelLogLikelihood::evaluate_mean_function( Vector const& linear_combinations, Vector const& outcomes ) const {
			assert( linear_combinations.size() == outcomes.size() ) ;
			Vector exps = linear_combinations.array().exp() ;
			Vector ones = Vector::Ones( linear_combinations.size() ) ;
			return ( ones.array() + outcomes.array() * ( exps.array() - ones.array() ) )  / ( ones + exps ).array() ;
		}

		double NullModelLogLikelihood::get_value_of_function() const {
			return m_value_of_function ;
		}

		NullModelLogLikelihood::Vector NullModelLogLikelihood::get_value_of_first_derivative() const {
			return m_value_of_first_derivative ;
		}

		NullModelLogLikelihood::Matrix NullModelLogLikelihood::get_value_of_second_derivative() const {
			return m_value_of_second_derivative ;
		}
	}
}
