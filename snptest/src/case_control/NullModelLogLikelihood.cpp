#include <iostream>
#include <set>
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
			deal_with_exclusions( std::vector< std::size_t >() ) ;
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
			deal_with_exclusions( excluded_samples ) ;
		}

		NullModelLogLikelihood::Matrix NullModelLogLikelihood::calculate_design_matrix( Matrix const& covariates ) const {
			Matrix result( m_phenotypes.rows(), covariates.cols() + 1 ) ;
			result.leftCols( 1 ).setOnes() ;
			if( covariates.cols() > 0 ) {
				result.rightCols( covariates.cols() ) = covariates ;
			}
			return result ;
		}

		void NullModelLogLikelihood::deal_with_exclusions( std::vector< std::size_t > const& exclusions ) {
			std::set< std::size_t > sorted_exclusions( exclusions.begin(), exclusions.end() ) ;
			sorted_exclusions.insert( m_phenotypes.size() ) ;
			m_exclusions.assign( sorted_exclusions.begin(), sorted_exclusions.end() ) ;
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
			{
				m_value_of_function = 0.0 ;
				int start_row = 0 ;
				for( std::size_t i = 0; i < m_exclusions.size(); ++i ) {
					int const end_row = m_exclusions[i] ;
					m_value_of_function += V.segment( start_row, end_row - start_row ).rowwise().sum().array().log().sum() ;
					start_row = end_row + 1 ;
				}
			}

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
				m_value_of_first_derivative = Vector::Zero( D ) ;
				int start_row = 0 ;
				for( std::size_t i = 0; i < m_exclusions.size(); ++i ) {
					int const end_row = m_exclusions[i] ;
					std::cerr << "i = " << i << ", V.rows() = " << V.rows() << ", start_row = " << start_row << ", end_row = " << end_row << ".\n" ;
					m_value_of_first_derivative += (
						( signs.segment( start_row, end_row - start_row ).array() * V.segment( start_row, end_row - start_row ).array() ).matrix().asDiagonal()
						* m_design_matrix.block( start_row, 0, end_row - start_row, m_design_matrix.cols () )
					).colwise().sum() ;
					start_row = end_row + 1 ;
				}
			}
			
			// ...and its second derivative...

			m_value_of_second_derivative = Matrix::Zero( m_design_matrix.cols(), m_design_matrix.cols() ) ;
			Vector coeffs = V.array() * ( Vector::Ones( N ) - 2.0 * m_outcome_probabilities ).array() ;
			{
				int sample = 0 ;
				for( int i = 0; i < m_exclusions.size(); ++i ) {
					for( ; sample < m_exclusions[i]; ++sample ) {
						m_value_of_second_derivative +=
							coeffs( sample ) * ( m_design_matrix.row( sample ).transpose() * m_design_matrix.row( sample ) ) ;
						RowVector second_term = V( sample ) * m_design_matrix.row( sample ) ;
						m_value_of_second_derivative -= second_term.transpose() * second_term ;
					}
					++sample ; // skip excluded sample.
				}
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
