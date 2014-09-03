
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <set>
#include <vector>
#include <boost/iterator/counting_iterator.hpp>
#include "snptest/case_control/NullModelLogLikelihood.hpp"
#include "genfile/Error.hpp"

//#define DEBUG_LOGLIKELIHOOD 1

namespace snptest {
	namespace case_control {
		NullModelLogLikelihood::NullModelLogLikelihood() {}

		NullModelLogLikelihood& NullModelLogLikelihood::set_phenotypes( Vector const& phenotypes ) {
			if( m_genotype_call_probabilities.rows() > 0 || m_genotype_call_probabilities.cols() > 0 ) {
				if( phenotypes.rows() != m_genotype_call_probabilities.rows() ) {
					throw genfile::BadArgumentError( "snptest::case_control::NullModelLogLikelihood::set_phenotypes()", "phenotypes" ) ;
				}
			}
			
			if( m_covariates.rows() > 0 || m_covariates.cols() > 0 ) {
				if( phenotypes.rows() != m_covariates.rows() ) {
					throw genfile::BadArgumentError( "snptest::case_control::NullModelLogLikelihood::set_phenotypes()", "phenotypes" ) ;
				}
			} else {
				m_covariates = Matrix( phenotypes.rows(), 0 ) ;
			}

			m_phenotypes = phenotypes ;
			add_exclusions( m_phenotypes ) ;


			return *this ;
		}
		
		void NullModelLogLikelihood::add_exclusions( Vector const& v ) {
			// add covariate exclusions...
			for( int i = 0; i < v.size(); ++i ) {
				if( v(i) != v(i) ) {
					std::vector< int >::iterator where = std::lower_bound( m_exclusions.begin(), m_exclusions.end(), i ) ;
					if( where == m_exclusions.end() || *where != i ) {
						m_exclusions.insert( where, i ) ;
					}
				}
			}
			if( !std::binary_search( m_exclusions.begin(), m_exclusions.end(), v.size() ) ) {
				m_exclusions.push_back( v.size() ) ;
			}
		}

		NullModelLogLikelihood& NullModelLogLikelihood::set_covariates( Matrix const& covariates ) {
			if( m_genotype_call_probabilities.rows() > 0 || m_genotype_call_probabilities.cols() > 0 ) {
				if( covariates.rows() != m_genotype_call_probabilities.rows() ) {
					throw genfile::BadArgumentError( "snptest::case_control::NullModelLogLikelihood::set_covariates()", "phenotypes" ) ;
				}
			}
			
			if( m_phenotypes.size() > 0 ) {
				if( covariates.rows() != m_phenotypes.size() ) {
					throw genfile::BadArgumentError( "snptest::case_control::NullModelLogLikelihood::set_covariates()", "phenotypes" ) ;
				}
			}
			m_covariates = covariates ;
			m_design_matrix = calculate_design_matrix( m_covariates ) ;
			add_exclusions( m_covariates ) ;
			return *this ;
		}
		
		NullModelLogLikelihood& NullModelLogLikelihood::set_predictor_probs( Matrix const& genotypes, Vector const& levels ) {
			if( m_covariates.rows() > 0 || m_covariates.cols() > 0 ) {
				if( genotypes.rows() != m_covariates.rows() ) {
					throw genfile::BadArgumentError( "snptest::case_control::NullModelLogLikelihood::set_predictor_probs()", "phenotypes" ) ;
				}
			}
			
			if( m_phenotypes.size() > 0 ) {
				if( genotypes.rows() != m_phenotypes.size() ) {
					throw genfile::BadArgumentError( "snptest::case_control::NullModelLogLikelihood::set_predictor_probs()", "phenotypes" ) ;
				}
			}

			assert( levels.size() == genotypes.cols() ) ;
			m_genotype_call_probabilities = genotypes ;
			m_genotype_levels = levels ;
			
			// treat all-zero genotype calls as really missing...
			for( int i = 0; i < m_genotype_call_probabilities.rows(); ++i ) {
				if( m_genotype_call_probabilities.row( i ).maxCoeff() <= 0.0 ) {
					m_genotype_call_probabilities.row( i ).setConstant( std::numeric_limits< double >::quiet_NaN() ) ;
				}
			}
			
			m_exclusions.clear() ;
			add_exclusions( m_phenotypes ) ;
			add_exclusions( m_covariates ) ;
			add_exclusions( m_genotype_call_probabilities ) ;
			
			if( m_design_matrix.rows() != m_phenotypes.size() ) {
				m_design_matrix = calculate_design_matrix( m_covariates ) ;
			}
			return *this ;
		}
		
		void NullModelLogLikelihood::add_exclusions( Matrix const& matrix ) {
			// add covariate exclusions...
			for( int i = 0; i < matrix.rows(); ++i ) {
				double const rowSum = matrix.row( i ).sum() ;
				if( rowSum != rowSum ) {
					std::vector< int >::iterator where = std::lower_bound( m_exclusions.begin(), m_exclusions.end(), i ) ;
					if( where == m_exclusions.end() || *where != i ) {
						m_exclusions.insert( where, i ) ;
					}
				}
			}
			if( !std::binary_search( m_exclusions.begin(), m_exclusions.end(), matrix.rows() )) {
				m_exclusions.push_back( matrix.rows() ) ;
			}
		}

		NullModelLogLikelihood& NullModelLogLikelihood::add_exclusions( std::vector< int > const& exclusions ) {
			for( std::size_t i = 0; i < exclusions.size(); ++i ) {
				std::vector< int >::iterator where = std::lower_bound( m_exclusions.begin(), m_exclusions.end(), exclusions[i] ) ;
				if( where != m_exclusions.end() && *where != exclusions[i] ) {
					m_exclusions.insert( where, exclusions[i] ) ;
				}
			}
			return *this ;
		}

		NullModelLogLikelihood::Matrix NullModelLogLikelihood::calculate_design_matrix( Matrix const& covariates ) const {
			Matrix result( m_phenotypes.rows(), covariates.cols() + 1 ) ;
			result.leftCols( 1 ).setOnes() ;
			if( covariates.cols() > 0 ) {
				result.rightCols( covariates.cols() ) = covariates ;
			}
			return result ;
		}

		void NullModelLogLikelihood::evaluate_at( Vector const& parameters ) {
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
					m_value_of_function += V.segment( start_row, end_row - start_row ).array().log().sum() ;
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
					//std::cerr << "i = " << i << ", V.rows() = " << V.rows() << ", start_row = " << start_row << ", end_row = " << end_row << ".\n" ;
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
