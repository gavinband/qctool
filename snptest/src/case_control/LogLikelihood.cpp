#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <boost/iterator/counting_iterator.hpp>
#include "snptest/case_control/LogLikelihood.hpp"
#include "genfile/Error.hpp"

namespace snptest {
	namespace case_control {
		LogLikelihood::LogLikelihood() {}

		LogLikelihood& LogLikelihood::set_phenotypes( Vector const& phenotypes ) {
			if( m_genotype_call_probabilities.rows() > 0 || m_genotype_call_probabilities.cols() > 0 ) {
				if( !phenotypes.rows() == m_genotype_call_probabilities.rows() ) {
					throw genfile::BadArgumentError( "snptest::case_control::LogLikelihood::set_phenotypes()", "phenotypes" ) ;
				}
			}
			
			if( m_covariates.rows() > 0 || m_covariates.cols() > 0 ) {
				if( !phenotypes.rows() == m_covariates.rows() ) {
					throw genfile::BadArgumentError( "snptest::case_control::LogLikelihood::set_phenotypes()", "phenotypes" ) ;
				}
			}

			m_phenotypes = phenotypes ;
			add_exclusions( m_phenotypes ) ;

			calculate_design_matrix( m_covariates ) ;

			return *this ;
		}
		
		void LogLikelihood::add_exclusions( Vector const& v ) {
			// add covariate exclusions...
			for( int i = 0; i < v.size(); ++i ) {
				if( v(i) != v(i) ) {
					std::vector< int >::iterator where = std::lower_bound( m_exclusions.begin(), m_exclusions.end(), i ) ;
					if( *where != i ) {
						m_exclusions.insert( where, i ) ;
					}
				}
			}
			if( !std::binary_search( m_exclusions.begin(), m_exclusions.end(), v.size() ) ) {
				m_exclusions.push_back( v.size() ) ;
			}
		}

		LogLikelihood& LogLikelihood::set_covariates( Matrix const& covariates ) {
			if( m_genotype_call_probabilities.rows() > 0 || m_genotype_call_probabilities.cols() > 0 ) {
				if( !covariates.rows() == m_genotype_call_probabilities.rows() ) {
					throw genfile::BadArgumentError( "snptest::case_control::LogLikelihood::set_covariates()", "phenotypes" ) ;
				}
			}
			
			if( m_phenotypes.size() > 0 ) {
				if( !covariates.rows() == m_phenotypes.rows() ) {
					throw genfile::BadArgumentError( "snptest::case_control::LogLikelihood::set_covariates()", "phenotypes" ) ;
				}
			}
			m_covariates = covariates ;
			m_design_matrix = calculate_design_matrix( m_covariates ) ;
			add_exclusions( m_covariates ) ;
			compute_tensor_squares_of_rows( m_design_matrix, m_design_matrix_row_tensor_squares ) ;
			return *this ;
		}
		
		void LogLikelihood::add_exclusions( Matrix const& matrix ) {
			// add covariate exclusions...
			for( int i = 0; i < matrix.rows(); ++i ) {
				double const rowSum = matrix.row( i ).sum() ;
				if( rowSum != rowSum ) {
					std::vector< int >::iterator where = std::lower_bound( m_exclusions.begin(), m_exclusions.end(), i ) ;
					if( *where != i ) {
						m_exclusions.insert( where, i ) ;
					}
				}
			}
			if( !std::binary_search( m_exclusions.begin(), m_exclusions.end(), matrix.rows() )) {
				m_exclusions.push_back( matrix.rows() ) ;
			}
		}

		LogLikelihood& LogLikelihood::set_genotypes( Matrix const& genotypes, Vector const& levels ) {
			if( m_covariates.rows() > 0 || m_covariates.cols() > 0 ) {
				if( !genotypes.rows() == m_covariates.rows() ) {
					throw genfile::BadArgumentError( "snptest::case_control::LogLikelihood::set_genotypes()", "phenotypes" ) ;
				}
			}
			
			if( m_phenotypes.size() > 0 ) {
				if( !genotypes.rows() == m_phenotypes.size() ) {
					throw genfile::BadArgumentError( "snptest::case_control::LogLikelihood::set_genotypes()", "phenotypes" ) ;
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
			add_exclusions( genotypes ) ;
			
			// Check if we need to re-compute the design matrix row tensor squares.
			// (We'll need to re-compute if the number of genotype levels has changed.)
			if( m_design_matrix_row_tensor_squares.cols() != ( m_design_matrix.cols() * m_genotype_levels.size() ) ) {
				compute_tensor_squares_of_rows( m_design_matrix, m_design_matrix_row_tensor_squares ) ;
			}
			
			return *this ;
		}

		LogLikelihood::Matrix LogLikelihood::calculate_design_matrix( Matrix const& covariates ) const {
			Matrix result( m_phenotypes.rows(), covariates.cols() + 2 ) ;
			result.leftCols( 1 ).setOnes() ;
			if( covariates.cols() > 0 ) {
				result.rightCols( covariates.cols() ) = covariates ;
			}
			return result ;
		}
		
		// Construct a big matrix of size (D*N) x (D*# of genotype levels).
		// The block at position (D*sample,D*g) is the tensor square of the design matrix
		// row for the given sample and given genotype.
		void LogLikelihood::compute_tensor_squares_of_rows( Matrix& design_matrix, Matrix& result ) {
			int const D = design_matrix.cols() ;
			result.resize( D * m_phenotypes.size(), D * m_genotype_levels.size() ) ;
			for( int g = 0; g < m_genotype_levels.size(); ++g ) {
				design_matrix.col(1).setConstant( m_genotype_levels( g ) ) ;
				for( int sample = 0; sample < design_matrix.rows(); ++sample ) {
					if( sample > 0 && design_matrix.row( sample ) == design_matrix.row( sample - 1 )) {
						result.block( D * sample, D * g, D, D ) = result.block( D * ( sample - 1 ), D * g, D, D ) ;
					}
					else {
						result.block( D * sample, D * g, D, D ) = design_matrix.row( sample ).transpose() * design_matrix.row( sample ) ;
					}
				}
			}
		}
		
		void LogLikelihood::evaluate_at( Point const& parameters ) {
			calculate_outcome_probabilities( parameters, m_phenotypes, m_design_matrix, &m_outcome_probabilities ) ;
			assert( m_outcome_probabilities.rows() == m_genotype_call_probabilities.rows() && m_outcome_probabilities.cols() == m_genotype_call_probabilities.cols() ) ;

#if DEBUG_LOGLIKELIHOOD
			int const N = m_phenotypes.size() ;
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
			m_value_of_function = 0.0 ;
			int start_row = 0 ;
			for( std::size_t i = 0; i < m_exclusions.size(); ++i ) {
				int end_row = m_exclusions[i] ;
				m_value_of_function += V.block( start_row, 0, end_row - start_row, V.cols() ).rowwise().sum().array().log().sum() ;
				start_row = end_row + 1 ;
			}
		#if DEBUG_LOGLIKELIHOOD
			std::cerr << "Function value = " << m_value_of_function << ".\n" ;
		#endif
		}

		void LogLikelihood::compute_value_of_first_derivative( Matrix& V ) {
			int const N = m_phenotypes.size() ;
			int const D = m_design_matrix.cols() ;
			
			// ...and its first derivative...
			{
				RowVector ones( RowVector::Ones( m_genotype_levels.size() )) ;
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
				for( int g = 0; g < m_genotype_levels.size(); ++g ) {
					m_design_matrix.col(1).setConstant( m_genotype_levels( g ) ) ;
					int start_row = 0 ;
					for( std::size_t i = 0; i < m_exclusions.size(); ++i ) {
						int end_row = m_exclusions[i] ;
						// multiply ith row of design matrix by sign times ith entry of column of V.
						m_value_of_first_derivative += (
							(
								signs.segment( start_row, end_row - start_row ).array()
								* V.col( g ).segment( start_row, end_row - start_row ).array()
							).matrix().asDiagonal()
							* m_design_matrix.block( start_row, 0, end_row - start_row, m_design_matrix.cols() )
						).colwise().sum() ;
						start_row = end_row + 1 ;
					}
				}
			}
		}

		void LogLikelihood::compute_value_of_second_derivative( Matrix const& V ) {
			int const N = m_phenotypes.size() ;
			int const D = m_design_matrix.cols() ;

			// ...and its second derivative...
			m_value_of_second_derivative = Matrix::Zero( m_design_matrix.cols(), m_design_matrix.cols() ) ;
			Vector coeffs( N ) ;
			for( int g = 0; g < m_genotype_levels.size(); ++g ) {
				coeffs = V.col( g ).array() * ( Vector::Ones( N ) - 2.0 * m_outcome_probabilities.col( g ) ).array() ;
				m_design_matrix.col(1).setConstant( m_genotype_levels( g ) ) ;
				int sample = 0 ;
				for( std::size_t i = 0; i < m_exclusions.size(); ++i ) {
					for( ; sample < m_exclusions[i]; ++sample ) {
						double const coeff = coeffs( sample ) - V( sample, g ) * V( sample, g ) ;
						m_value_of_second_derivative += coeff * m_design_matrix_row_tensor_squares.block( D * sample, D * g, D, D ) ;
					}
					++sample ; // skip excluded sample.
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
			result->resize( phenotypes.size(), m_genotype_levels.size() ) ;
			Vector linear_combination ;
			for( int g = 0; g < m_genotype_levels.size(); ++g ) {
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
