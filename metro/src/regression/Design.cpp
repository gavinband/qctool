
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <utility>
#include <algorithm>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/bind.hpp>
#include "metro/regression/Design.hpp"
#include "metro/intersect_ranges.hpp"

// #define DEBUG_REGRESSIONDESIGN 1

namespace metro {
	namespace regression {
		Design::UniquePtr Design::create(
			Matrix const& outcome, NonmissingnessMatrix const& phenotype_nonmissingness, std::vector< std::string > const& outcome_names,
			Matrix const& covariates, Matrix const& covariate_nonmissingness, std::vector< std::string > const& covariate_names,
			std::vector< std::string > const& predictor_names,
			Transform transform,
			std::vector< int > const& interacting_covariates
		) {
			return Design::UniquePtr(
				new Design(
					outcome, phenotype_nonmissingness, outcome_names,
					covariates, covariate_nonmissingness, covariate_names,
					predictor_names,
					transform,
					interacting_covariates
				)
			) ;
		}

		Design::Design(
			Matrix const& outcome, NonmissingnessMatrix const& phenotype_nonmissingness,  std::vector< std::string > const& outcome_names,
			Matrix const& covariates, NonmissingnessMatrix const& nonmissing_covariates, std::vector< std::string > const& covariate_names,
			std::vector< std::string > const& predictor_names,
			Transform transform,
			std::vector< int > const& interacting_covariates
		):
			// outcome variables
		 	m_outcome( outcome ),
			m_outcome_names( outcome_names ),
		 	m_nonmissing_outcome( phenotype_nonmissingness ),
			// predictor variables
			m_predictor_level_probabilities( Matrix::Zero( outcome.rows(), 0 )),
			m_predictor_levels( Matrix::Zero( 0, predictor_names.size() )),
			m_predictor_names( predictor_names ),
			m_nonmissing_predictors( SampleRanges() ),
			// predictor variables
			m_covariate_names( covariate_names ),
			m_nonmissing_covariates( nonmissing_covariates ),
			// Design matrix column transform
			m_transform( transform )
		{
			assert( covariates.rows() == m_outcome.rows() ) ;
			assert( m_outcome_names.size() == m_outcome.cols() ) ;
			calculate_design_matrix(
				m_outcome.rows(),
				covariates, nonmissing_covariates,
				m_predictor_levels.cols(),
				interacting_covariates,
				&m_design_matrix,
				&m_design_matrix_interaction_columns,
				&m_design_matrix_column_names
			) ;
			// 
			m_nonmissing_samples = compute_nonmissing_samples(
				m_nonmissing_outcome, m_nonmissing_predictors, m_nonmissing_covariates
			) ;
			assert( m_design_matrix_column_names.size() == m_design_matrix.cols() ) ;
		}

		std::string const& Design::get_predictor_name( std::size_t i ) const {
			return m_design_matrix_column_names[i] ;
		}

		std::vector< std::string > const& Design::design_matrix_column_names() const {
			return m_design_matrix_column_names ;
		}

		std::string const& Design::get_outcome_name( std::size_t i ) const {
			return m_outcome_names[i] ;
		}

		std::string Design::get_summary() const {
			std::ostringstream out ;

			int const N = matrix().rows() ;
			int const D = matrix().cols() ;
			int const P = number_of_predictors() ;
	
			std::vector< std::size_t > widths( D + 1 ) ;
			widths[0] = 9 ;
			out << "        phenotype   " ;
			for( int j = 0; j < D; ++j ) {
				std::string const name = m_design_matrix_column_names[j]  ;
				widths[j+1] = std::max( name.size(), 5ul ) ;
				out << std::setw( widths[j+1] ) << name << " " ;
			}
			out << "\n" ;
			for( int i = 0; i < N; ++i ) {
				out << std::setw(5) << i ;
				out << std::setw(3) << "   " ;
				if( outcome()(i) == outcome()(i) ) {
					out << std::setw(widths[0]) << outcome()(i) ;
				} else {
					out << std::setw(widths[0]) << "NA" ;
				}
				out << (( i == 6 ) ? " ~ " : "   " ) ;
				for( int j = 0; j < D; ++j ) {
					if( j > 0 ) {
						out << std::setw(1) << " " ;
					}
					out << std::setw( widths[j+1] ) ;
					if( j > 0 & j <= P ) {
						out << "?" ;
					} else if( matrix()(i,j) == matrix()(i,j) ){
						out << matrix()(i,j) ;
					} else {
						out << "NA" ;
					}
				}
				out << "\n" ;
			
				if( i == 6 && N > 8 ) {
					for( int i2 = 0; i2 < 3; ++i2 ) {
						out << std::setw(5) << " " << std::setw(3) << "   " ;
						out << std::setw(widths[0]) << "." ;
						out << "   " ;
						for( int j2 = 0; j2 < D; ++j2 ) {
							if( j2 > 0 ) {
								out << std::setw(1) << " " ;
							}
							out << std::setw( widths[j2+1] ) << "." ;
						}
						out << "\n" ;
					}
					i = N - 5 ;
				}
			}
			return out.str() ;
		}

		Design::SampleRanges Design::compute_nonmissing_samples(
			NonmissingnessMatrix const& nonmissing_outcome,
			SampleRanges const& nonmissing_predictors,
			NonmissingnessMatrix const& nonmissing_covariates
		) const {
			SampleRanges result ;
			int inclusion_range_start = 0 ;
			bool in_inclusion_range = true ;
			for( int i = 0; i < m_outcome.rows(); ++i ) {
				if(
					nonmissing_outcome.row(i).minCoeff() == 0
					|| ( nonmissing_covariates.cols() > 0 && nonmissing_covariates.row( i ).minCoeff() == 0 )
				) {
					// missing outcome or covariates => end the current range of nonmissing samples
					if( in_inclusion_range ) {
						result.push_back( metro::SampleRange( inclusion_range_start, i ) ) ;
					}
					in_inclusion_range = false ;
					inclusion_range_start = i + 1 ;
				} else {
					in_inclusion_range = true ;
				}
			}
			if( in_inclusion_range ) {
				result.push_back( metro::SampleRange( inclusion_range_start, int( m_outcome.rows() ))) ;
			}
			result = impl::intersect_ranges( result, nonmissing_predictors ) ;
			return result ;
		}

		void Design::calculate_design_matrix(
			int const number_of_samples,
			Matrix const& covariates, NonmissingnessMatrix const& covariate_nonmissingness,
			int const number_of_predictors,
			std::vector< int > const& interacting_covariates,
			Matrix* result,
			std::vector< int >* design_matrix_interaction_cols,
			std::vector< std::string >* design_matrix_column_names
		) const {
			assert( result != 0 ) ;
			assert( design_matrix_interaction_cols != 0 ) ;
			assert( design_matrix_interaction_cols->size() == 0 ) ;
			// design matrix is layed out as:
			// baseline column of 1's
			// columns for predictors with values set later
			// columns for predictor vs. covariate interactions with values set later
			// columns for covariates with values set here.
			result->setZero(
				number_of_samples,
				1 + number_of_predictors * ( 1 + interacting_covariates.size() )
				+ covariates.cols()
			) ;
			result->leftCols( 1 ).setOnes() ;
			design_matrix_column_names->push_back( "baseline" ) ;
			design_matrix_column_names->insert(
				design_matrix_column_names->end(),
				m_predictor_names.begin(),
				m_predictor_names.end()
			) ;

			if( interacting_covariates.size() > 0 ) {
				for( std::size_t i = 0; i < interacting_covariates.size(); ++i ) {
					design_matrix_interaction_cols->push_back( interacting_covariates[i] + result->cols() - covariates.cols() ) ;
					for( std::size_t predictor = 0; predictor < m_predictor_names.size(); ++predictor ) {
						design_matrix_column_names->push_back( m_predictor_names[predictor] + "x" + m_covariate_names[interacting_covariates[i]] ) ;
					}
				}
			}

			if( covariates.cols() > 0 ) {
				result->rightCols( covariates.cols() ) = covariates ;
				design_matrix_column_names->insert(
					design_matrix_column_names->end(),
					m_covariate_names.begin(),
					m_covariate_names.end()
				) ;
			}
	#if DEBUG_REGRESSIONDESIGN
			std::cerr << "Design::calculate_design_matrix():"
				<< "result has " << result->cols() << " columns, and " << design_matrix_column_names->size() << " names.\n" ;
	#endif
		
		
			assert( design_matrix_column_names->size() == result->cols() ) ;
		}
	
		Design& Design::set_outcome(
			Matrix const& outcome,
			NonmissingnessMatrix const& nonmissingness,
			std::vector< std::string > const& names
		) {
			assert( outcome.rows() == m_design_matrix.rows() ) ;
			assert( nonmissingness.size() == outcome.rows() ) ;
			assert( names.size() == outcome.cols() ) ;

			m_outcome = outcome ;
			m_outcome_names = names ;
			if( nonmissingness != m_nonmissing_outcome ) {
				m_nonmissing_outcome = nonmissingness ;
				m_nonmissing_samples = compute_nonmissing_samples( m_nonmissing_outcome, m_nonmissing_predictors, m_nonmissing_covariates ) ;
			}
			return *this ;
		}
	
		Design& Design::set_predictors(
			Matrix const& levels,
			Matrix const& probabilities,
			std::vector< metro::SampleRange > const& nonmissingness
		) {
			assert( levels.cols() == m_predictor_levels.cols() ) ;
			assert( probabilities.cols() == levels.rows() ) ;
			assert( probabilities.rows() == m_design_matrix.rows() ) ;

			// Set probabilities
			m_predictor_level_probabilities = probabilities ;

			// Set levels
			{
				// This is complexified by the need to handle possible interactions. 
				// There is one column per predictor and one column per predictor per interacting covariate.
				// The levels go in a matrix of dimension N x ((#predictors & interactions) x (#levels)).
				int const N = probabilities.rows() ;
				int const number_of_predictors_and_interactions = ( 1 + m_design_matrix_interaction_columns.size() ) * m_predictor_levels.cols() ;
				m_predictor_levels = Matrix::Zero( N, number_of_predictors_and_interactions * levels.rows() ) ;

				// Lay this out as follows.
				// For predictor level level_i, we fill in the block of m_predictor_levels with indices
				// [0,N) x [ level_i*npi,(level_i+1)*npi )
				// where npi is the total number of predictor and interaction terms.
				for( int level_i = 0; level_i < levels.rows(); ++level_i ) {
					m_predictor_levels.block(
						0, level_i * number_of_predictors_and_interactions,
						N, levels.cols()
					).rowwise() = levels.row(level_i) ;
					for( std::size_t i = 0; i < m_design_matrix_interaction_columns.size(); ++i ) {
						m_predictor_levels.block(
							0, level_i * number_of_predictors_and_interactions + (i+1)*m_predictor_levels.cols(),
							N, m_predictor_levels.cols()
						) = m_design_matrix.col( m_design_matrix_interaction_columns[i] ) * levels.row(level_i) ;
					}
				}
				
				if( m_transform == eMeanCentre ) {
					mean_centre_predictor_levels(
						m_predictor_level_probabilities,
						&m_predictor_levels,
						nonmissingness
					) ;
				}
			}

			// Set nonmissingness if necessary
			if( nonmissingness != m_nonmissing_predictors ) {
				m_nonmissing_predictors = nonmissingness ;
				m_nonmissing_samples = compute_nonmissing_samples( m_nonmissing_outcome, m_nonmissing_predictors, m_nonmissing_covariates ) ;
			}

			
			return *this ;
		}

		void Design::mean_centre_predictor_levels(
			Matrix const& probs,
			Matrix* predictor_levels,
			SampleRanges const& included_samples
		) const {
			assert( probs.cols() == m_predictor_levels.rows() ) ;
		
			int const number_of_predictors_and_interactions = ( 1 + m_design_matrix_interaction_columns.size() ) * m_predictor_levels.cols() ;
			RowVector means = RowVector::Zero( number_of_predictors_and_interactions ) ;
			double total = 0 ;
			for( std::size_t i = 0; i < included_samples.size(); ++i ) {
				int const start = included_samples[i].begin() ;
				int const end = included_samples[i].end() ;
				for( int level_i = 0; level_i < m_predictor_levels.rows(); ++level_i ) {
					means += (
						( probs.col( level_i ).segment( start, end - start ) ).asDiagonal()
						* predictor_levels->block(
							start, number_of_predictors_and_interactions * level_i,
							end - start, number_of_predictors_and_interactions
						)
					).colwise().sum() ;
				}
			
				total += probs.block( start, 0, end - start, probs.cols() ).sum() ;
			}
			means /= total ;
			for( int level_i = 0; level_i < m_predictor_levels.rows(); ++level_i ) {
				predictor_levels->block(
					0, number_of_predictors_and_interactions * level_i,
					predictor_levels->rows(), number_of_predictors_and_interactions
				).rowwise() -= means ;
			}
		}

		Design& Design::set_predictor_level( int level ) {
			if( m_predictor_levels.cols() > 0 ) {
				int const number_of_predictors_and_interactions = ( 1 + m_design_matrix_interaction_columns.size() ) * m_predictor_levels.cols() ;
				m_design_matrix.block(
					0, 1,
					m_design_matrix.rows(), number_of_predictors_and_interactions
				) = m_predictor_levels.block(
					0, level * number_of_predictors_and_interactions,
					m_design_matrix.rows(), number_of_predictors_and_interactions
				) ;
			}
			return *this ;
		}
	}
}
