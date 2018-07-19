
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
			Matrix const& outcome, SampleRanges const& nonmissing_outcome, std::vector< std::string > const& outcome_names,
			Matrix const& covariates, SampleRanges const& nonmissing_covariate, std::vector< std::string > const& covariate_names,
			std::vector< std::string > const& predictor_names,
			Transform transform,
			std::vector< int > const& interacting_covariates
		) {
			return Design::UniquePtr(
				new Design(
					outcome, nonmissing_outcome, outcome_names,
					covariates, nonmissing_covariate, covariate_names,
					predictor_names,
					transform,
					interacting_covariates
				)
			) ;
		}

		Design::UniquePtr Design::create(
			Matrix const& outcome, SampleRanges const& nonmissing_outcome, std::vector< std::string > const& outcome_names,
			std::vector< std::string > const& predictor_names,
			Transform transform
		) {
			return Design::UniquePtr(
				new Design(
					outcome, nonmissing_outcome, outcome_names,
					predictor_names,
					transform
				)
			) ;
		}

		Design::Design(
			Matrix const& outcome, SampleRanges const& nonmissing_outcome,  std::vector< std::string > const& outcome_names,
			Matrix const& covariates, SampleRanges const& nonmissing_covariates, std::vector< std::string > const& covariate_names,
			std::vector< std::string > const& predictor_names,
			Transform transform,
			std::vector< int > const& interacting_covariates
		):
			// outcome variables
		 	m_outcome( outcome ),
			m_outcome_names( outcome_names ),
		 	m_nonmissing_outcome( nonmissing_outcome ),
			// predictor variables
			m_predictor_level_probabilities( Matrix::Zero( outcome.rows(), 0 )),
			m_predictor_levels( Matrix::Zero( 0, predictor_names.size() )),
			m_predictor_names( predictor_names ),
			m_nonmissing_predictors( SampleRanges() ),
			// covariate variables
			m_covariates( 1, covariates ),
			m_covariate_names( covariate_names ),
			m_nonmissing_covariates( nonmissing_covariates ),
			// predictor-covariate interactions
			m_predictor_covariate_interactions( interacting_covariates ),
			// Design matrix column transform
			m_transform( transform )
		{
			assert( covariates.rows() == m_outcome.rows() ) ;
			assert( m_outcome_names.size() == m_outcome.cols() ) ;
			recalculate() ;
		}

		Design::Design(
			Matrix const& outcome, SampleRanges const& nonmissing_outcome,  std::vector< std::string > const& outcome_names,
			std::vector< std::string > const& predictor_names,
			Transform transform
		):
			// outcome variables
		 	m_outcome( outcome ),
			m_outcome_names( outcome_names ),
		 	m_nonmissing_outcome( nonmissing_outcome ),
			// predictor variables
			m_predictor_level_probabilities( Matrix::Zero( outcome.rows(), 0 )),
			m_predictor_levels( Matrix::Zero( 0, predictor_names.size() )),
			m_predictor_names( predictor_names ),
			m_nonmissing_predictors( SampleRanges() ),
			// covariate variables
			m_covariates(),
			m_covariate_names(),
			m_nonmissing_covariates( 1, SampleRange( 0, m_outcome.rows() )),
			// predictor-covariate interactions
			m_predictor_covariate_interactions(),
			// Design matrix column transform
			m_transform( transform )
		{
			assert( m_outcome_names.size() == m_outcome.cols() ) ;
			recalculate() ;
		}
		
		// add a single continuous covariate
		// missing values are encoded by NaN
		void Design::add_single_covariate(
			std::string const& name,
			boost::function< double( std::size_t ) > const& values
		) {
			Vector data = Vector::Constant( m_outcome.rows(), 0 ) ;
			std::vector< SampleRange > nonmissingness ;
			int last_nonmissing_sample_i = 0 ;
			for( int i = 0; i < m_outcome.rows(); ++i ) {
				double v = values(i) ;
				if( v == v ) {
					data(i) = v ;
				} else {
					if( i > last_nonmissing_sample_i ) {
						nonmissingness.push_back( SampleRange( last_nonmissing_sample_i, i )) ;
					}
					last_nonmissing_sample_i = i+1 ;
				}
			}
			if( last_nonmissing_sample_i < m_outcome.rows() ) {
				nonmissingness.push_back( SampleRange( last_nonmissing_sample_i, m_outcome.rows() )) ;
			}

			// Ok push it.
			m_covariates.push_back( data ) ;
			m_covariate_names.push_back( name ) ;
			m_nonmissing_covariates = impl::intersect_ranges( m_nonmissing_covariates, nonmissingness ) ;
			
			recalculate() ;
		}
		
		// add a single discrete covariate (which expands into multiple columns of zeroes and ones)
		// levels are specified by non-negative integers.
		// missing values are encoded by -1
		void Design::add_discrete_covariate(
			std::string const& name,
			boost::function< int( std::size_t ) > const& values,
			boost::function< std::string( int ) > const& levelNames,
			int numberOfLevels
		) {
			// column for level 0 is not provided as will be a linear combination of the others.
			Matrix data = Matrix::Constant( m_outcome.rows(), numberOfLevels - 1, 0.0 ) ;
			std::vector< SampleRange > nonmissingness ;
			int last_nonmissing_sample_i = 0 ;
			std::vector< std::string > covariateLevelNames ;
			for( int i = 0; i < m_outcome.rows(); ++i ) {
				int level = values( i ) ;
				if( level > 0 ) {
					data( i, level-1 ) = 1.0 ;
				} else if( level == -1 ) {
					if( i > last_nonmissing_sample_i ) {
						nonmissingness.push_back( SampleRange( last_nonmissing_sample_i, i )) ;
					}
					last_nonmissing_sample_i = i+1 ;
				}
			}
			if( last_nonmissing_sample_i < m_outcome.rows() ) {
				nonmissingness.push_back( SampleRange( last_nonmissing_sample_i, m_outcome.rows() )) ;
			}
			m_covariates.push_back( data ) ;
			for( int i = 1; i < numberOfLevels; ++i ) {
				m_covariate_names.push_back( name + "=" + levelNames(i) ) ;
			}
			m_nonmissing_covariates = impl::intersect_ranges( m_nonmissing_covariates, nonmissingness ) ;

			recalculate() ;
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
			int const P = number_of_predictors() + number_of_interaction_terms() ;

			// Compute mean predictors
			Eigen::MatrixXd mean_predictors = Eigen::MatrixXd::Zero( N, P ) ;
			{
				for( int level = 0; level < m_predictor_level_probabilities.cols(); ++level ) {
					mean_predictors += m_predictor_level_probabilities.col(level).asDiagonal() * m_predictor_levels.block( 0, P*level, N, P ) ;
				}
			}

			std::vector< std::size_t > widths( D + 1 ) ;
			widths[0] = 9 ;
			out << "         outcome   " ;
			for( int j = 0; j < D; ++j ) {
				std::string const name = m_design_matrix_column_names[j]  ;
				widths[j+1] = std::max( name.size(), 5ul ) ;
				out << std::setw( widths[j+1] ) << name << " " ;
			}
			out << "\n" ;
			for( int i = 0; i < N; ++i ) {
				out << std::setw(5) << i ;
				out << std::setw(3) << "   " ;
				if( outcome().row(i).sum() == outcome().row(i).sum() ) {
					out << outcome()(i,0) << " " << outcome()(i,1) ;
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
						out << mean_predictors(i,j-1) ;
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
			SampleRanges const& nonmissing_outcome,
			SampleRanges const& nonmissing_predictors,
			SampleRanges const& nonmissing_covariates
		) const {
			SampleRanges result = impl::intersect_ranges( nonmissing_outcome, nonmissing_predictors ) ;
			return impl::intersect_ranges( result, nonmissing_covariates ) ;
		}

		void Design::recalculate() {
			calculate_design_matrix(
				m_outcome.rows(),
				m_covariates,
				m_covariate_names,
				m_predictor_names,
				m_predictor_covariate_interactions,
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

		void Design::calculate_design_matrix(
			int const number_of_samples,
			std::vector< Matrix > const& covariates,
			std::vector< std::string > const& covariate_names,
			std::vector< std::string > const& predictor_names,
			std::vector< int > const& interacting_covariates,
			Matrix* result,
			std::vector< int >* design_matrix_interaction_cols,
			std::vector< std::string >* design_matrix_column_names
		) const {
			assert( result != 0 ) ;
			assert( design_matrix_interaction_cols != 0 ) ;
			assert( design_matrix_interaction_cols->size() == 0 ) ;
			std::size_t const number_of_predictors = predictor_names.size() ;
			// Compute the number of covariates
			std::size_t numberOfCovariates = 0 ;
			for( std::size_t i = 0; i < covariates.size(); ++i ) {
				numberOfCovariates += covariates[i].cols() ;
			}
			
			// design matrix is layed out as:
			// baseline column of 1's
			// columns for predictors with values set later
			// columns for predictor vs. covariate interactions with values set later
			// columns for covariates with values set here.
			std::size_t numberOfDesignMatrixColumns
				= 1 + number_of_predictors * ( 1 + interacting_covariates.size() ) + numberOfCovariates ;

			result->setZero( number_of_samples, numberOfDesignMatrixColumns ) ;
			design_matrix_column_names->clear() ;

			// baseline column
			result->leftCols( 1 ).setOnes() ;
			design_matrix_column_names->push_back( "baseline" ) ;

			// predict columns
			design_matrix_column_names->insert(
				design_matrix_column_names->end(),
				m_predictor_names.begin(),
				m_predictor_names.end()
			) ;

			// predictor-covariate interactions
			if( interacting_covariates.size() > 0 ) {
				for( std::size_t i = 0; i < interacting_covariates.size(); ++i ) {
					design_matrix_interaction_cols->push_back( interacting_covariates[i] + result->cols() - numberOfCovariates ) ;
					for( std::size_t predictor = 0; predictor < number_of_predictors; ++predictor ) {
						design_matrix_column_names->push_back( m_predictor_names[predictor] + "x" + covariate_names[interacting_covariates[i]] ) ;
					}
				}
			}

			// covariates
			if( numberOfCovariates > 0 ) {
				int col = numberOfDesignMatrixColumns - numberOfCovariates ;
				for( std::size_t i = 0; i < covariates.size(); ++i ) {
					result->block(
						0, col,
						result->rows(), covariates[i].cols()
					) = covariates[i] ;
					col += covariates[i].cols() ;
				}
				design_matrix_column_names->insert(
					design_matrix_column_names->end(),
					covariate_names.begin(),
					covariate_names.end()
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
			SampleRanges const& nonmissingness,
			std::vector< std::string > const& names
		) {
			assert( outcome.rows() == m_design_matrix.rows() ) ;
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
			assert( levels.cols() == number_of_predictors() ) ;
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
				int const number_of_predictors_and_interactions = number_of_predictors() + number_of_interaction_terms() ;
				m_predictor_levels = Matrix::Zero( N, number_of_predictors_and_interactions * levels.rows() ) ;

				// Lay this out as follows.
				// For predictor level level_i, we fill in the block of m_predictor_levels with indices
				// [0,N) x [ level_i*npi,(level_i+1)*npi )
				// where npi is the total number of predictor and interaction terms.
				for( int level_i = 0; level_i < levels.rows(); ++level_i ) {
					m_predictor_levels.block(
						0, level_i * number_of_predictors_and_interactions,
						N, number_of_predictors()
					).rowwise() = levels.row(level_i) ;
					for( std::size_t i = 0; i < m_design_matrix_interaction_columns.size(); ++i ) {
						m_predictor_levels.block(
							0, level_i * number_of_predictors_and_interactions + (i+1) * number_of_predictors(),
							N, number_of_predictors()
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
		
			int const number_of_predictors_and_interactions = number_of_predictors() + number_of_interaction_terms() ;
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
				int const number_of_predictors_and_interactions = ( 1 + m_design_matrix_interaction_columns.size() ) * number_of_predictors() ;
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
