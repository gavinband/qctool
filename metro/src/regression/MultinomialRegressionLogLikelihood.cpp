
//			Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//	  (See accompanying file LICENSE_1_0.txt or copy at
//			http://www.boost.org/LICENSE_1_0.txt)

#define EIGEN_RUNTIME_NO_MALLOC
#include <iostream>
#include <numeric>
#include <iomanip>
#include <vector>
#include <utility>
#include <algorithm>
#include <boost/iterator/counting_iterator.hpp>
#include "metro/regression/MultinomialRegressionLogLikelihood.hpp"
#include "metro/intersect_ranges.hpp"
#include "genfile/Error.hpp"
#include "genfile/string_utils/string_utils.hpp"

// #define DEBUG_LOGLIKELIHOOD 1

namespace metro {
	namespace regression {
		
		namespace {
			template< typename T >
			std::ostream& print_matrix( std::string const& name, T const& matrix, std::ostream& out ) {
				return out << name << " = (" << matrix.rows() << "x" << matrix.cols() << ")\n"
					<< matrix << "\n" ;
			}
		}
		
		MultinomialRegressionLogLikelihood::UniquePtr MultinomialRegressionLogLikelihood::create(
			regression::Design::UniquePtr design
		) {
			return MultinomialRegressionLogLikelihood::UniquePtr(
				new MultinomialRegressionLogLikelihood( design )
			) ;
		}
		
		MultinomialRegressionLogLikelihood::MultinomialRegressionLogLikelihood(
			regression::Design::UniquePtr design
		):
			m_design( design ),
			m_number_of_samples( m_design->outcome().rows() ),
			m_parameter_vector( Vector::Zero( m_design->matrix().cols() * ( m_design->outcome().cols() - 1 ) )),
			m_parameter_matrix( Matrix::Zero( m_design->matrix().cols(), m_design->outcome().cols() - 1 ) ),
			m_parameter_identity( compute_parameter_identity() ),
			m_numberOfDerivativesComputed( 0 )
		{
			assert( m_design.get() != 0 ) ;
		}

		std::string MultinomialRegressionLogLikelihood::get_summary() const {
			using genfile::string_utils::to_string ;
			return "MultinomialRegressionLogLikelihood( " + to_string( m_number_of_samples ) + " samples, " + to_string( m_design->outcome().cols() ) + " outcomes )" ;
		}

		void MultinomialRegressionLogLikelihood::set_predictor_levels(
			Matrix const& levels,
			Matrix const& probabilities,
			std::vector< metro::SampleRange > const& included_samples
		) {
			m_design->set_predictors( levels, probabilities, included_samples ) ;

#if DEBUG_LOGLIKELIHOOD
			int const print_N = std::min( m_design->outcome().rows(), Vector::Index( 10 ) ) ;
			int const number_of_outcomes = m_design->outcome().cols() ;
			std::cerr << "MultinomialRegressionLogLikelihood::number of outcomes = " << number_of_outcomes << "\n" ;
			std::cerr << "MultinomialRegressionLogLikelihood::set_predictor_levels(): Included samples:" ;
			for( std::size_t i = 0; i < included_samples.size(); ++i ) {
				std::cerr << " " << included_samples[i].begin() << "-" << included_samples[i].end() ;
			}
			std::cerr << "\n" ;
			Matrix const& design_matrix = m_design->set_predictor_level(0).matrix() ;
			std::cerr << "MultinomialRegressionLogLikelihood::design matrix = (" << design_matrix.rows() << "x" << design_matrix.cols() << ")\n"
				<< design_matrix.block( 0, 0, print_N, design_matrix.cols() ) << "\n" ;
			std::cerr << "MultinomialRegressionLogLikelihood::predictor levels = (" << levels.rows() << "x" << levels.cols() << ")\n"
				<< levels << "\n" ;
#endif
			compute_psi( m_design->outcome(), levels.rows(), included_samples, &m_psi ) ;
			compute_rearranger( levels.rows(), &m_outcome_wise_to_predictor_wise_rearranger ) ;
			m_numberOfDerivativesComputed = 0 ;

#if DEBUG_LOGLIKELIHOOD
			std::cerr << "MultinomialRegressionLogLikelihood:" ;
			std::cerr << "outcome =\n" << m_design->outcome().head( print_N ) << "\n" ;
			print_matrix( "psi", m_psi.block( 0, 0, print_N, m_psi.cols() ), std::cerr ) ;
#endif
		}
		
		void MultinomialRegressionLogLikelihood::compute_psi(
			Matrix const& outcome,
			int const number_of_predictor_levels,
			std::vector< metro::SampleRange > const& included_samples,
			Matrix* result
		) const {
			// psi is N x LM, thought of as M blocks of NxL.
			// For each sample, we put a row of L 1's in the columns corresponding to the outcome for that sample.
			assert( result != 0 ) ;
			
			int const number_of_outcomes = outcome.cols() ;
			result->setZero( m_number_of_samples, number_of_outcomes * number_of_predictor_levels ) ;
			for( int i = 0; i < included_samples.size(); ++i ) {
				int const start_row = included_samples[i].begin() ;
				int const end_row = included_samples[i].end() ;
				for( int sample = start_row; sample < end_row; ++sample ) {
					// look for the 1 identifying the outcome for this sample
					Matrix::Index outcome_i = 0 ;
					double const maxCoeff = outcome.row(sample).array().maxCoeff(&outcome_i) ;
					result->block( sample, outcome_i * number_of_predictor_levels, 1, number_of_predictor_levels ).setConstant( 1.0 ) ;
				}
			}
		}
		
		void MultinomialRegressionLogLikelihood::compute_rearranger( int const number_of_predictor_levels, PermutationMatrix* result ) const {
			// rearranger is MLxML.	 It is a square matrix which when multiplied on the right turns a matrix like
			// f1(x1) f2(x1) ... fM(x1) f1(x2)...etc (L blocks of NxM)
			// into
			// f1(x1) f1(x2) ... f1(xl) f2(x1)...etc (M blocks of NxL).
			int const number_of_outcomes = m_design->outcome().cols() ;
			result->setIdentity(
				number_of_predictor_levels * number_of_outcomes
			) ;
			int column = 0 ;
			for( int j = 0; j < number_of_outcomes; ++j ) {
				for( int level_i = 0; level_i < number_of_predictor_levels; ++level_i, ++column ) {
//					result->indices().data()[( level_i * number_of_outcomes + j )] = column ;
					result->indices().data()[ column ] = ( level_i * number_of_outcomes ) + j ;
				}
			}
#if DEBUG_LOGLIKELIHOOD			
			std::cerr << "Rearranger2:\n" << result->toDenseMatrix() << "\n" ;
#endif
		}
		
		MultinomialRegressionLogLikelihood::Vector const& MultinomialRegressionLogLikelihood::get_parameters() const {
			return m_parameter_vector ;
		}
		
		void MultinomialRegressionLogLikelihood::evaluate_at(
			Vector const& parameters,
			int const numberOfDerivatives
		) {
#if DEBUG_LOGLIKELIHOOD
			std::cerr << "MultinomialRegressionLogLikelihood::evaluate_at(): parameters = " << parameters.transpose() << "\n" ;
#endif
			int const number_of_outcomes = m_design->outcome().cols() ;
			assert( parameters.size() == m_design->matrix().cols() * ( number_of_outcomes - 1 ) ) ;
			evaluate_at_impl( parameters, numberOfDerivatives ) ;
		}

		void MultinomialRegressionLogLikelihood::set_parameter_naming_scheme( MultinomialRegressionLogLikelihood::GetParameterName get_parameter_name ) {
			assert( get_parameter_name ) ;
			m_get_parameter_name = get_parameter_name ;
		}

		std::string MultinomialRegressionLogLikelihood::get_parameter_name( std::size_t i ) const {
			IntegerMatrix const identity = identify_parameters() ;
			return m_get_parameter_name( m_design->get_predictor_name( identity(i,1)), identity(i,0) ) ;
		}

		LogLikelihood::IntegerMatrix MultinomialRegressionLogLikelihood::identify_parameters() const {
			return m_parameter_identity ;
		}

		int MultinomialRegressionLogLikelihood::number_of_outcomes() const {
			return m_design->outcome().cols() ;
		}

		LogLikelihood::IntegerMatrix MultinomialRegressionLogLikelihood::compute_parameter_identity() const {
			// each parameter corresponds to a column of the design matrix and a non-baseline phenotype level.
			IntegerMatrix result = IntegerMatrix::Zero( m_parameter_vector.size(), 2 ) ;
			int const D = m_design->matrix().cols() ;
			int const number_of_outcomes = m_design->outcome().cols() ;
			for( int j = 0; j < D; ++j ) {
				for( int i = 1; i < number_of_outcomes; ++i ) {
					result( (i-1)*D + j, 0 ) = i ;
					result( (i-1)*D + j, 1 ) = j ;
				}
			}
			return result ;
		}

		void MultinomialRegressionLogLikelihood::evaluate_at_impl(
			Vector const& parameters,
			int const numberOfDerivatives
		) {
			std::vector< metro::SampleRange > const& included_samples = m_design->nonmissing_samples() ;
			if( m_parameter_vector == parameters && m_numberOfDerivativesComputed > numberOfDerivatives ) {
				return ;
			}
			m_parameter_vector = parameters ;
			m_numberOfDerivativesComputed = 0 ;
			int const number_of_outcomes = m_design->outcome().cols() ;
			assert( m_parameter_vector.size() == m_parameter_identity.rows() ) ;
			assert( numberOfDerivatives < 3 ) ;
			m_parameter_matrix.resize( m_design->matrix().cols(), number_of_outcomes - 1 ) ; 
			for( int i = 0; i < m_parameter_vector.size(); ++i ) {
				m_parameter_matrix( m_parameter_identity( i, 1 ), m_parameter_identity( i, 0 ) - 1 ) = m_parameter_vector(i) ;
			}

#if DEBUG_LOGLIKELIHOOD
			std::cerr << "parameters = " << m_parameter_vector.transpose() << ".\n" ;
			print_matrix( "parameter matrix", m_parameter_matrix, std::cerr ) ;
#endif
			compute_F( m_parameter_matrix, &m_F ) ;
			compute_A_and_function_value( m_F, m_design->get_predictor_level_probabilities(), included_samples, &m_A, &m_value_of_function ) ;
#if DEBUG_LOGLIKELIHOOD
			int const print_N = std::min( m_design->outcome().rows(), Vector::Index( 10 ) ) ;
			print_matrix( "F", m_F.block( 0, 0, print_N, m_F.cols() ), std::cerr ) ;
			print_matrix( "A", m_A.block( 0, 0, print_N, m_A.cols() ), std::cerr ) ;
			print_matrix( "psi", m_psi.block( 0, 0, print_N, m_psi.cols() ), std::cerr ) ;
			std::cerr << "function value = \n" << m_value_of_function << "\n" ;
#endif
			++m_numberOfDerivativesComputed ;
			if( numberOfDerivatives > 0 ) {
				compute_B( m_A, m_F, m_psi, &m_B ) ;
				compute_value_of_first_derivative( m_B, included_samples, &m_value_of_first_derivative, &m_first_derivative_terms ) ;
				++m_numberOfDerivativesComputed ;
#if DEBUG_LOGLIKELIHOOD
				print_matrix( "B", m_B.block( 0, 0, print_N, m_B.cols() ), std::cerr ) ;
				print_matrix( "C", m_C.block( 0, 0, print_N, m_C.cols() ), std::cerr ) ;
				std::cerr << "derivative = \n" << m_value_of_first_derivative << "\n" ;
#endif
				if( numberOfDerivatives > 1 ) {
					compute_C( m_A, m_F, m_psi, included_samples, &m_C ) ;
					compute_value_of_second_derivative( m_B, m_C, included_samples, &m_value_of_second_derivative ) ;
					++m_numberOfDerivativesComputed ;
#if DEBUG_LOGLIKELIHOOD
					print_matrix( "C", m_C.block( 0, 0, print_N, m_C.cols() ), std::cerr ) ;
				
					std::cerr << "2nd derivative = \n"
						<< m_value_of_second_derivative << ".\n" ;
#endif
				}
			}
		}
		
		void MultinomialRegressionLogLikelihood::compute_F( Matrix const& parameters, Matrix* result ) const {
			// F is NxLM, though of as (M+1) blocks of L.
			// We compute it as L blocks of (M+1), since this leads to simpler code, and then rearrange columns.
			int const number_of_levels = m_design->get_number_of_predictor_levels() ;
			int const number_of_outcomes = m_design->outcome().cols() ;
			result->setZero( m_number_of_samples, number_of_levels *  number_of_outcomes ) ;
			Vector sums = Vector::Zero( m_number_of_samples ) ;
			for( int level_i = 0; level_i < number_of_levels; ++level_i ) {
				Matrix const& design_matrix = m_design->set_predictor_level( level_i ).matrix() ;
				int const block_first_column = level_i * number_of_outcomes ;
				MatrixBlock block = result->block( 0, block_first_column, m_number_of_samples, number_of_outcomes ) ;
				// numerator for baseline outcome is 1
				block.col(0).setConstant( 1.0 ) ;
				// compute numerator for non-baseline outcomes
				MatrixBlock alternative_outcome_block = result->block( 0, block_first_column + 1, m_number_of_samples, number_of_outcomes - 1 ) ;
				alternative_outcome_block = ( design_matrix * parameters ).array().exp() ;
				// Normalise all outcome probabilities
				sums = alternative_outcome_block.rowwise().sum() ;
				for( int col = 0; col < block.cols(); ++col ) {
					block.col( col ).array() /= ( Vector::Constant( m_number_of_samples, 1.0 ) + sums ).array() ;
				}
			}
			// We've computed F in the order
			// f1(x1) f2(x1) ... fM(x1) f1(x2) ...
			// We now re-arrange so that it comes in the order
			// f1(x1) f1(x2) ... f1(xl) f2(x1)...
			// as this fits better with the later computations.
			// Split into seperate function for profiling.
			rearrange_F( result ) ;
		}

		void MultinomialRegressionLogLikelihood::rearrange_F( Matrix* result ) const {
			// We've computed F in the order
			// f1(x1) f2(x1) ... fM(x1) f1(x2) ...
			// We now re-arrange so that it comes in the order
			// f1(x1) f1(x2) ... f1(xl) f2(x1)...
			// as this fits better with the later computations.
			(*result) *= m_outcome_wise_to_predictor_wise_rearranger ;
		}
		
		void MultinomialRegressionLogLikelihood::compute_A_and_function_value(
			Matrix const& F,
			Matrix const& predictor_probabilities,
			std::vector< metro::SampleRange > const& included_samples,
			Matrix* result,
			double* value_of_function
		) const {
			// A is NxL.
			int const number_of_predictor_levels = predictor_probabilities.cols() ;
			*value_of_function = 0.0 ;
			result->setZero( m_number_of_samples, number_of_predictor_levels ) ;
			for( std::size_t i = 0; i < included_samples.size(); ++i ) {
				int const start_row = included_samples[i].begin() ;
				int const end_row = included_samples[i].end() ;
				for( int sample = start_row; sample < end_row; ++sample ) {
					// Look for the 1 indicating which outcome this sample has
					Matrix::Index outcome_i = 0 ;
					m_design->outcome().row( sample ).array().maxCoeff( &outcome_i ) ;
					// Find the F values corresponding to this sample's outcome
					ConstMatrixBlock block_of_F = F.block( sample, outcome_i * number_of_predictor_levels, 1, number_of_predictor_levels ) ;
					result->row( sample ) = predictor_probabilities.row( sample ).array() * block_of_F.array() ;
				}
				// Compute the log-likelihood.
				*value_of_function += result->block( start_row, 0, end_row - start_row, number_of_predictor_levels ).rowwise().sum().array().log().sum() ;

				// Now renormalise this block of A.
				for( int sample = start_row; sample < end_row; ++sample ) {
					result->row( sample ) /= result->row( sample ).sum() ;
				}
			}
		}

		void MultinomialRegressionLogLikelihood::compute_B( Matrix const& A, Matrix const& F, Matrix const& psi, Matrix* result ) const {
			// B is NxLM, thought of as M blocks of L.
			int const number_of_levels = m_design->get_number_of_predictor_levels() ;
			int const M = m_design->outcome().cols() - 1 ;
			result->setZero( m_number_of_samples, number_of_levels * M ) ;
			for( int j = 0; j < M; ++j ) {
				MatrixBlock block = result->block( 0, j * number_of_levels, m_number_of_samples, number_of_levels ) ;
				// We use j+1 here to skip the baseline columns in F and psi.
				ConstMatrixBlock block_of_F = F.block( 0, (j+1) * number_of_levels, m_number_of_samples, number_of_levels ) ;
				ConstMatrixBlock block_of_psi = psi.block( 0, (j+1) * number_of_levels, m_number_of_samples, number_of_levels ) ;
				block = A.array() * ( block_of_psi - block_of_F ).array() ;
			}
		}

		void MultinomialRegressionLogLikelihood::compute_value_of_first_derivative(
			Matrix const& B,
			std::vector< metro::SampleRange > const& included_samples,
			Vector* result,
			Matrix* terms
		) const {
			int const D = m_design->matrix().cols() ;
			int const number_of_levels = m_design->get_number_of_predictor_levels() ;
			int const M = m_design->outcome().cols() - 1 ;
			result->setZero( D * M ) ;
			terms->setZero( m_design->matrix().rows(), D*M ) ;
			for( int level_i = 0; level_i < number_of_levels; ++level_i ) {
				Matrix const& design_matrix = m_design->set_predictor_level( level_i ).matrix() ;
				for( int i = 0; i < included_samples.size(); ++i ) {
					int const start_row = included_samples[i].begin() ;
					int const end_row = included_samples[i].end() ;
					for( int j = 0; j < M; ++j ) {
						MatrixBlock block = terms->block( start_row, j*D, end_row - start_row, D ) ;
						// multiply each row of the design matrix by the relevant term in the relevant column of B.
						block += (
							B.col( j * number_of_levels + level_i ).segment( start_row, end_row - start_row ).asDiagonal()
							* design_matrix.block( start_row, 0, end_row - start_row, design_matrix.cols() )
						) ;
					}
				}
			}
			(*result) = m_first_derivative_terms.colwise().sum() ;
		}

		void MultinomialRegressionLogLikelihood::compute_C(
			Matrix const& A,
			Matrix const& F,
			Matrix const& psi,
			std::vector< metro::SampleRange > const& included_samples,
			Matrix* result
		) const {
			// B is NxLM^2, thought of as M blocks of ML, each block is M blocks of L.
			int const number_of_levels = m_design->get_number_of_predictor_levels() ;
			int const M = m_design->outcome().cols() - 1 ;
			int const outer_blocksize = number_of_levels * M ;
			result->setZero( m_number_of_samples, outer_blocksize * M ) ;
			for( int j = 0; j < M; ++j ) {
				ConstMatrixBlock jth_block_of_F = F.block( 0, (j+1) * number_of_levels, m_number_of_samples, number_of_levels ) ;
				ConstMatrixBlock jth_block_of_psi = psi.block( 0, (j+1) * number_of_levels, m_number_of_samples, number_of_levels ) ;
				for( int k = 0; k < M; ++k ) {
					ConstMatrixBlock kth_block_of_F = F.block( 0, (k+1) * number_of_levels, m_number_of_samples, number_of_levels ) ;
					ConstMatrixBlock kth_block_of_psi = psi.block( 0, (k+1) * number_of_levels, m_number_of_samples, number_of_levels ) ;
					Matrix Zjk = -kth_block_of_F ;
					if( j == k ) {
						Zjk.array() = 1.0 + Zjk.array() ;
					}
					result->block( 0, j * outer_blocksize + k * number_of_levels, m_number_of_samples, number_of_levels ).array()
						= A.array()
						* (
							(
								( jth_block_of_psi - jth_block_of_F ).array()
								* ( kth_block_of_psi - kth_block_of_F ).array()
							)
							- (
								jth_block_of_F.array() * Zjk.array()
							)
						) ;
				}
			}
		}

		void MultinomialRegressionLogLikelihood::compute_value_of_second_derivative(
			Matrix const& B,
			Matrix const& C,
			std::vector< metro::SampleRange > const& included_samples,
			Matrix* result
		) const {
			int const D = m_design->matrix().cols() ;
			int const number_of_levels = m_design->get_number_of_predictor_levels() ;
			int const M = m_design->outcome().cols() - 1 ;
			int const outer_blocksize = M * number_of_levels ;
			result->setZero( D*M, D*M ) ;
			m_design_matrix_tensor_square_rows.resize( m_design->matrix().rows(), D ) ;
#if DEBUG_LOGLIKELIHOOD
			std::cerr << "MultinomialRegressionLogLikelihood::compute_value_of_second_derivative(): C =\n"
				<< C.block( 0, 0, std::min( int( C.rows() ), 10 ), C.cols() ) << "\n" ;
#endif
			for( int level_i = 0; level_i < number_of_levels; ++level_i ) {
				Matrix const& design_matrix = m_design->set_predictor_level( level_i ).matrix() ;
				for( int row_i = 0; row_i < D; ++row_i ) {
#ifdef EIGEN_RUNTIME_NO_MALLOC
					Eigen::internal::set_is_malloc_allowed( false ) ;
#endif
					m_design_matrix_tensor_square_rows = design_matrix.col( row_i ).asDiagonal() * design_matrix ;
#ifdef EIGEN_RUNTIME_NO_MALLOC
					Eigen::internal::set_is_malloc_allowed( true ) ;
#endif
					for( int i = 0; i < included_samples.size(); ++i ) {
						int const start_row = included_samples[i].begin() ;
						int const end_row = included_samples[i].end() ;
						MatrixBlock design_matrix_tensor_square_row_block
							= m_design_matrix_tensor_square_rows.block( start_row, 0, end_row - start_row, design_matrix.cols() ) ;
						for( int j = 0; j < M; ++j ) {
							for( int k = j; k < M; ++k ) {
								MatrixBlock resultRow = result->block( j*D+row_i, k*D, 1, D ) ;
									resultRow
										+= (
											C.block( start_row, j * outer_blocksize + k*number_of_levels + level_i, end_row - start_row, 1 ).transpose()
											* design_matrix_tensor_square_row_block
										)
									;
							}
						}
					}
				}
			}
#ifdef EIGEN_RUNTIME_NO_MALLOC
			Eigen::internal::set_is_malloc_allowed( false ) ;
#endif
			for( int i = 0; i < included_samples.size(); ++i ) {
				int const start_row = included_samples[i].begin() ;
				int const end_row = included_samples[i].end() ;
				for( int row_i = 0; row_i < M*D; ++row_i ) {
					MatrixBlock resultRow = result->block( row_i, 0, 1, M*D ) ;
					resultRow.noalias() -= (
						m_first_derivative_terms.col( row_i ).segment( start_row, end_row - start_row ).asDiagonal()
							* m_first_derivative_terms.block( start_row, 0, end_row - start_row, M*D )
						)
						.colwise().sum() ;
				}
			}
#ifdef EIGEN_RUNTIME_NO_MALLOC
			Eigen::internal::set_is_malloc_allowed( true ) ;
#endif

			// fill lower diagonal
			for( int j = 1; j < M; ++j ) {
				for( int k = 0; k < j; ++k ) {
					MatrixBlock target = result->block( j*D, k*D, D, D ) ;
					MatrixBlock resultBlock = result->block( k*D, j*D, D, D ) ;
					target.noalias() = resultBlock.transpose() ;
				}
			}
		}

		double MultinomialRegressionLogLikelihood::get_value_of_function() const {
			return m_value_of_function ;
		}

		MultinomialRegressionLogLikelihood::Vector MultinomialRegressionLogLikelihood::get_value_of_first_derivative() const {
			return m_value_of_first_derivative ;
		}

		MultinomialRegressionLogLikelihood::Matrix MultinomialRegressionLogLikelihood::get_value_of_second_derivative() const {
			return m_value_of_second_derivative ;
		}
	}
}
