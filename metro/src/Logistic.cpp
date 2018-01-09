
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
#include "metro/regression/Logistic.hpp"
#include "metro/intersect_ranges.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/Error.hpp"

// #define DEBUG_LOGLIKELIHOOD 2

namespace metro {
	namespace regression {
		
		namespace {
			std::string standard_parameter_name( std::string const& predictor_name, std::string const& outcome_name ) {
				return predictor_name + "/" + outcome_name ;
			}
		}
		
		Logistic::UniquePtr Logistic::create(
			regression::Design& design
		) {
			return Logistic::UniquePtr(
				new Logistic( design )
			) ;
		}

		Logistic::UniquePtr Logistic::create(
			regression::Design::UniquePtr design
		) {
			return Logistic::UniquePtr(
				new Logistic( design )
			) ;
		}

		namespace {
			void check_design( regression::Design const& design ) {
				// Verify design looks sane.
				// Right now we handle only a bernoulli outcome, i.e.
				// each sample must have a 1 in column 0 or 1 of the outcome matrix.
				Eigen::MatrixXd const& outcome = design.outcome() ;
				for( int i = 0; i < outcome.rows(); ++i ) {
					if(
						(outcome.row(i).sum() == outcome.row(i).sum())
						&& (
							( outcome.row(i).sum() != 1.0 )
							|| ( outcome(i,0) != 1.0 && outcome(i,1) != 1.0 )
						)
					) {
						std::cerr << "outcome:\n" << outcome << "\n" ;
						throw genfile::BadArgumentError(
							"metro::case_control::Logistic::check_design()",
							"design",
							( boost::format( "Outcome (%f for individual %d) has non-missing level other than zero or one." ) % outcome(i) % (i+1) ).str()
						) ;
					}
				}
			}
		}

		Logistic::Logistic(
			regression::Design::UniquePtr design
		):
			m_design( design.release() ),
			m_design_owned( true ),
			m_get_parameter_name( &standard_parameter_name ),
			m_parameters( Vector::Zero( m_design->matrix().cols() ) ),
			m_state( e_Uncomputed )
		{
			check_design( *m_design ) ;
		}

		Logistic::Logistic(
			regression::Design& design
		):
			m_design( &design ),
			m_design_owned( false ),
			m_get_parameter_name( &standard_parameter_name ),
			m_parameters( Vector::Zero( m_design->matrix().cols() ) ),
			m_state( e_Uncomputed )
		{
			assert( m_design != 0 ) ;
			check_design( *m_design ) ;
		}
		
		Logistic::~Logistic() {
			if( m_design_owned ) {
				delete m_design ;
			}
		}

		std::string Logistic::get_summary() const {
			using genfile::string_utils::to_string ;
			return "Logistic( "
				+ to_string( m_design->outcome().rows() )
				+ " samples )" ;
		}

		void Logistic::set_parameter_naming_scheme( Logistic::GetParameterName get_parameter_name ) {
			assert( get_parameter_name ) ;
			m_get_parameter_name = get_parameter_name ;		
		}
	
		std::string Logistic::get_parameter_name( std::size_t i ) const {
			IntegerMatrix const identity = identify_parameters() ;
			return m_get_parameter_name( m_design->get_predictor_name( identity(i,1) ), m_design->get_outcome_name(identity(i,0)) ) ;
		}
	
		LogLikelihood::IntegerMatrix Logistic::identify_parameters() const {
			IntegerMatrix result = IntegerMatrix::Zero( m_design->matrix().cols(), 2 ) ;
			// Each parameter corresponds to an effect for the non-baseline outcome level
			// and a column of the design matrix.
			result.col(0).setConstant( 1.0 ) ;
			// There is one parameter per design matrix column.
			for( int i = 0; i < m_design->matrix().cols(); ++i ) {
				result( i, 1 ) = i ;
			}
			return result ;
		}

		int Logistic::number_of_outcomes() const {
			return 2 ;
		}

		void Logistic::set_predictor_levels(
			Matrix const& levels,
			Matrix const& probabilities,
			std::vector< metro::SampleRange > const& included_samples
		) {
			m_design->set_predictor_levels( levels, probabilities, included_samples ) ;
			m_state = e_Uncomputed ;
#if DEBUG_LOGLIKELIHOOD
			std::cerr << "Logistic::set_predictor_probs(): Included samples:" ;
			for( std::size_t i = 0; i < included_samples.size(); ++i ) {
				std::cerr << " " << included_samples[i].begin() << "-" << included_samples[i].end() ;
			}
			std::cerr << "\n" ;
#endif
		}

		void Logistic::evaluate_at( Point const& parameters, int const numberOfDerivatives ) {
			evaluate_at_impl( parameters, m_design->per_predictor_included_samples(), numberOfDerivatives ) ;
		}

		void Logistic::evaluate_at_impl(
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
			int const N = m_design->outcome().rows() ;
			std::cerr << "==== Logistic::evaluate_at() ====\n" ;
			std::cerr << "  m_state: " << ( boost::format( "0x%x" ) % m_state ).str() << ", "
				<< "Number of samples: " << N << ", "
				<< "Included samples:" ;
			for( std::size_t i = 0; i < included_samples.size(); ++i ) {
				std::cerr << " " << included_samples[i].begin() << "-" << included_samples[i].end() ;
			}
			std::cerr << "\n" ;
			std::cerr << "  Params: " << parameters.transpose() << ".\n" ;
			std::cerr << "  Phenotypes: " << m_design->outcome().col(1).head( std::min( N, 5 ) ).transpose() << "...\n";
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
		
		Logistic::Point const& Logistic::get_parameters() const {
			return m_parameters ;
		}

		// Calculate P( outcome | genotype, covariates, parameters ).
		void Logistic::calculate_outcome_probabilities(
			Vector const& parameters,
			Matrix const& outcome,
			Logistic::Matrix* result
		) const {
			assert( outcome.cols() == 2 ) ;
			result->resize( outcome.rows(), m_design->get_number_of_predictor_levels() ) ;
			Vector linear_combination ;
			int const number_of_predictor_levels = m_design->get_number_of_predictor_levels() ;
			for( int g = 0; g < number_of_predictor_levels; ++g ) {
				Matrix const& design_matrix = m_design->set_predictor_level( g ).matrix() ;
				result->col( g ) = evaluate_mean_function( design_matrix * parameters, outcome ) ;
			}
		}
		
		Logistic::Vector Logistic::evaluate_mean_function( Vector const& linear_combinations, Matrix const& outcomes ) const {
			assert( linear_combinations.size() == outcomes.rows() ) ;
			// Only the 2nd column of the outcome (i.e. outcome is not baseline) is used in this calculation.
			Vector exps = linear_combinations.array().exp() ;
			Vector ones = Vector::Ones( linear_combinations.size() ) ;
			return ( ones.array() + outcomes.col(1).array() * ( exps.array() - ones.array() ) )  / ( ones + exps ).array() ;
		}

		void Logistic::compute_value_of_function( Matrix const& V, std::vector< metro::SampleRange > const& included_samples ) {
			m_value_of_function = 0.0 ;
			for( std::size_t i = 0; i < included_samples.size(); ++i ) {
				int const start_row = included_samples[i].begin() ;
				int const end_row = included_samples[i].end() ;
				m_value_of_function += (
					V.block( start_row, 0, end_row - start_row, V.cols() ).rowwise().sum().array().log().sum()
				) ;
			}
		}

		void Logistic::compute_value_of_first_derivative( Matrix const& A, std::vector< metro::SampleRange > const& included_samples, Matrix* B ) {
			int const N = m_design->outcome().rows() ;
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
				Vector signs = (( m_design->outcome().col(1) * 2.0 ) - Vector::Ones( N ) ) ; // == (-1)^{phenotype + 1}
				for( int g = 0; g < N_predictor_levels; ++g ) {
					Matrix const& design_matrix = m_design->set_predictor_level( g ).matrix() ;
#if DEBUG_LOGLIKELIHOOD
					std::cerr << "Logistic::compute_value_of_first_derivative(): design_matrix for predictor level " << g << " =\n"
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

		void Logistic::compute_value_of_second_derivative( Matrix const& B, std::vector< metro::SampleRange > const& included_samples ) {
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
		
		double Logistic::get_value_of_function() const {
			return m_value_of_function ;
		}

		Logistic::Vector Logistic::get_value_of_first_derivative() const {
			return m_value_of_first_derivative ;
		}

		Logistic::Matrix Logistic::get_value_of_second_derivative() const {
			return m_value_of_second_derivative ;
		}
	}
}
