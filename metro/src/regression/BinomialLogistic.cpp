
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
#include "metro/regression/BinomialLogistic.hpp"
#include "metro/intersect_ranges.hpp"
#include "metro/summation.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/Error.hpp"

// #define DEBUG_LOGLIKELIHOOD 1

namespace metro {
	namespace regression {
		namespace {
			std::string standard_parameter_name( std::string const& predictor_name, std::string const& outcome_name ) {
				return predictor_name + "/" + outcome_name ;
			}
		}
		
		BinomialLogistic::UniquePtr BinomialLogistic::create(
			regression::Design& design
		) {
			return BinomialLogistic::UniquePtr(
				new BinomialLogistic( design )
			) ;
		}

		BinomialLogistic::UniquePtr BinomialLogistic::create(
			regression::Design::UniquePtr design
		) {
			return BinomialLogistic::UniquePtr(
				new BinomialLogistic( design )
			) ;
		}

		namespace {
			void check_design( regression::Design const& design ) {
				// Verify design looks sane.
				BinomialLogistic::Matrix const& outcome = design.outcome() ;
				// We insist that outcomes are smallish non-negative integers.
				// Check this by comparing to absolute valued, floored version.
				if( outcome.rows() > 0 && outcome.cols() > 0
					&& (outcome.array() - outcome.array().abs().floor()).abs().maxCoeff() != 0 ) {
					throw genfile::BadArgumentError(
						"metro::case_control::BinomialLogistic::check_design()",
						"design",
						"Outcome does not consist solely of nonnegative integers."
					) ;
				}
			}
		}

		BinomialLogistic::BinomialLogistic(
			regression::Design::UniquePtr design
		):
			m_design( design.release() ),
			m_design_owned( true ),
			m_get_parameter_name( &standard_parameter_name ),
			m_parameters( Vector::Zero( m_design->matrix().cols() ) ),
			m_included_samples( 1, SampleRange( 0, m_design->matrix().rows() )),
			m_numberOfDerivativesComputed( -1 ),
			m_numberOfCDLDerivativesComputed( -1 )
		{
			check_design( *m_design ) ;
		}

		BinomialLogistic::BinomialLogistic(
			regression::Design::UniquePtr design,
			std::vector< metro::SampleRange > included_samples
		):
			m_design( design.release() ),
			m_design_owned( true ),
			m_get_parameter_name( &standard_parameter_name ),
			m_parameters( Vector::Zero( m_design->matrix().cols() ) ),
			m_included_samples( included_samples ),
			m_numberOfDerivativesComputed( -1 ),
			m_numberOfCDLDerivativesComputed( -1 )
		{
			check_design( *m_design ) ;
		}

		BinomialLogistic::BinomialLogistic(
			regression::Design& design
		):
			m_design( &design ),
			m_design_owned( false ),
			m_get_parameter_name( &standard_parameter_name ),
			m_parameters( Vector::Zero( m_design->matrix().cols() ) ),
			m_included_samples( 1, SampleRange( 0, m_design->matrix().rows() )),
			m_numberOfDerivativesComputed( -1 ),
			m_numberOfCDLDerivativesComputed( -1 )
		{
			assert( m_design != 0 ) ;
			check_design( *m_design ) ;
		}

		BinomialLogistic::BinomialLogistic(
			regression::Design& design,
			std::vector< metro::SampleRange > included_samples
		):
			m_design( &design ),
			m_design_owned( false ),
			m_get_parameter_name( &standard_parameter_name ),
			m_parameters( Vector::Zero( m_design->matrix().cols() ) ),
			m_included_samples( included_samples ),
			m_numberOfDerivativesComputed( -1 ),
			m_numberOfCDLDerivativesComputed( -1 )
		{
			assert( m_design != 0 ) ;
			check_design( *m_design ) ;
		}
		
		BinomialLogistic::~BinomialLogistic() {
			if( m_design_owned ) {
				delete m_design ;
			}
		}

		std::string BinomialLogistic::get_summary() const {
			using genfile::string_utils::to_string ;
			std::string result = "BinomialLogistic( "
				+ to_string( m_design->outcome().rows() )
				+ " samples ): " ;
			result += "(" + m_design->get_outcome_name(1) + ") ~ " ;
			IntegerMatrix const identity = identify_parameters() ;
			for( int i = 0; i < identity.rows(); ++i ) {
				result +=
					((i>0) ? " + " : "")
					+ get_parameter_name(i) ;
			}
			return result ;
		}

		void BinomialLogistic::set_parameter_naming_scheme( BinomialLogistic::GetParameterName get_parameter_name ) {
			assert( get_parameter_name ) ;
			m_get_parameter_name = get_parameter_name ;		
		}
	
		std::string BinomialLogistic::get_parameter_name( std::size_t i ) const {
			IntegerMatrix const identity = identify_parameters() ;
			return m_get_parameter_name( m_design->get_predictor_name( identity(i,1) ), m_design->get_outcome_name(identity(i,0)) ) ;
		}
	
		LogLikelihood::IntegerMatrix BinomialLogistic::identify_parameters() const {
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

		int BinomialLogistic::number_of_outcomes() const {
			return 2 ;
		}

		void BinomialLogistic::evaluate_at( Point const& parameters, int const numberOfDerivatives ) {
			evaluate_at_impl(
				parameters,
				metro::impl::intersect_ranges(
					m_included_samples,
					m_design->nonmissing_samples()
				),
				numberOfDerivatives
			) ;
		}

		void BinomialLogistic::evaluate_at_impl(
			Point const& parameters,
			std::vector< metro::SampleRange > const& evaluation_samples,
			int const numberOfDerivatives
		) {
			assert( parameters.size() == m_design->matrix().cols() ) ;
			assert( numberOfDerivatives < 3 ) ;

			if( parameters != m_parameters || evaluation_samples != m_evaluated_samples ) {
				m_numberOfDerivativesComputed = -1 ;
				m_numberOfCDLDerivativesComputed = -1 ;
				m_parameters = parameters ;
				m_evaluated_samples = evaluation_samples ;
			}
			
#if DEBUG_LOGLIKELIHOOD
			std::cerr << std::fixed << std::setprecision(4) ;
			int const N = m_design->outcome().rows() ;
			std::cerr << "==== BinomialLogistic::evaluate_at() ====\n" ;
			std::cerr << "  m_numberOfDerivativesComputed: " << m_numberOfDerivativesComputed << ", "
				<< "Number of samples: " << N << ", "
				<< "Included samples:" ;
			for( std::size_t i = 0; i < evaluation_samples.size(); ++i ) {
				std::cerr << " " << evaluation_samples[i].begin() << "-" << evaluation_samples[i].end() ;
			}
			std::cerr << "\n" ;
			std::cerr << "   Params: " << parameters.transpose() << ".\n" ;
			std::cerr << "   Phenotypes: " << m_design->outcome().col(1).segment(0, std::min( N, 5 )).transpose() << "...\n";
			std::cerr << "   Design:\n" << m_design->matrix().block( 0, 0, std::min( N, 5 ), m_design->matrix().cols() ) << "...\n";
#endif
			
			evaluate( numberOfDerivatives ) ;
		}
		
		void BinomialLogistic::evaluate( int const numberOfDerivatives ) {
			bool computeFunction = ( m_numberOfDerivativesComputed < 0 ) ;
			bool compute1stDerivative = ( (numberOfDerivatives > 0) & (m_numberOfDerivativesComputed < 1) ) ;
			bool compute2ndDerivative = ( (numberOfDerivatives > 1) & (m_numberOfDerivativesComputed < 2) ) ;

#if DEBUG_LOGLIKELIHOOD
			std::cerr << "BinomialLogistic::evaluate(): numberOfDerivatives = " << numberOfDerivatives << ".\n" ;
			std::cerr << "BinomialLogistic::evaluate(): m_numberOfDerivativesComputed = " << m_numberOfDerivativesComputed << ".\n" ;
			std::cerr << "Computing: "
				<< (computeFunction ? "function " : "" )
				<< (compute1stDerivative ? "1st derivative " : "" )
				<< (compute2ndDerivative ? "2nd derivative " : "" )
				<< "\n" ;
			//std::cerr << "  Predictor levels: " << m_design->get_predictor_levels().transpose() << ".\n" ;
#endif

			if( numberOfDerivatives > m_numberOfCDLDerivativesComputed ) {
#if DEBUG_LOGLIKELIHOOD
				std::cerr << "BinomialLogistic::evaluate(): computing CDL\n" ;
#endif
				// TODO: compute fx, dfx, ddfx here.
				compute_complete_data_likelihood_and_derivatives(
					m_parameters,
					&m_hx,
					&m_normalisedDhx,
					&m_normalisedDdhx,
					m_evaluated_samples,
					2
				) ;
				// Convert from complete data likelihood into expected complete data likelihood
				// and derivatives
				m_hx.array() *= m_design->get_predictor_level_probabilities().array() ;
				m_normalisedDhx.array() *= m_design->get_predictor_level_probabilities().array() ;
				m_normalisedDdhx.array() *= m_design->get_predictor_level_probabilities().array() ;
		
				// normalise the derivatives by h, as required in deritative computations
				typedef Eigen::DiagonalMatrix< double, Eigen::Dynamic > DiagonalMatrix ;
				DiagonalMatrix one_over_h = (m_hx.rowwise().sum().asDiagonal().inverse()) ;
				m_normalisedDhx = one_over_h * m_normalisedDhx ;
				m_normalisedDdhx = one_over_h * m_normalisedDdhx ;

				m_numberOfCDLDerivativesComputed = 2 ;
			}
#if DEBUG_LOGLIKELIHOOD
			else {
				std::cerr << "BinomialLogistic::evaluate(): skipping CDL computation.\n" ;
			}
#endif
			
			if( computeFunction ) {
#if DEBUG_LOGLIKELIHOOD
				std::cerr << "BinomialLogistic::evaluate(): computing fx\n" ;
#endif
				compute_value_of_loglikelihood(
					m_hx,
					&m_value_of_function,
					m_evaluated_samples
				) ;
				m_numberOfDerivativesComputed = 0 ;
			}

			if( compute1stDerivative ) {
#if DEBUG_LOGLIKELIHOOD
				std::cerr << "BinomialLogistic::evaluate(): computing dfx\n" ;
#endif
				compute_value_of_first_derivative(
					m_normalisedDhx,
					&m_first_derivative_terms,
					&m_value_of_first_derivative,
					m_evaluated_samples
				) ;
				m_numberOfDerivativesComputed = 1 ;
			}

			if( compute2ndDerivative ) {
#if DEBUG_LOGLIKELIHOOD
				std::cerr << "BinomialLogistic::evaluate(): computing ddfx\n" ;
#endif
				compute_value_of_second_derivative(
					m_first_derivative_terms,
					m_normalisedDdhx,
					&m_value_of_second_derivative,
					m_evaluated_samples
				) ;
				m_numberOfDerivativesComputed = 2 ;
			}

#if DEBUG_LOGLIKELIHOOD
			std::cerr << "  " << ( computeFunction ? "(computed) " : "(memoized) " ) << "function value = " << m_value_of_function << "\n" ;
#endif

#if DEBUG_LOGLIKELIHOOD
			if( m_numberOfDerivativesComputed > 0 ) {
				std::cerr << "  " << ( compute1stDerivative ? "(computed) " : "(memoized) " ) << "derivative = " << m_value_of_first_derivative.transpose() << "\n" ;
			} else {
				std::cerr << "  (1st derivative not computed).\n" ;
			}
#endif

#if DEBUG_LOGLIKELIHOOD
			if( m_numberOfDerivativesComputed > 1 ) {
				std::cerr << "  " << ( compute2ndDerivative ? "(computed) " : "(memoized) " ) << "2nd derivative = \n"
					<< m_value_of_second_derivative << ".\n" ;
			} else {
				std::cerr << "  (2nd derivative not computed).\n" ;
			}
			std::cerr << "====\n" ;
#endif
		}
		
		void BinomialLogistic::compute_complete_data_likelihood_and_derivatives(
			Point const& parameters,
			Matrix* fx,
			Matrix* dfx,
			Matrix* ddfx,
			std::vector< metro::SampleRange > const& included_samples,
			int const numberOfDerivatives
		) {
			int const N = m_design->outcome().rows() ;
			int const D = m_design->matrix().cols() ;
			int const L = m_design->get_number_of_predictor_levels() ;
			Matrix const& outcome = m_design->outcome() ;

			fx->resize( N, L ) ;
			if( numberOfDerivatives > 0 ) {
				dfx->resize( N, L ) ;
				if( numberOfDerivatives > 1 ) {
					ddfx->resize( N, L ) ;
				}
			}

			// Outcome data is two columns, conceptually row sum is n
			// (the number of trials) and 2nd column is k (number of successes / non-baseline outcomes)
			// The computation requires us to compute f0 and f1, the prob of a single
			// baseline / non-baseline observation.
			assert( outcome.cols() == 2 ) ;
			m_f1.resize( N ) ;
			Vector linear_combination ;
			auto one = Vector::Ones(N) ;
			for( int g = 0; g < L; ++g ) {
				Matrix const& design_matrix = m_design->set_predictor_level( g ).matrix() ;
#if DEBUG_LOGLIKELIHOOD
				std::cerr << "design(" << g << ") =\n"
					<< design_matrix.block( 0, 0, std::min( 10, int(design_matrix.rows())), design_matrix.cols() ) << ".\n" ;
				std::cerr << "outcome(" << g << ") =\n"
					<< outcome.block( 0, 0, std::min( 10, int(outcome.rows())), outcome.cols() ) << ".\n" ;
#endif
				m_f1 = (design_matrix * parameters).array().exp() ;
				m_f1 = m_f1.array() * ((one + m_f1).array().inverse()) ;
				
				// Likelihood is obtained by raising f1 and f0 to appropriate powers
				
				
#if DEBUG_LOGLIKELIHOOD
				std::cerr << "(*fx) = \n"
					<< fx->block( 0, 0, std::min( 10, int(fx->rows())), fx->cols() ) << ".\n" ;
				std::cerr << "f1 = \n"
					<< m_f1.block( 0, 0, std::min( 10, int(m_f1.rows())), m_f1.cols() ) << ".\n" ;
#endif
				(*fx).col(g)
					= m_f1.array().pow( outcome.col(1).array() )
					* (one.array() - m_f1.array() ).pow( outcome.col(0).array() ) ;

#if DEBUG_LOGLIKELIHOOD
				std::cerr << "Computed f1 and fx, f1 = " << m_f1.segment( 0, std::min( 5, N )).transpose() << ", \n"
					<< "fx(" << g << ") = " << fx->col(g).segment(0, std::min( 5, N )).transpose() << ".\n" ;
#endif

				if( numberOfDerivatives > 0 ) {
					// compute coefficient of x^t in 1st derivative
					// To share computation with 2nd derivative,
					// first compute the shared term, we multiple by fx below.

#if DEBUG_LOGLIKELIHOOD
					std::cerr << "    dfx: " << dfx->rows() << "x" << dfx->cols() << ".\n" ;
					std::cerr << "outcome: " << outcome.rows() << "x" << outcome.cols() << ".\n" ;
#endif
					dfx->col(g)
						= ( outcome.col(1).array() * ( one - m_f1 ).array() - outcome.col(0).array() * m_f1.array() )
					;

#if DEBUG_LOGLIKELIHOOD
					std::cerr << "dfx(" << g << ") before fx term = "
						<< dfx->col(g).segment(0, std::min( 5, N )).transpose() << ".\n" ;
#endif

					if( numberOfDerivatives > 1 ) {
						// compute coefficient of x ⊗ x^t in 2nd derivative
						ddfx->col(g)
							= dfx->col(g).array().square()
							- (
									(outcome.col(0) + outcome.col(1)).array()
									* (m_f1.array() * ( one.array() - m_f1.array() ))
							).array()
						;
						ddfx->col(g).array() *= fx->col(g).array() ;

#if DEBUG_LOGLIKELIHOOD
						std::cerr << "ddfx(" << g << ") = "
							<< ddfx->col(g).segment(0, std::min( 5, N )).transpose() << ".\n" ;
#endif
					}
					dfx->col(g).array() *= fx->col(g).array() ;

#if DEBUG_LOGLIKELIHOOD
				std::cerr << "dfx(" << g << ") = "
					<< dfx->col(g).segment(0, std::min( 5, N )).transpose() << ".\n" ;
#endif
				}
			}
		}

		
		void BinomialLogistic::compute_value_of_loglikelihood(
			Matrix const& hx,
			double* result,
			std::vector< metro::SampleRange > const& included_samples
		) {
			(*result) = 0 ;
			for( std::size_t i = 0; i < included_samples.size(); ++i ) {
				int const start_row = included_samples[i].begin() ;
				int const end_row = included_samples[i].end() ;
				// Compute the log-likelihood for this block of samples as:
				// the sum over samples of...
				// (*result) += naive_sum(
				(*result) += neumaier_sum(
					// the likelihood terms for these samples...
					hx.block( start_row, 0, end_row - start_row, hx.cols() )
					 // ...summed over predictor levels...
					.rowwise().sum()
					// ...and logarithm-ed...
					.array().log()
				) ;
			}
		}

		void BinomialLogistic::compute_value_of_first_derivative(
			Matrix const& normalisedDhx,
			Matrix* result_terms,
			Vector* result,
			std::vector< metro::SampleRange > const& included_samples
		) {
			int const N = m_design->outcome().rows() ;
			int const D = m_design->matrix().cols() ;
			int const N_predictor_levels = m_design->get_number_of_predictor_levels() ;

			// Currently we assume
			// dimension of design matrix = dimension of parameters = dimension of derivative
			// More generally we may have more parameters, e.g. multinomial logistic or normal linear model
			// in which case something different must be done.  (Probably, identify_parameters can be
			// used to generalise to this case?)
			assert( D == m_parameters.size() ) ;

			result->setZero(D) ;
			result_terms->setZero( N, D ) ;

			for( int g = 0; g < N_predictor_levels; ++g ) {
				Matrix const& design_matrix = m_design->set_predictor_level( g ).matrix() ;
				for( std::size_t i = 0; i < included_samples.size(); ++i ) {
					int const start_row = included_samples[i].begin() ;
					int const end_row = included_samples[i].end() ;

					// Eigen expressions for relevant blocks of matrices
					MatrixBlock termsBlock = result_terms->block(
						start_row, 0,
						end_row - start_row, result_terms->cols()
					) ;
					
					ConstMatrixBlock normalisedDhxBlock = normalisedDhx.block(
						start_row, 0,
						end_row - start_row, normalisedDhx.cols()
					) ;

					ConstMatrixBlock design_matrix_block = design_matrix.block(
						start_row, 0,
						end_row - start_row, design_matrix.cols()
					) ;

					// Compute the contribution of these samples & this predictor level
					termsBlock +=
						normalisedDhxBlock.col(g)
						.asDiagonal()
						* design_matrix_block
					;
					
				}
			}
			// Now the rows of terms are the contribution of each sample to the derivative
			(*result) = result_terms->colwise().sum() ;
		}
		
		void BinomialLogistic::compute_value_of_second_derivative(
			Matrix const& first_derivative_terms,
			Matrix const& normalisedDdhx,
			Matrix* result,
			std::vector< metro::SampleRange > const& included_samples
		) {
			int const N = m_design->outcome().rows() ;
			int const D = m_design->matrix().cols() ;
			int const N_predictor_levels = m_design->get_number_of_predictor_levels() ;

			// Currently we assume
			// dimension of design matrix = dimension of parameters = dimension of derivative
			// More generally we may have more parameters, e.g. multinomial logistic or normal linear model
			// in which case something different must be done.  (Probably, identify_parameters can be
			// used to generalise to this case?)
			assert( D == m_parameters.size() ) ;

			result->setZero(D,D) ;

			for( int g1 = 0; g1 < N_predictor_levels; ++g1 ) {
				Matrix const& design_matrix = m_design->set_predictor_level( g1 ).matrix() ;
				for( std::size_t sample_i = 0; sample_i < included_samples.size(); ++sample_i ) {
					int const start_row = included_samples[sample_i].begin() ;
					int const end_row = included_samples[sample_i].end() ;

					ConstMatrixBlock design_matrix_block
						= design_matrix.block( start_row, 0, end_row - start_row, design_matrix.cols() ) ;

					// We are implementing the formula: sum_i (sum_x DD^t h_i)/(sum_x h_i) × x ⊗ x^t
					// Or, written in the way we compute it here, sum_g ( sum_i Z_i × x_i(g) ⊗ x_i(g)^t )
					// where Z_i = (DD^t h_i(g)) / h_i and x_i(g)^t is the design matrix row for sample i given predictor level g.
					// We let matrix multiplication do the sum over i for us, and compute this as
					// diag(Z) × X(g)^t ⊗ X(g)
					(*result) += (
							normalisedDdhx.col(g1).segment( start_row, end_row - start_row ).asDiagonal()
							* design_matrix_block
						)
						.transpose()
						* design_matrix_block
					;
				}
			}
			(*result) -= ( first_derivative_terms.transpose() * first_derivative_terms ) ;

#if DEBUG_LOGLIKELIHOOD
			std::cerr << "D = " << D << "\n" ;
			std::cerr << "first_derivative_terms: " << first_derivative_terms.rows() << "x" << first_derivative_terms.cols() << "\n" ;
			std::cerr << "result: " << result->rows() << "x" << result->cols() << "\n" ;
#endif
		}
		
		BinomialLogistic::Point const& BinomialLogistic::parameters() const {
			return m_parameters ;
		}

		double BinomialLogistic::get_value_of_function() const {
			return m_value_of_function ;
		}

		BinomialLogistic::Vector BinomialLogistic::get_value_of_first_derivative() const {
			return m_value_of_first_derivative ;
		}

		BinomialLogistic::Matrix BinomialLogistic::get_value_of_second_derivative() const {
			return m_value_of_second_derivative ;
		}
	}
}
