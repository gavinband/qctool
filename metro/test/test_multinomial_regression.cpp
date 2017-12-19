#include <vector>
#include <iostream>
#include "test_case.hpp"
#include "metro/SampleRange.hpp"
#include "metro/case_control/RegressionDesign.hpp"
#include "metro/case_control/LogisticRegressionLogLikelihood.hpp"
#include "metro/case_control/MultinomialRegressionLogLikelihood.hpp"

// #define DEBUG_TESTS 1
namespace {
	void check_close_or_small(
	    double expected, double observed, 
	    double small, double pct_tol
	) {
	    if (std::fabs(expected) < small) {
	        BOOST_CHECK_SMALL(observed, small);
	    } else {
	        BOOST_CHECK_CLOSE(expected, observed, pct_tol);
	    }
	}
}

BOOST_AUTO_TEST_SUITE( test_multinomial ) ;

AUTO_TEST_CASE( test_multinomialregression_one_sample ) {
	using metro::SampleRange ;
	using namespace metro::case_control ;
	typedef RegressionDesign::Matrix Matrix ;
	typedef RegressionDesign::Vector Vector ;
	using std::exp ;
	using std::log ;
	
	int const nSamples = 1 ;
	std::vector< metro::SampleRange > included_samples( 1, metro::SampleRange( 0, 1 ) ) ;
	Matrix predictor_levels = Matrix::Zero( 2, 1 ) ;
	predictor_levels(0,0) = 1 ;
	predictor_levels(1,0) = 2 ;

	for( int nOutcomes = 3; nOutcomes < 5; ++nOutcomes ) {
		Vector outcome = Vector::Zero(nSamples) ;
		Vector outcome_nonmissingness = Vector::Constant( nSamples, 1.0 ) ;
		Matrix covariates = Matrix::Zero( nSamples, 0 ) ;
		Matrix covariate_nonmissingness = Matrix::Constant( nSamples, 0, 1.0 ) ;

//		for( int outcome_i = 0; outcome_i < nOutcomes; ++outcome_i ) {
		for( int outcome_i = 1; outcome_i < nOutcomes; ++outcome_i ) {
			outcome(0) = outcome_i ;
			
			MultinomialRegressionLogLikelihood ll(
				RegressionDesign::create(
					outcome, outcome_nonmissingness, "case",
					covariates, covariate_nonmissingness, std::vector< std::string >(),
					std::vector< std::string >( 1, "predictor" )
				),
				nOutcomes
			) ;
			
			// Set up a single predictor, equal to 1 with certainty.
			Matrix predictor_probabilities = Matrix::Zero( 1, 2 ) ;
			predictor_probabilities(0, 0) = 1.0 ;
			ll.set_predictor_levels( predictor_levels, predictor_probabilities, included_samples ) ;

			int const nParameters = 2 * ( nOutcomes - 1 ) ;
			int const M = nOutcomes - 1 ;
			// Now test values.
			{
				Vector parameters = Vector::Zero( nParameters ) ;
				ll.evaluate_at( parameters ) ;
				// All probs are equal to 1/(M+1).  So outcome probability equals this.
				Matrix parameter_matrix = parameters ;
				parameter_matrix.resize( 2, nOutcomes - 1 ) ;
#if DEBUG_TESTS
				std::cerr << "parameters = " << parameters.transpose() << ".\n" ;
				std::cerr << "outcome = " << outcome.transpose() << ".\n" ;
#endif
			
				double const expected_ll = std::log( 1.0 / double( nOutcomes ) ) ;
				BOOST_CHECK_CLOSE( expected_ll, ll.get_value_of_function(), 0.00001 ) ;
				
				Matrix F = Matrix::Constant( 1, nOutcomes, 1.0 / double( nOutcomes ) ) ;
				
				Matrix B = -F.block( 0, 1, 1, M ) ;
				if( outcome_i > 0 ) {
					B( 0, outcome_i - 1 ) += 1.0 ;
				}
				
				Matrix design_matrix = Matrix::Constant( 1, 2, 1.0 ) ;
				Matrix expected_first_derivative = Matrix::Zero( 1, M * 2 ) ;
				for( int j = 0; j < M ; ++j ) {
					expected_first_derivative.block( 0, j*2, 1, 2 ) = design_matrix.row(0) * B( 0, j ) ;
				}

				Vector const first_derivative = ll.get_value_of_first_derivative() ;
				for( int j = 0; j < M*2; ++j ) {
					BOOST_CHECK_CLOSE( expected_first_derivative(0,j), first_derivative(j), 0.00001 ) ;
				}

				Matrix C = Matrix::Zero( M, M ) ;
				for( int j = 0; j < M ; ++j ) {
					double Z_ij = ( j == outcome_i - 1 ) ? ( 1 - F(j+1) ) : ( - F(j+1) ) ;
					for( int k = 0; k < M; ++k ) {
						double Z_ik = ( k == outcome_i - 1 ) ? ( 1 - F(k+1) ) : ( - F(k+1) ) ;
						double Z_jk = ( j == k ) ? ( 1 - F(k+1) ) : ( - F(k+1) ) ;
//						std::cerr << "j:" << j << ", k:" << k << ", Zij:" << Z_ij << ", Z_ik:" << Z_ik << ", Z_jk:" << Z_jk << ".\n" ;
						C( j, k ) = ( Z_ij * Z_ik ) - ( F( j+1 ) * Z_jk ) ;
					}
				}

				Matrix expected_2nd_derivative = Matrix::Zero( M*2, M*2 ) ;
				for( int j = 0; j < M ; ++j ) {
					for( int k = 0; k < M ; ++k ) {
						expected_2nd_derivative.block( j*2, k*2, 2, 2 )
							= ( C( j, k ) - B( 0, j ) * B( 0, k ) ) * design_matrix.row(0).transpose() * design_matrix.row(0) ;
					}
				}

				Matrix const second_derivative = ll.get_value_of_second_derivative() ;
				for( int j = 0; j < M*2; ++j ) {
					for( int k = 0; k < M*2; ++k ) {
						BOOST_CHECK_CLOSE( expected_2nd_derivative(j,k), second_derivative(j,k), 0.00001 ) ;
					}
				}
				
#if DEBUG_TESTS
				std::cerr << "test_multinomialregression_one_sample:\n" ;
				std::cerr << "nOutcomes = " << nOutcomes << "\n" ;
				std::cerr << "F = " << F << "\n" ;
				std::cerr << "B = " << B << "\n" ;
				std::cerr << "C = \n" << C << "\n" ;
				std::cerr << "expected 2nd derivative = \n"
					<< expected_2nd_derivative << "\n" ;
				std::cerr << "actual 2nd derivative = \n"
					<< second_derivative << "\n" ;
#endif
			}

			for( int nonzero_i = 0; nonzero_i < nParameters; ++nonzero_i )
			{

				Vector parameters = Vector::Zero( 2 * ( nOutcomes - 1 ) ) ;
				parameters( nonzero_i ) = log(2.0) ;
#if DEBUG_TESTS
				std::cerr << "##############\n" ;
				std::cerr << "nOutcomes = " << nOutcomes << ".\n" ;
				std::cerr << "parameters = " << parameters.transpose() << ".\n" ;
				std::cerr << "outcome = " << outcome.transpose() << ".\n" ;
				std::cerr << "nonzero = " << nonzero_i << ".\n" ;
#endif
				ll.evaluate_at( parameters ) ;

				// f_i will equal 1/(M+2), or 2/(M+2) for i=floor(nonzero_i/2).
				Matrix parameter_matrix = parameters ;
				parameter_matrix.resize( 2, nParameters ) ;
			
				// So if the outcome equals floor(nonzero_i/2) we have
				int const special_outcome = nonzero_i / 2 ;
				double expected_ll ;
				if( outcome_i == 0 ) {
					expected_ll = std::log( 1.0 / ( M + 2.0 ) ) ;
				} else if( ( outcome_i - 1 ) == special_outcome ) {
					expected_ll = std::log( 2.0 / ( M + 2.0 ) ) ;
				} else {
					expected_ll = std::log( 1.0 / ( M + 2.0 ) ) ;
				}
				BOOST_CHECK_CLOSE( expected_ll, ll.get_value_of_function(), 0.00001 ) ;

				Matrix F = Matrix::Constant( 1, nOutcomes, 1.0 / ( M + 2.0 ) ) ;
				F( special_outcome + 1 ) = 2.0 / ( M + 2.0 ) ;
				
				Matrix B = -F.block( 0, 1, 1, M ) ;
				if( outcome_i > 0 ) {
					B( 0, outcome_i - 1 ) += 1.0 ;
				}
				Matrix design_matrix = Matrix::Constant( 1, 2, 1.0 ) ;
				Matrix expected_first_derivative = Matrix::Zero( 1, M * 2 ) ;
				for( int j = 0; j < M ; ++j ) {
					expected_first_derivative.block( 0, j*2, 1, 2 ) = design_matrix.row(0) * B( 0, j ) ;
				}

				Vector const first_derivative = ll.get_value_of_first_derivative() ;
				for( int j = 0; j < M*2; ++j ) {
					BOOST_CHECK_CLOSE( expected_first_derivative(0,j), first_derivative(j), 0.00001 ) ;
				}

				Matrix C = Matrix::Zero( M, M ) ;
				for( int j = 0; j < M ; ++j ) {
					double Z_ij = ( j == outcome_i - 1 ) ? ( 1 - F(j+1) ) : ( - F(j+1) ) ;
					for( int k = 0; k < M; ++k ) {
						double Z_ik = ( k == outcome_i - 1 ) ? ( 1 - F(k+1) ) : ( - F(k+1) ) ;
						double Z_jk = ( j == k ) ? ( 1 - F(k+1) ) : ( - F(k+1) ) ;
						C( j, k ) = ( Z_ij * Z_ik ) - ( F( j+1 ) * Z_jk ) ;
						//std::cerr << "j:" << j << ", k:" << k << ", Zij:" << Z_ij << ", Z_ik:" << Z_ik << ", Z_jk:" << Z_jk << ".\n" ;
					}
				}
				Matrix expected_2nd_derivative = Matrix::Zero( M*2, M*2 ) ;
				for( int j = 0; j < M ; ++j ) {
					for( int k = 0; k < M ; ++k ) {
						expected_2nd_derivative.block( j*2, k*2, 2, 2 )
							= ( C( j, k ) - B( 0, j ) * B( 0, k ) ) * design_matrix.row(0).transpose() * design_matrix.row(0) ;
					}
				}

				Matrix const second_derivative = ll.get_value_of_second_derivative() ;
				
#if DEBUG_TESTS
				std::cerr << "test_multinomialregression_one_sample:\n" ;
				std::cerr << "nOutcomes = " << nOutcomes << "\n" ;
				std::cerr << "F = " << F << "\n" ;
				std::cerr << "B = " << B << "\n" ;
				std::cerr << "C = \n" << C << "\n" ;
				std::cerr << "expected 2nd derivative = \n"
					<< expected_2nd_derivative << "\n" ;
				std::cerr << "actual 2nd derivative = \n"
					<< second_derivative << "\n" ;
#endif
				for( int j = 0; j < M*2; ++j ) {
					for( int k = 0; k < M*2; ++k ) {
						BOOST_CHECK_CLOSE( expected_2nd_derivative(j,k), second_derivative(j,k), 0.00001 ) ;
					}
				}
			}
		}
	}
}

AUTO_TEST_CASE( test_multinomialregression_two_outcomes_two_samples ) {
	using metro::SampleRange ;
	using namespace metro::case_control ;
	typedef RegressionDesign::Matrix Matrix ;
	typedef RegressionDesign::Vector Vector ;
	using std::exp ;
	using std::log ;
	
	{
		int const nSamples = 2 ;
		Vector outcome = Vector::Zero(nSamples) ;
		outcome(1) = 1 ;
		Vector outcome_nonmissingness = Vector::Constant( nSamples, 1.0 ) ;
		Matrix covariates = Matrix::Zero( nSamples, 0 ) ;
		Matrix covariate_nonmissingness = Matrix::Constant( nSamples, 0, 1.0 ) ;

		MultinomialRegressionLogLikelihood ll(
			RegressionDesign::create(
				outcome, outcome_nonmissingness, "outcome",
				covariates, covariate_nonmissingness, std::vector< std::string >(),
				std::vector< std::string >( 1, "predictor" )
			),
			2
		) ;
		
		Matrix predictor_levels = Matrix::Zero( 2, 1 ) ;
		predictor_levels(0,0) = 0 ;
		predictor_levels(1,0) = 1 ;

		// Set half the sample to certainly 0, half to certainly 1.
		Matrix predictor_probabilities = Matrix::Zero( nSamples, 2 ) ;
		predictor_probabilities( 0, 0 ) = 1.0 ;
		predictor_probabilities( 1, 1 ) = 1.0 ;

		std::vector< metro::SampleRange > const included_samples = std::vector< metro::SampleRange >( 1, metro::SampleRange( 0, nSamples ) ) ;
		
		ll.set_predictor_levels( predictor_levels, predictor_probabilities, included_samples ) ;
		
		{
			Vector parameters = Vector::Zero( 2 ) ;
			ll.evaluate_at( parameters ) ;
			// all probabilities are one-half.
			Vector f1( 2 ) ;
			f1 << 0.5, 0.5 ;
			Vector F = f1 ; // F=f because probs are 1/2

			double expected_ll = nSamples * std::log( 0.5 ) ;
			BOOST_CHECK_CLOSE( ll.get_value_of_function(), expected_ll, 0.0000001 ) ;

			Matrix design_matrix = Matrix::Zero( nSamples, 2 ) ;
			design_matrix.col(0).setConstant( 1.0 ) ;
			design_matrix( 1, 1 ) = 1.0 ;
			
			Vector B = F.array() * ( Vector::Constant( 2, 1.0 ) - F ).array() ;
			for( int i = 0; i < B.size(); ++i ) {
				B(i) = B(i) / F(i) ;
			}
			B(0) = -B(0) ;
			Vector const expected_first_derivative = ( B.transpose() * design_matrix ) ;
			Vector const first_derivative = ll.get_value_of_first_derivative() ;
			BOOST_CHECK_CLOSE( expected_first_derivative(0), first_derivative(0), 0.0000001 ) ;
			BOOST_CHECK_CLOSE( expected_first_derivative(1), first_derivative(1), 0.0000001 ) ;

			// Compute second derivative
			Vector C = ( Vector::Constant( 2, 1.0 ) - 2.0 * F ).array() * B.array() ;
#if DEBUG_TEST
			std::cerr << "expected C = " << C << ".\n" ;
#endif
			Matrix DD = Matrix::Zero( 2, 2 ) ;
			for( int i = 0; i < design_matrix.rows(); ++i ) {
				DD += ( C(i) - B(i) * B(i) ) * design_matrix.row(i).transpose() * design_matrix.row(i) ;
			}
#if DEBUG_TEST
			std::cerr << "expected second derivative = " << D << ".\n" ;
#endif
			Matrix const& second_derivative = ll.get_value_of_second_derivative() ;
			BOOST_CHECK_CLOSE( DD(0,0), second_derivative(0,0), 0.0000001 ) ;
			BOOST_CHECK_CLOSE( DD(0,1), second_derivative(0,1), 0.0000001 ) ;
			BOOST_CHECK_CLOSE( DD(1,0), second_derivative(1,0), 0.0000001 ) ;
			BOOST_CHECK_CLOSE( DD(1,1), second_derivative(1,1), 0.0000001 ) ;
		}

		{
			// Set a nonzero predictor effect parameter.
			Vector parameters = Vector::Zero( 2 ) ;
			parameters << 0.0, 1.0 ;
			ll.evaluate_at( parameters ) ;
			// Now we have exp( x^t theta ) = 1, e for the two samples.  Thus f1 is 1/2, e/1+e
			Vector f1( 2 ) ;
			f1 << 0.5, ( std::exp( 1 ) / ( 1.0 + std::exp(1.0) ) ) ;
			Vector F = f1 ; // f_0 = f_1 for first sample because the first prob is 0.5.
			double const expected_ll = F.array().log().sum() ;
			BOOST_CHECK_CLOSE( ll.get_value_of_function(), expected_ll, 0.0000001 ) ;

			// f1 = 1/2, e/1+e
			Matrix design_matrix = Matrix::Zero( nSamples, 2 ) ;
			design_matrix.col(0).setConstant( 1.0 ) ;
			design_matrix( 1, 1 ) = 1.0 ;
			
			Vector B = F.array() * ( Vector::Constant( 2, 1.0 ) - F ).array() ;
			for( int i = 0; i < B.size(); ++i ) {
				B(i) = B(i) / F(i) ;
			}
			B(0) = -B(0) ;
			Vector const expected_first_derivative = ( B.transpose() * design_matrix ) ;
			Vector const first_derivative = ll.get_value_of_first_derivative() ;
			BOOST_CHECK_CLOSE( expected_first_derivative(0), first_derivative(0), 0.0000001 ) ;
			BOOST_CHECK_CLOSE( expected_first_derivative(1), first_derivative(1), 0.0000001 ) ;

			// Compute second derivative
			Vector C = ( Vector::Constant( 2, 1.0 ) - 2.0 * F ).array() * B.array() ;
			
#if DEBUG_TESTS
			std::cerr << "test_multinomialregression_two_samples: B = " << B.transpose() << "\n" ;
			std::cerr << "test_multinomialregression_two_samples: C = " << C.transpose() << "\n" ;
#endif
			
			Matrix DD = Matrix::Zero( 2, 2 ) ;
			for( int i = 0; i < design_matrix.rows(); ++i ) {
				DD -= ( B(i) * B(i) ) * design_matrix.row(i).transpose() * design_matrix.row(i) ;
			}
#if DEBUG_TESTS
			std::cerr << "test_multinomialregression_two_samples: after B terms, expected 2nd derivative =\n" << DD << ".\n" ;
#endif
			for( int i = 0; i < design_matrix.rows(); ++i ) {
				DD += C(i) * design_matrix.row(i).transpose() * design_matrix.row(i) ;
			}
			Matrix const& second_derivative = ll.get_value_of_second_derivative() ;
#if DEBUG_TESTS
			std::cerr << "test_multinomialregression_two_samples: expected 2nd derivative =\n" << DD << ".\n" ;
#endif
			BOOST_CHECK_CLOSE( DD(0,0), second_derivative(0,0), 0.0000001 ) ;
			BOOST_CHECK_CLOSE( DD(0,1), second_derivative(0,1), 0.0000001 ) ;
			BOOST_CHECK_CLOSE( DD(1,0), second_derivative(1,0), 0.0000001 ) ;
			BOOST_CHECK_CLOSE( DD(1,1), second_derivative(1,1), 0.0000001 ) ;
		}
	}
}


AUTO_TEST_CASE( test_multinomialregression_two_outcomes_certain_predictors ) {
	using metro::SampleRange ;
	using namespace metro::case_control ;
	typedef RegressionDesign::Matrix Matrix ;
	typedef RegressionDesign::Vector Vector ;
	using std::exp ;
	using std::log ;
	
	std::vector< double > parameter_values ;
	parameter_values.push_back( -10 ) ;
	parameter_values.push_back( -3.141592654 ) ;
	parameter_values.push_back( -1 ) ;
	parameter_values.push_back( -0.25 ) ;
	parameter_values.push_back( 0 ) ;
	parameter_values.push_back( 0.25 ) ;
	parameter_values.push_back( 0.5 ) ;
	parameter_values.push_back( 1 ) ;
	parameter_values.push_back( 3.141592654 ) ;
	parameter_values.push_back( 10 ) ;
	
	for( int nCases = 0; nCases < 20; ++nCases ) {
		for( int nControls = 0; nControls < 20; ++nControls ) {
			int nSamples = nCases + nControls ;

			// contruct outcome
			Vector outcome = Vector::Zero(nSamples) ;
			for( int i = 0; i < nCases; ++i ) {
				outcome(i) = 1 ;
			}
			Vector outcome_nonmissingness = Vector::Constant( nSamples, 1.0 ) ;

			Matrix covariates = Matrix::Zero( nSamples, 0 ) ;
			Matrix covariate_nonmissingness = Matrix::Constant( nSamples, 0, 1.0 ) ;

			MultinomialRegressionLogLikelihood ll(
				RegressionDesign::create(
					outcome, outcome_nonmissingness, "outcome",
					covariates, covariate_nonmissingness, std::vector< std::string >(),
					std::vector< std::string >( 1, "predictor" )
				),
				2
			) ;
		
			// construct a 0/1 predictor.
			Matrix predictor_levels = Matrix::Zero( 2, 1 ) ;
			predictor_levels(0,0) = 0 ;
			predictor_levels(1,0) = 1 ;

			// Set half the sample to certainly 0, half to certainly 1.
			Matrix predictor_probabilities = Matrix::Zero( nSamples, 2 ) ;
			int const nOnes = nSamples / 2 ;
			for( int i = 0; i < nSamples; ++i ) {
				if( i < nOnes ) {
					predictor_probabilities( i, 1 ) = 1 ;
				} else {
					predictor_probabilities( i, 0 ) = 1 ;
				}
			}
			
			std::vector< metro::SampleRange > const included_samples = std::vector< metro::SampleRange >( 1, metro::SampleRange( 0, nSamples ) ) ;
		
			ll.set_predictor_levels( predictor_levels, predictor_probabilities, included_samples ) ;
		
			// Ok we have our likelihood.  Now compute with it.
			{
				ll.evaluate_at( Vector::Zero( 2 ) ) ;
				// all probabilities are one-half.
				double expected_ll = nSamples * std::log( 0.5 ) ;
				BOOST_CHECK_CLOSE( ll.get_value_of_function(), expected_ll, 0.0000001 ) ;

				Vector one_minus_f1 = ( outcome.array() - 0.5 ) ;
				Matrix design_matrix = Matrix::Zero( nSamples, 2 ) ;
				design_matrix.col(0).setConstant( 1.0 ) ;
				design_matrix.col(1).segment( 0, nOnes ).setConstant( 1.0 ) ;
				
				Vector const expected_first_derivative = ( one_minus_f1.transpose() * design_matrix ) ;
				Vector const first_derivative = ll.get_value_of_first_derivative() ;
				BOOST_CHECK_CLOSE( expected_first_derivative(0), first_derivative(0), 0.0000001 ) ;
				BOOST_CHECK_CLOSE( expected_first_derivative(1), first_derivative(1), 0.0000001 ) ;
				
				Matrix DD = Matrix::Zero( 2, 2 ) ;
				for( int sample = 0; sample < nSamples; ++sample ) {
					DD -= 0.25 * design_matrix.row( sample ).transpose() * design_matrix.row( sample ) ;
				}
				Matrix const& second_derivative = ll.get_value_of_second_derivative() ;
				for( int i = 0; i < 2; ++i ) {
					for( int j = 0; j < 2; ++j ) {
						BOOST_CHECK_CLOSE( DD(i,j), second_derivative(i,j), 0.000001 ) ;
					}
				}
			}

			Vector parameters = Vector::Zero( 2 ) ;
			for( std::size_t pi0 = 0; pi0 < parameter_values.size(); ++pi0 ) {
				for( std::size_t pi1 = 0; pi1 < parameter_values.size(); ++pi1 ) {
					parameters(0) = parameter_values[ pi0 ] ;
					parameters(1) = parameter_values[ pi1 ] ;

					ll.evaluate_at( parameters ) ;

					Matrix design_matrix = Matrix::Zero( nSamples, 2 ) ;
					design_matrix.col(0).setConstant( 1.0 ) ;
					design_matrix.col(1).segment( 0, nOnes ).setConstant( 1.0 ) ;

					Vector F1 = ( design_matrix * parameters ).array().exp() ;
					Vector ones = Vector::Constant( nSamples, 1.0 ) ;
					Vector const F0 = ones.array() / ( 1.0 + F1.array() ) ;
					F1 = F1.array() / ( 1.0 + F1.array() ) ;

					Vector Fpsi = F1 ;
					Fpsi.segment( nCases, nControls ) = F0.segment( nCases, nControls ) ;

					// all probabilities are one-half.
					double expected_ll = Fpsi.array().log().sum() ;
					check_close_or_small( expected_ll, ll.get_value_of_function(), 0.00000000001, 0.0001 ) ;

					Vector outcome_minus_f1 = ( outcome - F1 ) ;

					Vector const expected_first_derivative = ( outcome_minus_f1.transpose() * design_matrix ) ;
					Vector first_derivative = ll.get_value_of_first_derivative() ;
					check_close_or_small( expected_first_derivative(0), first_derivative(0), 0.00000000001, 0.0001 ) ;
					check_close_or_small( expected_first_derivative(1), first_derivative(1), 0.00000000001, 0.0001 ) ;

					Vector A = ( 1.0 - Fpsi.array() ) ;
					Vector B = ( 1.0 - Fpsi.array() ) * ( 1.0 - 2.0 * Fpsi.array() ) ;
					Matrix DD = Matrix::Zero( 2, 2 ) ;
					for( int sample = 0; sample < nSamples; ++sample ) {
						DD += ( B(sample) - A(sample) * A(sample) ) * design_matrix.row( sample ).transpose() * design_matrix.row( sample ) ;
					}
					Matrix const& second_derivative = ll.get_value_of_second_derivative() ;
					
					for( int i = 0; i < 2; ++i ) {
						for( int j = 0; j < 2; ++j ) {
							check_close_or_small( DD(i,j), second_derivative(i,j), 0.00000000001, 0.0001 ) ;
						}
					}
				}
			}
		}
	}
}

BOOST_AUTO_TEST_SUITE_END() ;
