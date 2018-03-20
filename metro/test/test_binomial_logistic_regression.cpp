#include <vector>
#include <iostream>
#include "test_case.hpp"
#include "metro/SampleRange.hpp"
#include "metro/regression/Design.hpp"
#include "metro/regression/BinomialLogistic.hpp"

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

typedef std::vector< std::string > Names ;

BOOST_AUTO_TEST_SUITE( test_binomial_logistic ) ;

namespace {
	void check_likelihood_etc_for_design1(
		int trials,
		int successes,
		double f0,
		double f1,
		double const ll,
		Eigen::MatrixXd const dll,
		Eigen::MatrixXd const ddll,
		int line
	) {
#if DEBUG_TESTS
		std::cerr << "Checking ll on line " << line << "...\n" ;
		std::cerr << "k = " << successes << ", n = " << trials << "\n"
			<< "f0 = " << f0 << ", f1 = " << f1 << ".\n" ;
#endif
		typedef Eigen::MatrixXd Matrix ;
		Matrix design_matrix = Matrix::Constant( 1, 2, 1.0 ) ;

		double const expected_l = std::pow( f1, successes ) * std::pow( f0, trials - successes ) ;
		double const expected_ll = std::log( expected_l ) ;
		BOOST_CHECK_CLOSE( expected_ll, ll, 0.00001 ) ;

#if DEBUG_TESTS
		std::cerr << "expected_l = " << expected_l << ", expected_ll = " << expected_ll << ", ll = " << ll << ".\n" ;
#endif
		Matrix const expected_dll
			= design_matrix.row(0)
			// * expected_l // this gets multiplied by and then divided again
			* ( successes * f0  - (trials - successes) * f1)
		;
		for( int j = 0; j < design_matrix.cols(); ++j ) {
			check_close_or_small( expected_dll(0,j), dll(j), 0.00001, 0.0001 ) ;
		}

		double const ddll_first_term
			= (
				std::pow( successes * f0  - (trials - successes) * f1, 2 )
				- trials * f0 * f1
			) ;
		double const ddll_second_term
			= std::pow(
				( successes * f0  - (trials - successes) * f1 ),
				2
			) ;
			
		Matrix const expected_ddll
			= (ddll_first_term - ddll_second_term)
			* design_matrix.row(0).transpose() * design_matrix.row(0)
		;

		for( int j = 0; j < design_matrix.cols(); ++j ) {
			for( int k = 0; k < design_matrix.cols(); ++k ) {
				check_close_or_small( expected_ddll(j,k), ddll(j,k), 0.00001, 0.0001 ) ;
			}
		}
		
#if DEBUG_TESTS
				std::cerr << "test_binomial_logistic_regression_one_sample:\n" ;
				
				std::cerr << "expected 2nd derivative = \n"
					<< expected_ddll << "\n" ;
				std::cerr << "actual 2nd derivative = \n"
					<< ddll << "\n" ;
#endif
	}
}

using metro::SampleRange ;
using namespace metro::regression ;
using metro::regression::Design ;
typedef Design::Matrix Matrix ;
typedef Design::Vector Vector ;
using std::exp ;
using std::log ;

// Test that when all samples are missing, the answer is zero.
AUTO_TEST_CASE( test_binomial_logisticregression_missing_outcome ) {
	typedef std::vector< metro::SampleRange > SampleRanges ;

	SampleRanges included_samples ;
	Matrix predictor_levels = Matrix::Zero( 2, 1 ) ;
	predictor_levels(0,0) = 1 ;
	predictor_levels(1,0) = 2 ;
	int const nOutcomes = 2 ;
	for( int nSamples = 0; nSamples < 100; ++nSamples ) {
		SampleRanges const all_samples( SampleRanges( 1, metro::SampleRange( 0, nSamples ))) ;
		Matrix covariates = Matrix::Zero( nSamples, 0 ) ;

		for( int nMissing = 0; nMissing < nSamples; ++nMissing ) {
			Matrix outcome = Matrix::Zero( nSamples, nOutcomes ) ;
			outcome.col(0).setConstant(1) ;
			outcome.col(1).setZero() ;
			SampleRanges nonmissing_samples ;
			nonmissing_samples.push_back( metro::SampleRange( nMissing, nSamples )) ;
			Names outcome_names ;
			for( int j = 0; j < nOutcomes; ++j ) {
				outcome_names.push_back( "outcome" + std::to_string( j )) ;
			}
			BinomialLogistic ll(
				Design::create(
					outcome, nonmissing_samples, outcome_names,
					covariates, all_samples, Names(),
					Names( 1, "predictor" )
				)
			) ;
		
#if DEBUG_TESTS
			std::cerr << "nSamples = " << nSamples << ", nMissing = " << nMissing << ", nOutcomes = " << nOutcomes << ".\n" ;
			std::cerr << "ll =\n" << ll.get_summary() << "\n" ;
#endif
		
			Matrix predictor_probabilities = Matrix::Zero( nSamples, 2 ) ;
			predictor_probabilities.col(1).setConstant( 1.0 ) ;
			ll.design().set_predictors( predictor_levels, predictor_probabilities, included_samples ) ;
			// Now test values.
			Vector parameters = Vector::Zero( 2 ) ;
			parameters(0) = 1.0 ;
			parameters(1) = -1.0 ;
			ll.evaluate_at( parameters ) ;
		
			double const expected_ll = 0.0 ;
			Vector const expected_dll = Vector::Zero(2) ;
			Matrix const expected_ddll = Vector::Zero(2,2) ;
			BOOST_CHECK_EQUAL( ll.get_value_of_function(), expected_ll ) ;
			BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), expected_dll ) ;
			BOOST_CHECK_EQUAL( ll.get_value_of_second_derivative(), expected_ddll ) ;
		}
	}
}

AUTO_TEST_CASE( test_binomial_logisticregression_zero_samples ) {
	std::vector< metro::SampleRange > included_samples ;
	Matrix predictor_levels = Matrix::Zero( 2, 1 ) ;
	predictor_levels(0,0) = 1 ;
	predictor_levels(1,0) = 2 ;
	int const nOutcomes = 2 ;
	for( int nSamples = 0; nSamples < 10; ++nSamples ) {
		Matrix outcome = Matrix::Zero( nSamples, nOutcomes ) ;
		outcome.col(0).setConstant(1) ;
		outcome.col(1).setZero() ;

		Matrix covariates = Matrix::Zero( nSamples, 0 ) ;
		Names outcome_names ;
		for( int j = 0; j < nOutcomes; ++j ) {
			outcome_names.push_back( "outcome" + std::to_string( j )) ;
		}

		BinomialLogistic ll(
			Design::create(
				outcome, included_samples, outcome_names,
				covariates, included_samples, Names(),
				Names( 1, "predictor" )
			)
		) ;
	
		Matrix predictor_probabilities = Matrix::Zero( nSamples, 2 ) ;
		ll.design().set_predictors( predictor_levels, predictor_probabilities, included_samples ) ;
		// Now test values.
		Vector parameters = Vector::Zero( 2 ) ;
		ll.evaluate_at( parameters ) ;
	
		double const expected_ll = 0.0 ;
		Vector const expected_dll = Vector::Zero(2) ;
		Matrix const expected_ddll = Vector::Zero(2,2) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), expected_ll ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), expected_dll ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_second_derivative(), expected_ddll ) ;
	}
}

AUTO_TEST_CASE( test_binomial_logistic_regression_one_sample ) {
	using metro::SampleRange ;
	using namespace metro::regression ;
	using metro::regression::Design ;
	typedef Design::Matrix Matrix ;
	typedef Design::Vector Vector ;
	using std::exp ;
	using std::log ;
	
	int const nSamples = 1 ;
	std::vector< metro::SampleRange > all_samples( 1, metro::SampleRange( 0, 1 ) ) ;
	Matrix predictor_levels = Matrix::Zero( 2, 1 ) ;
	predictor_levels(0,0) = 1 ;
	predictor_levels(1,0) = 2 ;
	int const nOutcomes = 2 ;
	for( int trials = 0; trials < 100; ++trials ) {
		for( int successes = 0; successes <= trials; ++successes ) {
			Matrix outcome = Matrix::Zero( nSamples, nOutcomes ) ;
			outcome(0,0) = trials - successes ;
			outcome(0,1) = successes ;

			Matrix outcome_nonmissingness = Matrix::Constant( nSamples, nOutcomes, 1.0 ) ;
			Matrix covariates = Matrix::Zero( nSamples, 0 ) ;
			Matrix covariate_nonmissingness = Matrix::Constant( nSamples, 0, 1.0 ) ;
			Names outcome_names ;
			for( int j = 0; j < nOutcomes; ++j ) {
				outcome_names.push_back( "outcome" + std::to_string( j )) ;
			}

			BinomialLogistic ll(
				Design::create(
					outcome, all_samples, outcome_names,
					covariates, all_samples, Names(),
					Names( 1, "predictor" )
				)
			) ;
			
			// Set up a single predictor, equal to 1 with certainty.
			Matrix predictor_probabilities = Matrix::Zero( 1, 2 ) ;
			predictor_probabilities(0, 0) = 1.0 ;
			ll.design().set_predictors( predictor_levels, predictor_probabilities, all_samples ) ;

			int const nParameters = 2 * ( nOutcomes - 1 ) ;
			// Now test values.
			{
				Vector parameters = Vector::Zero( nParameters ) ;
				ll.evaluate_at( parameters ) ;
				// All probs are equal to 1/(M+1).	So outcome probability equals this.
				Matrix parameter_matrix = parameters ;
				parameter_matrix.resize( 2, nOutcomes - 1 ) ;
#if DEBUG_TESTS
				std::cerr << "parameters = " << parameters.transpose() << ".\n" ;
				std::cerr << "outcome = " << outcome.transpose() << ".\n" ;
#endif
	
				// params are zero so probs are 1/2
				double const f0 = 0.5 ;
				double const f1 = 0.5 ;

#if DEBUG_TESTS
				std::cerr << "##############\n" ;
				std::cerr << "(successes, trials) = " << successes << ", " << trials << "\n" ;
				std::cerr << "parameters = " << parameters.transpose() << ".\n" ;
				std::cerr << "outcome = " << outcome.transpose() << ".\n" ;
				std::cerr << "design = " << ll.design().matrix() << ".\n" ;
#endif
				check_likelihood_etc_for_design1(
					trials, successes,
					f0, f1,
					ll.get_value_of_function(),
					ll.get_value_of_first_derivative(),
					ll.get_value_of_second_derivative(),
					__LINE__
				) ;
			}

			for( int nonzero_i = 0; nonzero_i < nParameters; ++nonzero_i )
			{

				Vector parameters = Vector::Zero( 2 * ( nOutcomes - 1 ) ) ;
				parameters( nonzero_i ) = log(2.0) ;
				ll.evaluate_at( parameters ) ;

				// f_i will equal 1/(M+2), or 2/(M+2) for i=floor(nonzero_i/2).
				Matrix parameter_matrix = parameters ;
				parameter_matrix.resize( 2, nParameters ) ;

				// params are zero so probs are 1/2
				Matrix design_matrix = Matrix::Constant( 1, 2, 1.0 ) ;
				double const f0 = 1.0 / (1.0 + std::exp( design_matrix.row(0) * parameters )) ;
				double const f1 = 1 - f0 ;

#if DEBUG_TESTS
				std::cerr << "##############\n" ;
				std::cerr << "(successes, trials) = " << successes << ", " << trials << "\n" ;
				std::cerr << "nonzero = " << nonzero_i << ".\n" ;
				std::cerr << "parameters = " << parameters.transpose() << ".\n" ;
				std::cerr << "outcome = " << outcome.transpose() << ".\n" ;
				std::cerr << "design = " << ll.design().matrix() << ".\n" ;
#endif
				check_likelihood_etc_for_design1(
					trials, successes,
					f0, f1,
					ll.get_value_of_function(),
					ll.get_value_of_first_derivative(),
					ll.get_value_of_second_derivative(),
					__LINE__
				) ;
			}
		}
	}
}


AUTO_TEST_CASE( test_binomial_logisticregression_two_samples ) {
	using metro::SampleRange ;
	using namespace metro::regression ;
	using metro::regression::Design ;
	typedef Design::Matrix Matrix ;
	typedef Design::Vector Vector ;
	using std::exp ;
	using std::log ;
	
	typedef std::vector< metro::SampleRange > SampleRanges ;
	SampleRanges const all_samples( 1, metro::SampleRange( 0, 2 )) ;
	
	{
		int const nSamples = 2 ;
		Matrix outcome = Matrix::Zero(nSamples, 2) ;
		outcome <<
			1, 0,
			0, 1 ;
		Matrix covariates = Matrix::Zero( nSamples, 0 ) ;

		BinomialLogistic ll(
			Design::create(
				outcome, all_samples, Names({"outcome=0", "outcome=1"}),
				covariates, all_samples, Names(),
				Names({"predictor"})
			)
		) ;
		
		Matrix predictor_levels = Matrix::Zero( 2, 1 ) ;
		predictor_levels(0,0) = 0 ;
		predictor_levels(1,0) = 1 ;

		// Set half the sample to certainly 0, half to certainly 1.
		Matrix predictor_probabilities = Matrix::Zero( nSamples, 2 ) ;
		predictor_probabilities( 0, 0 ) = 1.0 ;
		predictor_probabilities( 1, 1 ) = 1.0 ;

		ll.design().set_predictors( predictor_levels, predictor_probabilities, all_samples ) ;
		
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
			
			// Compute first derivative
			Vector B = Vector::Constant( 2, 1.0 ) - F ;
			B(0) = -B(0) ;
			Vector const expected_first_derivative = ( B.transpose() * design_matrix ) ;
			Vector const first_derivative = ll.get_value_of_first_derivative() ;
			BOOST_CHECK_CLOSE( expected_first_derivative(0), first_derivative(0), 0.0000001 ) ;
			BOOST_CHECK_CLOSE( expected_first_derivative(1), first_derivative(1), 0.0000001 ) ;
			
			// Compute second derivative
			Vector C = ( Vector::Constant( 2, 1.0 ) - 2.0 * F ).array() * B.array() ;
			Matrix DD = Matrix::Zero( 2, 2 ) ;
			for( int i = 0; i < design_matrix.rows(); ++i ) {
				DD += ( C(i) - B(i) * B(i) ) * design_matrix.row(i).transpose() * design_matrix.row(i) ;
			}
			
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
			Matrix DD = Matrix::Zero( 2, 2 ) ;
			for( int i = 0; i < design_matrix.rows(); ++i ) {
				DD += ( C(i) - B(i) * B(i) ) * design_matrix.row(i).transpose() * design_matrix.row(i) ;
			}
			Matrix const& second_derivative = ll.get_value_of_second_derivative() ;
			BOOST_CHECK_CLOSE( DD(0,0), second_derivative(0,0), 0.0000001 ) ;
			BOOST_CHECK_CLOSE( DD(0,1), second_derivative(0,1), 0.0000001 ) ;
			BOOST_CHECK_CLOSE( DD(1,0), second_derivative(1,0), 0.0000001 ) ;
			BOOST_CHECK_CLOSE( DD(1,1), second_derivative(1,1), 0.0000001 ) ;
		}
	}
}


AUTO_TEST_CASE( test_binomial_logisticregression_certain_predictors ) {
	using metro::SampleRange ;
	using namespace metro::regression ;
	using metro::regression::Design ;
	typedef Design::Matrix Matrix ;
	typedef Design::Vector Vector ;
	typedef std::vector< metro::SampleRange > SampleRanges ;
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
	
	for( int nCases = 0; nCases < 50; ++nCases ) {
		for( int nControls = 0; nControls < 50; ++nControls ) {
#if DEBUG_TESTS
			std::cerr << "nCases = " << nCases << ", nControls = " << nControls << ".\n" ;
#endif
			int nSamples = nCases + nControls ;

			SampleRanges const all_samples( 1, metro::SampleRange( 0, nSamples )) ;

			// contruct outcome
			Matrix outcome = Matrix::Zero(nSamples,2) ;
			{
				int i = 0 ;
				for( ; i < nCases; ++i ) {
					outcome(i,1) = 1 ;
				}
				for( ; i < nSamples; ++i ) {
					outcome(i,0) = 1 ;
				}
			}

			Matrix covariates = Matrix::Zero( nSamples, 0 ) ;

			BinomialLogistic ll(
				Design::create(
					outcome, all_samples, Names({"outcome=0", "outcome=1"}),
					covariates, all_samples, std::vector< std::string >(),
					std::vector< std::string >( 1, "predictor" )
				)
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
			
			ll.design().set_predictors( predictor_levels, predictor_probabilities, all_samples ) ;
		
			// Ok we have our likelihood.  Now compute with it.
			{
				ll.evaluate_at( Vector::Zero( 2 ) ) ;
				// all probabilities are one-half.
				double expected_ll = nSamples * std::log( 0.5 ) ;
				BOOST_CHECK_CLOSE( ll.get_value_of_function(), expected_ll, 0.0000001 ) ;

				Vector one_minus_f1 = ( outcome.col(1).array() - 0.5 ) ;
				Matrix design_matrix = Matrix::Zero( nSamples, 2 ) ;
				design_matrix.col(0).setConstant( 1.0 ) ;
				design_matrix.col(1).segment( 0, nOnes ).setConstant( 1.0 ) ;
				
				Vector const expected_first_derivative = ( one_minus_f1.transpose() * design_matrix ) ;
				Vector const first_derivative = ll.get_value_of_first_derivative() ;
				BOOST_CHECK_CLOSE( expected_first_derivative(0), first_derivative(0), 0.0000001 ) ;
				BOOST_CHECK_CLOSE( expected_first_derivative(1), first_derivative(1), 0.0000001 ) ;
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

					Vector outcome_minus_f1 = ( outcome.col(1) - F1 ) ;

					Vector const expected_first_derivative = ( outcome_minus_f1.transpose() * design_matrix ) ;
					Vector first_derivative = ll.get_value_of_first_derivative() ;
					check_close_or_small( expected_first_derivative(0), first_derivative(0), 0.00000000001, 0.0001 ) ;
					check_close_or_small( expected_first_derivative(1), first_derivative(1), 0.00000000001, 0.0001 ) ;
				}
			}
		}
	}
}

BOOST_AUTO_TEST_SUITE_END() ;
