#include <vector>
#include <iostream>
#include "test_case.hpp"
#include "metro/SampleRange.hpp"
#include "metro/regression/Design.hpp"
#include "metro/regression/Logistic.hpp"

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

BOOST_AUTO_TEST_SUITE( test_logistic ) ;

AUTO_TEST_CASE( test_logisticregression_two_samples ) {
	using metro::SampleRange ;
	using namespace metro::regression ;
	using metro::regression::Design ;
	typedef Design::Matrix Matrix ;
	typedef Design::Vector Vector ;
	using std::exp ;
	using std::log ;
	
	{
		int const nSamples = 2 ;
		Matrix outcome = Matrix::Zero(nSamples, 2) ;
		outcome <<
			1, 0,
			0, 1 ;
		Matrix outcome_nonmissingness = Vector::Constant( nSamples, 1.0 ) ;
		Matrix covariates = Matrix::Zero( nSamples, 0 ) ;
		Matrix covariate_nonmissingness = Matrix::Constant( nSamples, 0, 1.0 ) ;

		Logistic ll(
			Design::create(
				outcome, outcome_nonmissingness, Names({"outcome=0", "outcome=1"}),
				covariates, covariate_nonmissingness, Names(),
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


AUTO_TEST_CASE( test_logisticregression_certain_predictors ) {
	using metro::SampleRange ;
	using namespace metro::regression ;
	using metro::regression::Design ;
	typedef Design::Matrix Matrix ;
	typedef Design::Vector Vector ;
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
			Matrix outcome_nonmissingness = Matrix::Constant( nSamples, 2, 1.0 ) ;

			Matrix covariates = Matrix::Zero( nSamples, 0 ) ;
			Matrix covariate_nonmissingness = Matrix::Constant( nSamples, 0, 1.0 ) ;

			Logistic ll(
				Design::create(
					outcome, outcome_nonmissingness, Names({"outcome=0", "outcome=1"}),
					covariates, covariate_nonmissingness, std::vector< std::string >(),
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
			
			std::vector< metro::SampleRange > const included_samples = std::vector< metro::SampleRange >( 1, metro::SampleRange( 0, nSamples ) ) ;
		
			ll.set_predictor_levels( predictor_levels, predictor_probabilities, included_samples ) ;
		
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
