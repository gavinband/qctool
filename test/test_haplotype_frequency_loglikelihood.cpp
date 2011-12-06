#include <Eigen/Core>
#include "HaplotypeFrequencyLogLikelihood.hpp"
#include "../../config.hpp"
#include "test_case.hpp"

using std::log ;

AUTO_TEST_CASE( test_haplotype_frequency_loglikelihood_value ) {
	typedef Eigen::MatrixXd Matrix ;
	typedef Eigen::VectorXd Vector ;
	
	double const tolerance = 0.0000000000001 ;
	
	Matrix table( 3, 3 ) ;

	{
		table <<
			1,	0,	0,
			0,	0,	0,
			0,	0,	0 ;

		HaplotypeFrequencyLogLikelihood ll( table ) ;

		Vector params( 3 ) ;
		params << 0.0, 0.0, 0.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), log( 1.0 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), Vector::Constant( 3, -2.0 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_second_derivative(), -2.0 * Vector::Constant( 3, -1.0 ) * Vector::Constant( 3, -1.0 ).transpose() ) ;
		params << 0.2, 0.0, 0.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), 2.0 * log( 0.8 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), Vector::Constant( 3, -2.0 / 0.8 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_second_derivative(), -( 2.0 / ( 0.8 * 0.8 )) * Vector::Constant( 3, -1.0 ) * Vector::Constant( 3, -1.0 ).transpose() ) ;
		params << 0.0, 0.2, 0.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), 2.0 * log( 0.8 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), Vector::Constant( 3, -2.0 / 0.8 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_second_derivative(), -( 2.0 / ( 0.8 * 0.8 )) * Vector::Constant( 3, -1.0 ) * Vector::Constant( 3, -1.0 ).transpose() ) ;
		params << 0.0, 0.0, 0.2 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), 2.0 * log( 0.8 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), Vector::Constant( 3, -2.0 / 0.8 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_second_derivative(), -( 2.0 / ( 0.8 * 0.8 )) * Vector::Constant( 3, -1.0 ) * Vector::Constant( 3, -1.0 ).transpose() ) ;
		params << 0.2, 0.2, 0.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), 2.0 * log( 0.6 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), Vector::Constant( 3, -2.0 / 0.6 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_second_derivative(), -( 2.0 / ( 0.6 * 0.6 )) * Vector::Constant( 3, -1.0 ) * Vector::Constant( 3, -1.0 ).transpose() ) ;
		params << 0.2, 0.0, 0.2 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), 2.0 * log( 0.6 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), Vector::Constant( 3, -2.0 / 0.6 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_second_derivative(), -( 2.0 / ( 0.6 * 0.6 )) * Vector::Constant( 3, -1.0 ) * Vector::Constant( 3, -1.0 ).transpose() ) ;
		params << 0.0, 0.2, 0.2 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), 2.0 * log( 0.6 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), Vector::Constant( 3, -2.0 / 0.6 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_second_derivative(), -( 2.0 / ( 0.6 * 0.6 )) * Vector::Constant( 3, -1.0 ) * Vector::Constant( 3, -1.0 ).transpose() ) ;
		params << 0.2, 0.2, 0.2 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_CLOSE( ll.get_value_of_function(), 2.0 * log( 0.4 ), tolerance ) ;
		BOOST_CHECK_SMALL( ( ll.get_value_of_first_derivative() - Vector::Constant( 3, -2.0 / 0.4 ) ).array().abs().maxCoeff(), tolerance ) ;
		BOOST_CHECK_SMALL( ( ll.get_value_of_second_derivative() - ( -( 2.0 / ( 0.4 * 0.4 )) * Vector::Constant( 3, -1.0 ) * Vector::Constant( 3, -1.0 ).transpose() ) ).array().abs().maxCoeff(), tolerance ) ;
		params << 1.0, 0.0, 0.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), -std::numeric_limits< double >::infinity() ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), Vector::Constant( 3, -std::numeric_limits< double >::infinity() ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_second_derivative(), -Matrix::Constant( 3, 3, std::numeric_limits< double >::infinity() ) ) ;
		params << 0.0, 1.0, 0.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), -std::numeric_limits< double >::infinity() ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), Vector::Constant( 3, -std::numeric_limits< double >::infinity() ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_second_derivative(), -Matrix::Constant( 3, 3, std::numeric_limits< double >::infinity() ) ) ;
		params << 0.0, 0.0, 1.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), -std::numeric_limits< double >::infinity() ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), Vector::Constant( 3, -std::numeric_limits< double >::infinity() ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_second_derivative(), -	Matrix::Constant( 3, 3, std::numeric_limits< double >::infinity() ) ) ;
	}

	{
		table <<
			0,	0,	0,
			0,	0,	0,
			0,	0,	1 ;

		HaplotypeFrequencyLogLikelihood ll( table ) ;

		Vector params( 3 ) ;
		params << 0.0, 0.0, 1.0 ;

		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), log( 1.0 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), 2.0 * Vector::Unit( 3, 2 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_second_derivative(), -2.0 * Vector::Unit( 3, 2 ) * Vector::Unit( 3, 2 ).transpose() ) ;
		params << 0.0, 0.0, 0.8 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), 2.0 * log( 0.8 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), 2.0 * Vector::Unit( 3, 2 ) / 0.8 ) ;
		params << 0.2, 0.0, 0.8 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), 2.0 * log( 0.8 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), 2.0 * Vector::Unit( 3, 2 ) / 0.8 ) ;
		params << 0.0, 0.2, 0.8 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), 2.0 * log( 0.8 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), 2.0 * Vector::Unit( 3, 2 ) / 0.8 ) ;
		params << 0.2, 0.0, 0.6 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), 2.0 * log( 0.6 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), 2.0 * Vector::Unit( 3, 2 ) / 0.6 ) ;
		params << 0.0, 0.2, 0.6 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), 2.0 * log( 0.6 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), 2.0 * Vector::Unit( 3, 2 ) / 0.6 ) ;
		params << 0.2, 0.2, 0.6 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), 2.0 * log( 0.6 ) ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), 2.0 * Vector::Unit( 3, 2 ) / 0.6 ) ;
		params << 0.2, 0.2, 0.4 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_CLOSE( ll.get_value_of_function(), 2.0 * log( 0.4 ), tolerance ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), 2.0 * Vector::Unit( 3, 2 ) / 0.4 ) ;
		params << 0.0, 0.0, 0.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), -std::numeric_limits< double >::infinity() ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_first_derivative(), Vector::Constant( 3, -std::numeric_limits< double >::infinity() )) ;
		params << 1.0, 0.0, 0.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), -std::numeric_limits< double >::infinity() ) ;
		params << 0.0, 1.0, 0.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), -std::numeric_limits< double >::infinity() ) ;
	}

	{
		table <<
			0,	0,	1,
			0,	0,	0,
			0,	0,	0 ;

		HaplotypeFrequencyLogLikelihood ll( table ) ;

		Vector params( 3 ) ;
		params << 1.0, 0.0, 0.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), log( 1.0 ) ) ;
		params << 0.8, 0.0, 0.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), 2.0 * log( 0.8 ) ) ;
		params << 0.8, 0.2, 0.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), 2.0 * log( 0.8 ) ) ;
		params << 0.8, 0.0, 0.2 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), 2.0 * log( 0.8 ) ) ;
		params << 0.6, 0.2, 0.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), 2.0 * log( 0.6 ) ) ;
		params << 0.6, 0.0, 0.2 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), 2.0 * log( 0.6 ) ) ;
		params << 0.6, 0.2, 0.2 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), 2.0 * log( 0.6 ) ) ;
		params << 0.4, 0.2, 0.4 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_CLOSE( ll.get_value_of_function(), 2.0 * log( 0.4 ), tolerance ) ;
		params << 0.0, 0.0, 0.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), -std::numeric_limits< double >::infinity() ) ;
		params << 0.0, 1.0, 0.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), -std::numeric_limits< double >::infinity() ) ;
		params << 0.0, 0.0, 1.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), -std::numeric_limits< double >::infinity() ) ;
	}

	{
		table <<
			0,	0,	0,
			0,	0,	0,
			1,	0,	0 ;

		HaplotypeFrequencyLogLikelihood ll( table ) ;

		Vector params( 3 ) ;
		params << 0.0, 1.0, 0.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), log( 1.0 ) ) ;
		params << 0.0, 0.8, 0.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), 2.0 * log( 0.8 ) ) ;
		params << 0.2, 0.8, 0.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), 2.0 * log( 0.8 ) ) ;
		params << 0.0, 0.8, 0.2 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), 2.0 * log( 0.8 ) ) ;
		params << 0.2, 0.6, 0.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), 2.0 * log( 0.6 ) ) ;
		params << 0.0, 0.6, 0.2 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), 2.0 * log( 0.6 ) ) ;
		params << 0.2, 0.6, 0.2 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), 2.0 * log( 0.6 ) ) ;
		params << 0.4, 0.4, 0.2 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_CLOSE( ll.get_value_of_function(), 2.0 * log( 0.4 ), tolerance ) ;
		params << 0.0, 0.0, 0.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), -std::numeric_limits< double >::infinity() ) ;
		params << 1.0, 0.0, 0.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), -std::numeric_limits< double >::infinity() ) ;
		params << 0.0, 0.0, 1.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), -std::numeric_limits< double >::infinity() ) ;
	}

	{
		table <<
			0,	0,	0,
			0,	0,	0,
			1,	0,	0 ;

		HaplotypeFrequencyLogLikelihood ll( table ) ;

		Vector params( 3 ) ;
		params << 0.0, 1.0, 0.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), log( 1.0 ) ) ;
		params << 0.0, 0.8, 0.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), 2.0 * log( 0.8 ) ) ;
		params << 0.2, 0.8, 0.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), 2.0 * log( 0.8 ) ) ;
		params << 0.0, 0.8, 0.2 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), 2.0 * log( 0.8 ) ) ;
		params << 0.2, 0.6, 0.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), 2.0 * log( 0.6 ) ) ;
		params << 0.0, 0.6, 0.2 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), 2.0 * log( 0.6 ) ) ;
		params << 0.2, 0.6, 0.2 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), 2.0 * log( 0.6 ) ) ;
		params << 0.4, 0.4, 0.2 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_CLOSE( ll.get_value_of_function(), 2.0 * log( 0.4 ), tolerance ) ;
		params << 0.0, 0.0, 0.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), -std::numeric_limits< double >::infinity() ) ;
		params << 1.0, 0.0, 0.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), -std::numeric_limits< double >::infinity() ) ;
		params << 0.0, 0.0, 1.0 ;
		ll.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( ll.get_value_of_function(), -std::numeric_limits< double >::infinity() ) ;
	}
}