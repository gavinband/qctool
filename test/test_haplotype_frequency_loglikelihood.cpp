#include <Eigen/Core>
#include "genfile/Error.hpp"
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

namespace {
	void test_haplotype_frequency_estimation( Eigen::MatrixXd const& genotypes, Eigen::VectorXd const& expected ) {
		HaplotypeFrequencyLogLikelihood ll( genotypes ) ;
		std::cerr << "test_haplotype_frequency_estimation: genotypes: " << genotypes
			<< "\nexpected: " << expected << "\n" ;
		std::cerr << "     got: " << ll.get_MLE_by_EM() << ".\n" ;
		BOOST_CHECK_SMALL( ( ll.get_MLE_by_EM() - expected ).maxCoeff(), 0.000000000001 ) ;
	}
}

AUTO_TEST_CASE( test_haplotype_frequency_exceptions ) {
	{
		Eigen::MatrixXd	genotypes = Eigen::MatrixXd::Zero( 3, 3 ) ;
		BOOST_CHECK_THROW( { HaplotypeFrequencyLogLikelihood ll( genotypes ) ; }, genfile::BadArgumentError ) ;
	}
}

AUTO_TEST_CASE( test_single_person_haplotype_frequency_estimation ) {
	Eigen::MatrixXd	genotypes( 3, 3 ) ;
	Eigen::VectorXd params( 3 ) ;

	genotypes.setZero() ; params.setZero() ;
	genotypes( 0, 0 ) = 1 ;
	test_haplotype_frequency_estimation( genotypes, params ) ;

	genotypes.setZero() ; params.setZero() ;
	genotypes( 0, 1 ) = 1 ;
	params(0) = 0.5 ;
	test_haplotype_frequency_estimation( genotypes, params ) ;

	genotypes.setZero() ; params.setZero() ;
	genotypes( 0, 2 ) = 1 ;
	params(0) = 1 ;
	test_haplotype_frequency_estimation( genotypes, params ) ;

	genotypes.setZero() ; params.setZero() ;
	genotypes( 1, 0 ) = 1 ;
	params(1) = 0.5 ;
	test_haplotype_frequency_estimation( genotypes, params ) ;

	genotypes.setZero() ; params.setZero() ;
	genotypes( 1, 1 ) = 1 ;
	// both SNPs heterozygores; half a chance of AB_ab and half of Ab_aB
	params(0) = 0.25 ;
	params(1) = 0.25 ;
	params(2) = 0.25 ;
	test_haplotype_frequency_estimation( genotypes, params ) ;

	genotypes.setZero() ; params.setZero() ;
	genotypes( 1, 2 ) = 1 ;
	params(0) = 0.5 ;
	params(2) = 0.5 ;
	test_haplotype_frequency_estimation( genotypes, params ) ;

	genotypes.setZero() ; params.setZero() ;
	genotypes( 2, 0 ) = 1 ;
	params(1) = 1 ;
	test_haplotype_frequency_estimation( genotypes, params ) ;

	genotypes.setZero() ; params.setZero() ;
	genotypes( 2, 1 ) = 1 ;
	params(1) = 0.5 ;
	params(2) = 0.5 ;
	test_haplotype_frequency_estimation( genotypes, params ) ;

	genotypes.setZero() ; params.setZero() ;
	genotypes( 2, 2 ) = 1 ;
	params(2) = 1 ;
	test_haplotype_frequency_estimation( genotypes, params ) ;
}

AUTO_TEST_CASE( test_N_person_haplotype_frequency_estimation ) {
	Eigen::MatrixXd	genotypes( 3, 3 ) ;
	Eigen::VectorXd params( 3 ) ;

	for( std::size_t N = 1; N < 10; ++N ) {
		genotypes.setZero() ; params.setZero() ;
		genotypes( 0, 0 ) = 1 ;
		test_haplotype_frequency_estimation( genotypes, params ) ;

		genotypes.setZero() ; params.setZero() ;
		genotypes( 0, 1 )= N ;
		params(0) = 0.5 ;
		test_haplotype_frequency_estimation( genotypes, params ) ;

		genotypes.setZero() ; params.setZero() ;
		genotypes( 0, 2 )= N ;
		params(0) = 1 ;
		test_haplotype_frequency_estimation( genotypes, params ) ;

		genotypes.setZero() ; params.setZero() ;
		genotypes( 1, 0 )= N ;
		params(1) = 0.5 ;
		test_haplotype_frequency_estimation( genotypes, params ) ;

		genotypes.setZero() ; params.setZero() ;
		genotypes( 1, 1 )= N ;
		// both SNPs heterozygores; half a chance of AB_ab and half of Ab_aB
		params(0) = 0.25 ;
		params(1) = 0.25 ;
		params(2) = 0.25 ;
		test_haplotype_frequency_estimation( genotypes, params ) ;

		genotypes.setZero() ; params.setZero() ;
		genotypes( 1, 2 )= N ;
		params(0) = 0.5 ;
		params(2) = 0.5 ;
		test_haplotype_frequency_estimation( genotypes, params ) ;

		genotypes.setZero() ; params.setZero() ;
		genotypes( 2, 0 )= N ;
		params(1) = 1 ;
		test_haplotype_frequency_estimation( genotypes, params ) ;

		genotypes.setZero() ; params.setZero() ;
		genotypes( 2, 1 )= N ;
		params(1) = 0.5 ;
		params(2) = 0.5 ;
		test_haplotype_frequency_estimation( genotypes, params ) ;

		genotypes.setZero() ; params.setZero() ;
		genotypes( 2, 2 )= N ;
		params(2) = 1 ;
		test_haplotype_frequency_estimation( genotypes, params ) ;
	}
}

AUTO_TEST_CASE( test_2_person_haplotype_frequency_estimation ) {
	Eigen::MatrixXd	genotypes( 3, 3 ) ;
	Eigen::VectorXd params( 3 ) ;

	genotypes.setZero() ; params.setZero() ;
	genotypes( 0, 0 ) = 1 ;
	genotypes( 2, 1 ) = 1 ;
	params(1) = 0.25 ;
	params(2) = 0.25 ;
	test_haplotype_frequency_estimation( genotypes, params ) ;

	genotypes.setZero() ; params.setZero() ;
	genotypes( 2, 2 ) = 1 ;
	genotypes( 2, 1 ) = 1 ;
	params(1) = 0.25 ;
	params(2) = 0.75 ;
	test_haplotype_frequency_estimation( genotypes, params ) ;

	genotypes.setZero() ; params.setZero() ;
	genotypes( 2, 2 ) = 1 ;
	genotypes( 1, 0 ) = 1 ;
	params(1) = 0.25 ;
	params(2) = 0.5 ;
	test_haplotype_frequency_estimation( genotypes, params ) ;

	genotypes.setZero() ; params.setZero() ;
	genotypes( 2, 2 ) = 1 ;
	genotypes( 1, 1 ) = 1 ;
	params(0) = 0.125 ;
	params(1) = 0.125 ;
	params(2) = 0.5 + ( 0.25 * 0.5 );
	test_haplotype_frequency_estimation( genotypes, params ) ;
}
