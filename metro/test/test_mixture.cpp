
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <Eigen/Core>
#include "test_case.hpp"
#include "metro/likelihood/Mixture.hpp"
#include "metro/IndependentObservationLogLikelihood.hpp"

typedef Eigen::MatrixXd Matrix ;
typedef Eigen::VectorXd Vector ;

double const infinity = std::numeric_limits< double >::infinity() ;

// #define DEBUG_MULTIVARIATE_T 1

BOOST_AUTO_TEST_SUITE( test_mixture )

namespace {
	struct ConstantLogLikelihood: public metro::IndependentObservationLogLikelihood< double, Eigen::VectorXd, Eigen::MatrixXd > {
		typedef std::auto_ptr< ConstantLogLikelihood > UniquePtr ;

		ConstantLogLikelihood( Eigen::MatrixXd const& data ):
			m_n( data.rows() ),
			m_parameters( Eigen::VectorXd::Constant( 1, 0 ) )
		{}

		void set_data( Eigen::MatrixXd const& data ) {
			m_n = data.rows() ;
		}

		void evaluate_at( Vector const& parameters ) {
			evaluate_at( parameters, metro::DataRange( 0, m_n )) ;
		}

		void evaluate_at( Vector const& parameters, metro::DataSubset const& subset ) {
			assert( parameters.size() == 1 ) ;
			m_parameters = parameters ;
		}
		
		double get_value_of_function() const {
			return m_n * m_parameters(0) ;
		}

		Eigen::VectorXd get_value_of_first_derivative() const {
			assert(0) ;
		}
		Eigen::MatrixXd get_value_of_second_derivative() const {
			assert(0) ;
		}
		
		Eigen::VectorXd parameters() const {
			return m_parameters ;
		}

		std::string get_spec() const {
			return "ConstantLogLikelihood" ;
		}

		virtual void get_terms_of_function( MatrixRef result ) const {
			assert( result.rows() == m_n ) ;
			result.setConstant( m_parameters(0) ) ;
		} ;
		
	private:
		int m_n ;
		Eigen::VectorXd m_parameters ;
	} ;
}

AUTO_TEST_CASE( test_parameters1 ) {
	// A single data point with probability 1 under one component and 100 under the other
	Eigen::MatrixXd data = Eigen::MatrixXd::Zero( 1, 0 ) ;

	metro::likelihood::Mixture< double, Eigen::VectorXd, Eigen::MatrixXd > mixture( data ) ;
	// ConstantLogLikelihood has a single parameter, by default equal to 0.
	BOOST_CHECK_EQUAL( mixture.parameters().size(), 0 ) ;
	mixture.add_component( "one", 1, ConstantLogLikelihood::UniquePtr( new ConstantLogLikelihood( data ) ) ) ;
	BOOST_CHECK_EQUAL( mixture.parameters().size(), 1 ) ;
	BOOST_CHECK_EQUAL( mixture.parameters()(0), 0 ) ;

	mixture.add_component( "two", 1, ConstantLogLikelihood::UniquePtr( new ConstantLogLikelihood( data ) ) ) ;
	BOOST_CHECK_EQUAL( mixture.parameters().size(), 3 ) ;
	BOOST_CHECK_EQUAL( mixture.parameters()(0), 0 ) ;
	BOOST_CHECK_EQUAL( mixture.parameters()(1), 0.5 ) ;
	BOOST_CHECK_EQUAL( mixture.parameters()(2), 0 ) ;

	mixture.add_component( "three", 1, ConstantLogLikelihood::UniquePtr( new ConstantLogLikelihood( data ) ) ) ;
	BOOST_CHECK_EQUAL( mixture.parameters().size(), 5 ) ;
	BOOST_CHECK_EQUAL( mixture.parameters()(0), 0 ) ;
	BOOST_CHECK_CLOSE( mixture.parameters()(1), 1/3.0, 1E-12 ) ;
	BOOST_CHECK_EQUAL( mixture.parameters()(2), 0 ) ;
	BOOST_CHECK_CLOSE( mixture.parameters()(3), 1/3.0, 1E-12 ) ;
	BOOST_CHECK_EQUAL( mixture.parameters()(4), 0 ) ;

	mixture.add_component( "three", 1, ConstantLogLikelihood::UniquePtr( new ConstantLogLikelihood( data ) ) ) ;
	BOOST_CHECK_EQUAL( mixture.parameters().size(), 7 ) ;
	BOOST_CHECK_EQUAL( mixture.parameters()(0), 0 ) ;
	BOOST_CHECK_CLOSE( mixture.parameters()(1), 1/4.0, 1E-12 ) ;
	BOOST_CHECK_EQUAL( mixture.parameters()(2), 0 ) ;
	BOOST_CHECK_CLOSE( mixture.parameters()(3), 1/4.0, 1E-12 ) ;
	BOOST_CHECK_EQUAL( mixture.parameters()(4), 0 ) ;
	BOOST_CHECK_CLOSE( mixture.parameters()(5), 1/4.0, 1E-12 ) ;
	BOOST_CHECK_EQUAL( mixture.parameters()(6), 0 ) ;
}

AUTO_TEST_CASE( test_parameters2 ) {
	// A single data point with probability 1 under one component and 100 under the other
	Eigen::MatrixXd data = Eigen::MatrixXd::Zero( 1, 0 ) ;

	metro::likelihood::Mixture< double, Eigen::VectorXd, Eigen::MatrixXd > mixture( data ) ;
	mixture.add_component( "one", 1, ConstantLogLikelihood::UniquePtr( new ConstantLogLikelihood( data ) ) ) ;
	{
		Eigen::VectorXd params( 1 ) ;
		params << 5 ;
		mixture.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( mixture.parameters()(0), 5 ) ;
	}
	mixture.add_component( "two", 1, ConstantLogLikelihood::UniquePtr( new ConstantLogLikelihood( data ) ) ) ;
	{
		Eigen::VectorXd params( 3 ) ;
		params << 5, 0.5, 3 ;
		mixture.evaluate_at( params ) ;
		BOOST_CHECK_EQUAL( mixture.parameters()(0), 5 ) ;
		BOOST_CHECK_EQUAL( mixture.parameters()(1), 0.5 ) ;
		BOOST_CHECK_EQUAL( mixture.parameters()(2), 3 ) ;
	}
}

AUTO_TEST_CASE( test_one_component_mixture ) {
	// A single data point with probability 1 under one component and 100 under the other
	Eigen::MatrixXd data = Eigen::MatrixXd::Zero( 1, 0 ) ;

	metro::likelihood::Mixture< double, Eigen::VectorXd, Eigen::MatrixXd > mixture( data ) ;
	mixture.add_component( "one", 1, ConstantLogLikelihood::UniquePtr( new ConstantLogLikelihood( data ) ) ) ;
	ConstantLogLikelihood ll( data ) ;
	ll.set_data( data ) ;
	mixture.set_data( data ) ;

	Eigen::VectorXd parameters = Eigen::VectorXd::Constant( 1, 0.5 ) ;
	mixture.evaluate_at( parameters ) ;
	ll.evaluate_at( parameters ) ;
	// should be 0.5 * 1 + 0.5 * 100 = 50.5
	BOOST_CHECK_CLOSE( mixture.get_value_of_function(), ll.get_value_of_function(), 1E-6 ) ;
	BOOST_CHECK_CLOSE( mixture.get_value_of_function(), 0.5, 1E-6 ) ;
}

AUTO_TEST_CASE( test_two_component_mixture ) {
	Eigen::MatrixXd data = Eigen::MatrixXd::Zero( 1, 0 ) ;

	// Equally-weighted mixture...
	metro::likelihood::Mixture< double, Eigen::VectorXd, Eigen::MatrixXd > mixture( data ) ;
	mixture.add_component( "one", 1, new ConstantLogLikelihood( data ) ) ;
	mixture.add_component( "two", 1, new ConstantLogLikelihood( data ) ) ;
	mixture.set_data( data ) ;
	
	Eigen::VectorXd parameters( 3 ) ;
	// 50:50 mixture of 1 and 100.
	parameters << std::log( 1 ), 0.5, std::log( 100 ) ;
	mixture.evaluate_at( parameters ) ;
	BOOST_CHECK_CLOSE( mixture.get_value_of_function(), std::log( 50.5 ), 1E-7 ) ;

	// 50:50 Mixture of 8 and 16.
	parameters << std::log( 8 ), 0.5, std::log( 16 ) ;
	mixture.evaluate_at( parameters ) ;
	BOOST_CHECK_CLOSE( mixture.get_value_of_function(), std::log( 12 ), 1E-7 ) ;

	// 0:100 Mixture of 8 and 16.
	parameters << std::log( 8 ), 0, std::log( 16 ) ;
	mixture.evaluate_at( parameters ) ;
	BOOST_CHECK_CLOSE( mixture.get_value_of_function(), std::log( 16 ), 1E-7 ) ;

	// 1:0 Mixture of 8 and 16.
	parameters << std::log( 8 ), 1, std::log( 16 ) ;
	mixture.evaluate_at( parameters ) ;
	BOOST_CHECK_CLOSE( mixture.get_value_of_function(), std::log( 8 ), 1E-7 ) ;
}

BOOST_AUTO_TEST_SUITE_END()
