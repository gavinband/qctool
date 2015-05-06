
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

		ConstantLogLikelihood():
			m_parameters( 1 )
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
		
		Eigen::VectorXd get_parameters() const {
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

AUTO_TEST_CASE( test_one_component_mixture ) {
	// A single data point with probability 1 under one component and 100 under the other
	Eigen::MatrixXd data = Eigen::MatrixXd::Zero( 1, 0 ) ;

	metro::likelihood::Mixture< double, Eigen::VectorXd, Eigen::MatrixXd > mixture( data ) ;
	mixture.add_component( "one", 1, ConstantLogLikelihood::UniquePtr( new ConstantLogLikelihood ) ) ;
	ConstantLogLikelihood ll ;
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
	mixture.add_component( "one", 1, new ConstantLogLikelihood ) ;
	mixture.add_component( "two", 1, new ConstantLogLikelihood ) ;
	mixture.set_data( data ) ;
	
	Eigen::VectorXd parameters( 2 ) ;
	// Mixture of 1 and 100.
	parameters << std::log( 1 ), std::log( 100 ) ;
	mixture.evaluate_at( parameters ) ;
	BOOST_CHECK_CLOSE( mixture.get_value_of_function(), std::log( 50.5 ), 1E-7 ) ;

	// Mixture of -100 and 20.
	parameters << std::log( 8 ), std::log( 16 ) ;
	mixture.evaluate_at( parameters ) ;
	BOOST_CHECK_CLOSE( mixture.get_value_of_function(), std::log( 12 ), 1E-7 ) ;
}

BOOST_AUTO_TEST_SUITE_END()
