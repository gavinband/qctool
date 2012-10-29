
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <map>
#include <cmath>
#include <cassert>
#include <boost/function.hpp>
#include <Eigen/Core>
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "integration/NewtonRaphson.hpp"
#include "components/CallComparerComponent/GenotypeFrequencyTestCallComparer.hpp"
#include "metro/likelihood/Multinomial.hpp"
#include "metro/likelihood/ProductOfMultinomials.hpp"
#include "metro/rBF.hpp"

#define Genotype_FREQUENCY_TEST_CALL_COMPARER_DEBUG 0

namespace {
	typedef Eigen::VectorXd Vector ;
	typedef Eigen::MatrixXd Matrix ;
	typedef metro::likelihood::Multinomial< double, Vector, Matrix > Multinomial ;
	typedef metro::likelihood::ProductOfMultinomials< double, Vector, Matrix > ProductOfIndependentMultinomials ;
}

GenotypeFrequencyTestCallComparer::GenotypeFrequencyTestCallComparer():
	m_threshhold( 0.9 ),
 	m_chi_squared( 2.0 )
{}

void GenotypeFrequencyTestCallComparer::compare(
	Eigen::MatrixXd const& left,
	Eigen::MatrixXd const& right,
	Callback callback
) const {
	// Make table of counts
	assert( left.rows() == right.rows() ) ;
	assert( left.cols() == 3 ) ;
	assert( right.cols() == 3 ) ;

	std::size_t const N = left.rows() ;
	
	Matrix table = Matrix::Zero( 2, 3 ) ;

	for( std::size_t i = 0; i < N; ++i ) {
		for( int g = 0; g < 3; ++g ) {
			if( left( i, g ) > m_threshhold ) {
				++table( 0, g ) ;
				break ;
			}
		}

		for( int g = 0; g < 3; ++g ) {
			if( right( i, g ) > m_threshhold ) {
				++table( 1, g ) ;
				break ;
			}
		}
	}

	ProductOfIndependentMultinomials alt_model( table ) ;
	Multinomial null_model( table.row(0) + table.row(1) ) ;

	null_model.evaluate_at( null_model.get_MLE() ) ;
	alt_model.evaluate_at( alt_model.get_MLE() ) ;

	double likelihood_ratio_statistic = 2.0 * ( alt_model.get_value_of_function() - null_model.get_value_of_function() ) ;
	double p_value = std::numeric_limits< double >::quiet_NaN() ;
	if( likelihood_ratio_statistic != likelihood_ratio_statistic || likelihood_ratio_statistic < 0.0 ) {
		likelihood_ratio_statistic = std::numeric_limits< double >::quiet_NaN() ;
	}
	else {
		p_value = boost::math::cdf(
			boost::math::complement(
				m_chi_squared,
				likelihood_ratio_statistic
			)
		) ;
	}
	
#if Genotype_FREQUENCY_TEST_CALL_COMPARER_DEBUG
		std::cerr << "Table is:\n" << table << ".\n" ;
		std::cerr << "null_ml is:\n" << null_ml << ".\n" ;
		std::cerr << "alt_ml is:\n" << null_ml << ".\n" ;
		std::cerr << "alt likelihood is " << alt0.get_value_of_function() << " x " << alt1.get_value_of_function() << ".\n" ;
		std::cerr << "null likelihood is " << null_model.get_value_of_function() << ".\n" ;
		std::cerr << "likelihood_ratio_statistic = " << likelihood_ratio_statistic << ".\n" ;
		std::cerr << "p_value = " << p_value << ".\n" ;
#endif

	
	callback( "likelihood_ratio_test_statistic", likelihood_ratio_statistic ) ;
	callback( "pvalue", p_value ) ;
	callback( "lambda=60/rBF", metro::compute_rBF( table, 60 )) ;
}
