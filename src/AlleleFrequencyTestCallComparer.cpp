#include <string>
#include <map>
#include <cmath>
#include <cassert>
#include <boost/function.hpp>
#include <Eigen/Core>
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "integration/NewtonRaphson.hpp"
#include "AlleleFrequencyTestCallComparer.hpp"

namespace {
	struct MultinomialLogLikelihood {
		typedef Eigen::VectorXd Point ;
		typedef Eigen::VectorXd Vector ;
		typedef Eigen::MatrixXd Matrix ;
	
		MultinomialLogLikelihood( Vector const& counts ):
			m_counts( counts ),
			m_parameters( Vector::Zero( counts.size() ) )
		{
		}

		MultinomialLogLikelihood( Vector const& counts, Vector const& parameters ):
			m_counts( counts ),
			m_parameters( parameters )
		{
			assert( parameters.size() == counts.size() ) ;
		}

		void evaluate_at( Vector const& parameters ) {
			assert( parameters.size() == m_counts.size() ) ;
			m_parameters = parameters ;
		}

		double get_value_of_function() const {
			double result = 0.0 ;
			for( int i = 0; i < m_counts.size(); ++i ) {
				result += m_counts(i) * std::log( m_parameters( i )) ;
			}
			return result ;
		}
	
		Vector get_value_of_first_derivative() const {
			Vector result( m_counts.size() ) ;
			for( int i = 0; i < result.size(); ++i ) {
				result(i) = m_counts(i) / m_parameters(i) ;
			}
			return result ;
		}

		Matrix get_value_of_second_derivative() const {
			Matrix result = Matrix::Zero( m_counts.size(), m_counts.size() ) ;
			for( int i = 0; i < result.size(); ++i ) {
				result( i, i ) = -m_counts(i) / ( m_parameters(i) * m_parameters(i) ) ;
			}
			return result ;
		}

	private:
		Vector const m_counts ;
		Vector m_parameters ;
	} ;
}

AlleleFrequencyTestCallComparer::AlleleFrequencyTestCallComparer():
	m_threshhold( 0.9 ),
 	m_chi_squared( 2.0 )
{}

std::map< std::string, genfile::VariantEntry > AlleleFrequencyTestCallComparer::compare(
	genfile::SingleSNPGenotypeProbabilities const& left,
	genfile::SingleSNPGenotypeProbabilities const& right
) const {
	// Make table of counts
	assert( left.size() == right.size() ) ;
	std::size_t const N = left.size() ;
	typedef MultinomialLogLikelihood::Vector Vector ;
	typedef MultinomialLogLikelihood::Matrix Matrix ;
	
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
	
	Vector null_ml = Vector::Zero( 3 ) ;
	Matrix alt_ml = Matrix::Zero( 2, 3 ) ;
	
	for( int g = 0; g < 3; ++g ) {
		null_ml(g) = table.col( g ).sum() / table.sum() ;
		for( int row = 0; row < 2; ++row ) {
			alt_ml(row,g) = table( row, g ) / table.row(row).sum() ;
		}
	}

	MultinomialLogLikelihood
		null_model( table.row(0) + table.row(1) ),
		alt0( table.row(0) ),
		alt1( table.row(1) ) ;

	null_model.evaluate_at( null_ml ) ;
	alt0.evaluate_at( alt_ml.row(0) ) ;
	alt1.evaluate_at( alt_ml.row(1) ) ;

	double likelihood_ratio_statistic = 2.0 * ( alt0.get_value_of_function() + alt1.get_value_of_function() - null_model.get_value_of_function() ) ;
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
	
	std::map< std::string, genfile::VariantEntry > result ;
	result[ "likelihood_ratio_test_statistic" ] = likelihood_ratio_statistic ;
	result[ "pvalue" ] = p_value ;
	return result ;
}
