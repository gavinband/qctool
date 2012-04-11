
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef BAYESFACTOR_HPP
#define BAYESFACTOR_HPP

#include <map>
#include <vector>
#include <boost/math/special_functions/beta.hpp>
#include "CaseControlStatistic.hpp"

struct BetaCalculator
{
	double beta( double a, double b ) const {
		return boost::math::beta( a, b ) ;
	}
} ;

struct CachingBetaCalculator
{
	double beta( double a, double b ) const {
		assert( m_beta_cache.size() < 10000000 ) ;

		std::pair< double, double >
			a_and_b( a, b ) ;

		beta_cache_t::const_iterator
			where = m_beta_cache.find( a_and_b ) ;
		double result ;
		if( where == m_beta_cache.end() ) {
			result = boost::math::beta( a, b ) ;
			m_beta_cache[ a_and_b ] = result ;
		}
		else {
			result = where->second ;
		}
		return result ;
	}
	
	std::size_t size() const { return m_beta_cache.size() ; }
	
private:
	typedef std::map< std::pair< double, double >, double > beta_cache_t ;
	mutable beta_cache_t m_beta_cache ;	
} ;


//
// Calculate a bayes factor for case/control data sets,
// as on p.15 of  the 'Supplementary Methods' paper for
// "A new multipoint method for genome-wide association studies...", Marchini et al, 2007.
//
template< typename BetaCalculatorT >
struct SimpleBayesFactor: public CaseControlStatistic, public BetaCalculatorT
{
	SimpleBayesFactor(
		double phi_0 = 1.0,
		double eta_0 = 1.0,
		double phi_1 = 1.0,
		double eta_1 = 1.0,
		double phi_2 = 1.0,
		double eta_2 = 1.0
	)
		: m_phi(3), m_eta(3), m_precomputed_beta_of_phi_and_eta(3)
	{
		m_phi[0] = phi_0 ;
		m_eta[0] = eta_0 ;
		m_phi[1] = phi_1 ;
		m_eta[1] = eta_1 ;
		m_phi[2] = phi_2 ;
		m_eta[2] = eta_2 ;
		
		for( std::size_t i = 0; i < 3u; ++i ) {
			m_precomputed_beta_of_phi_and_eta[i] = beta( m_phi[i], m_eta[i] ) ;
		}
	}
	
protected:
	
	double beta( double a, double b ) const {
		return BetaCalculatorT::beta( a, b ) ;
	}
	
	double calculate_value( case_status_statistic_map_t const& case_status_statistic_map ) const {
		// Assume we have map keys 0.0 and 1.0, and ignore all others.
		case_status_statistic_map_t::const_iterator
			control_i = case_status_statistic_map.find( 0.0 ),
			case_i = case_status_statistic_map.find( 1.0 ) ;

		GenotypeProbabilities const& control_genotype_amounts = control_i->second ;
		GenotypeProbabilities const& case_genotype_amounts = case_i->second ;
		return probability_of_data_given_M4( control_genotype_amounts, case_genotype_amounts )
			/ probability_of_data_given_M0( control_genotype_amounts, case_genotype_amounts ) ;
	}

	double probability_of_data_given_M4( GenotypeProbabilities const& control_genotype_amounts, GenotypeProbabilities const& case_genotype_amounts ) const {
		double result_for_AA_to_0_coding = 0.5;
		double result_for_BB_to_0_coding = 0.5 ;
		for( std::size_t i = 0; i < 3u; ++i ) {
			result_for_AA_to_0_coding *= beta( case_genotype_amounts[i] + m_phi[i], control_genotype_amounts[i] + m_eta[i] ) / m_precomputed_beta_of_phi_and_eta[i] ;
			result_for_BB_to_0_coding *= beta( case_genotype_amounts[2-i] + m_phi[i], control_genotype_amounts[2-i] + m_eta[i] ) / m_precomputed_beta_of_phi_and_eta[i] ;
		}
		return result_for_AA_to_0_coding + result_for_BB_to_0_coding ;
	}

	double probability_of_data_given_M0( GenotypeProbabilities const& control_genotype_amounts, GenotypeProbabilities const& case_genotype_amounts ) const {
		return beta( case_genotype_amounts.sum() + m_phi[0], control_genotype_amounts.sum() + m_eta[0] )
			/ m_precomputed_beta_of_phi_and_eta[0] ;
	}
	
private:
	
	std::vector< double > m_phi ;
	std::vector< double > m_eta ;
	std::vector< double > m_precomputed_beta_of_phi_and_eta ;
} ;



#endif
