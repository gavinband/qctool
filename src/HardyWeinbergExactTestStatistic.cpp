
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cassert>
#include "HardyWeinbergExactTestStatistic.hpp"
#include "floating_point_utils.hpp"
#include "gamma.hpp"


double HardyWeinbergExactTestStatistic::calculate_value( GenotypeAssayStatistics const& statistics ) const {
	double observed_na = std::floor( statistics.get_allele_amounts().minor() ) ;
	double observed_nb = std::floor( statistics.get_allele_amounts().major() ) ;
	double observed_nab = std::floor( statistics.get_genotype_amounts().AB() ) ;

	return calculate_probability_of_HWE( observed_na, observed_nb, observed_nab ) ;
}


double HardyWeinbergExactTestStatistic::calculate_log_of_invariant_part_of_formula1( double const observed_na, double const observed_nb ) const {
	double observed_N = (observed_na + observed_nb) / 2.0 ;

	return
		log_of_factorial( observed_N )
		+ log_of_factorial( observed_na ) 
		+ log_of_factorial( observed_nb ) 
		- log_of_factorial( 2.0 * observed_N ) ;
}


double HardyWeinbergExactTestStatistic::calculate_probability_of_nab_heterozygotes( double const observed_na, double const observed_nb, double const nab, double const log_of_nab_invariant_part ) const {
	double naa = (observed_na - nab) / 2.0 ;
	double nbb = (observed_nb - nab) / 2.0 ;

	assert( naa >= 0.0 ) ;
	assert( nbb >= 0.0 ) ;
	
	double log_of_nab_dependent_part
		= nab * std::log(2) - log_of_factorial( naa ) - log_of_factorial( nab ) - log_of_factorial( nbb ) ;

	return std::exp( log_of_nab_dependent_part + log_of_nab_invariant_part ) ;
}


double HardyWeinbergExactTestStatistic::calculate_probability_of_HWE( double const observed_na, double const observed_nb, double const observed_nab ) const {
	double test_nab ;
	if( check_if_even( observed_na )) {
		test_nab = 0.0 ;
	}
	else {
		test_nab = 1.0 ;
	}

//	std::cout << "N = " << m_observed_N << ".\n" ;
//	std::cout << "#AB = " << m_observed_nab << ", #AA = " << (m_observed_na - m_observed_nab) / 2.0 << ", #AB = " << (m_observed_nb - m_observed_nab) / 2.0 << ".\n";
	
	double total_probability = 0.0 ;
	double log_of_nab_invariant_part = calculate_log_of_invariant_part_of_formula1( observed_na, observed_nb ) ;
	double probability_of_observed_nab = calculate_probability_of_nab_heterozygotes( observed_na, observed_nb, observed_nab, log_of_nab_invariant_part ) ;
	
	while( test_nab <= observed_na ) {
		double probability_of_test_nab = calculate_probability_of_nab_heterozygotes( observed_na, observed_nb, test_nab, log_of_nab_invariant_part ) ;
		if( probability_of_test_nab <= probability_of_observed_nab ) {
			total_probability += probability_of_test_nab ;
		}

		test_nab += 2.0 ;
	}

	return total_probability ;
}
