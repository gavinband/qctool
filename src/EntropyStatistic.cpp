
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include "EntropyStatistic.hpp"
#include "fill_genotype_probabilities.hpp"
#include "rescale_genotype_probabilities.hpp"

// Calculate the Entropy statistic given genotypes.
// Individuals with missing genotypes must have completely missing data.
// The calculation is as follows.
// Let \theta_g be the estimated frequency of genotype g in the data,
//
//             \sum_i p_ig
// \theta_g =  -----------
//                   N
//
// Then this is the average Kulback-Leibler divergence of the genotype distribution from the estimated one,
// divided by the maximum this could be given the genotype frequency estimates.
//
//                                p_ig
//      \sum_i \sum_g p_ig log( -------- )
//                              \theta_g
// KL = ------------------------------
//                                 1
//       sum_g N \theta_g log( --------- )
//                             \theta_g
//
double EntropyStatistic::calculate_value( GenRow const& row ) const {

	if( row.number_of_samples() == 0 ) {
		return 0.0 ;
	}

	std::vector< double > theta( 3 ) ;
	double number_of_non_missing_samples = 0.0 ;
	std::vector< double > e( row.number_of_samples() ) ;
	std::vector< double > f( row.number_of_samples() ) ;

	GenRow::genotype_proportion_const_iterator
		begin_i = row.begin_genotype_proportions(),
		i = begin_i,
		end_i = row.end_genotype_proportions() ;

	// calculate theta
	for( ; i != end_i ; ++i ) {
		for( std::size_t g = 0; g < 3; ++g ) {
			theta[g] += (*i)[ g ] ;
		}
		number_of_non_missing_samples += i->sum() ;
	}

	for( std::size_t g = 0; g < 3; ++g ) {
		theta[g] /= number_of_non_missing_samples ;
	}
	
	double denominator = 0.0 ;
	for( std::size_t g = 0; g < 3; ++g ) {
		denominator += ( number_of_non_missing_samples * theta[g] ) * ( -std::log( theta[g] ) ) ;
	}

	double numerator = 0.0 ;
	for( i = begin_i; i < end_i; ++i ) {
		for( std::size_t g = 0; g < 3; ++g ) {
			numerator += (*i)[ g ] * std::log( (*i)[g] / theta[g] ) ;
		}
	}
	
	//std::cerr << "SNP " << row.RSID() << ": theta = " << theta[0] << ", " << theta[1] << ", " << theta[2] << ", numerator = " << numerator << ", denominator = " << denominator << ".\n" ;

	return numerator / denominator ;
}

double PlainEntropyStatistic::calculate_value( GenRowStatistics const& stats ) const {
	return EntropyStatistic::calculate_value( stats.row() ) ;
}

double FillingEntropyStatistic::calculate_value( GenRowStatistics const& stats ) const {
	InternalStorageGenRow row( stats.row() ) ;
	fill_genotype_probabilities( &row, 0.25, 0.5, 0.25 ) ;
	return EntropyStatistic::calculate_value( row ) ;
}

double ScalingEntropyStatistic::calculate_value( GenRowStatistics const& stats ) const {
	InternalStorageGenRow row( stats.row() ) ;
	rescale_genotype_probabilities( &row, 0.1 ) ;
	return EntropyStatistic::calculate_value( row ) ;
}

