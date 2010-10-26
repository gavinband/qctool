#include <cmath>
#include "HardyWeinbergEntropyStatistic.hpp"
#include "fill_genotype_probabilities.hpp"
#include "rescale_genotype_probabilities.hpp"

// Calculate the HardyWeinbergEntropy statistic given genotypes.
// Individuals with missing genotypes must have completely missing data.
// The calculation is as follows.
// Let \theta be the estimated allele frequency in the data,
//
//           \sum_i p_i1 + 2p_i2
// \theta =  -------------------
//                   2*N
//
// and let \theta_g be the corresponding estimate of genotype frequency at hardy-weinberg.
// Then this is the average Kulback-Leibler divergence of the genotype distribution from the estimated one (\theta_g),
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
double HardyWeinbergEntropyStatistic::calculate_value( GenRow const& row ) const {

	double result ;

	if( row.number_of_samples() == 0 ) {
		return 0.0 ;
	}

	double number_of_non_missing_samples = 0.0 ;
	double T2 = 0.0 ;
	std::vector< double > e( row.number_of_samples() ) ;
	std::vector< double > f( row.number_of_samples() ) ;

	GenRow::genotype_proportion_const_iterator
		begin_i = row.begin_genotype_proportions(),
		i = begin_i,
		end_i = row.end_genotype_proportions() ;

	double allele_frequency = 0.0 ;
	// calculate allele_frequency
	for( ; i != end_i ; ++i ) {
		allele_frequency += i->AB() + 2 * i->BB() ;
		number_of_non_missing_samples += i->sum() ;
	}
	
	allele_frequency /= 2 * number_of_non_missing_samples ;

	std::vector< double > theta( 3 ) ;
	theta[0] = (1-allele_frequency)*(1-allele_frequency) ;
	theta[1] = 2*allele_frequency*(1-allele_frequency) ;
	theta[2] = allele_frequency*allele_frequency ;

	double denominator = 0.0 ;
	for( std::size_t g = 0; g < 3; ++g ) {
		denominator += number_of_non_missing_samples * theta[g] * ( -std::log( theta[g] ) ) ;
	}

	double numerator = 0.0 ;
	for( i = begin_i; i < end_i; ++i ) {
		for( std::size_t g = 0; g < 3; ++g ) {
			numerator += (*i)[ g ] * std::log( (*i)[g] / theta[g] ) ;
		}
	}
	
	std::cerr << "SNP " << row.RSID() << ": theta = " << theta[0] << ", " << theta[1] << ", " << theta[2] << ", numerator = " << numerator << ", denominator = " << denominator << ".\n" ;

	return numerator / denominator ;
}

double PlainHardyWeinbergEntropyStatistic::calculate_value( GenRowStatistics const& stats ) const {
	return HardyWeinbergEntropyStatistic::calculate_value( stats.row() ) ;
}

double FillingHardyWeinbergEntropyStatistic::calculate_value( GenRowStatistics const& stats ) const {
	InternalStorageGenRow row( stats.row() ) ;
	fill_genotype_probabilities( &row, 0.25, 0.5, 0.25 ) ;
	return HardyWeinbergEntropyStatistic::calculate_value( row ) ;
}

double ScalingHardyWeinbergEntropyStatistic::calculate_value( GenRowStatistics const& stats ) const {
	InternalStorageGenRow row( stats.row() ) ;
	rescale_genotype_probabilities( &row, 0.1 ) ;
	return HardyWeinbergEntropyStatistic::calculate_value( row ) ;
}

