#include "MachRsquaredStatistic.hpp"
#include "fill_genotype_probabilities.hpp"
#include "rescale_genotype_probabilities.hpp"

// Calculate MachRsquared statistic as per Supplementary S3 of
// "Genotype imputation for genome-wide association studies", Nature Reviews Genetics 11 (2010).
double MachRsquaredStatistic::calculate_value( GenRow const& row ) const {

	double result ;

	if( row.number_of_samples() == 0 ) {
		return 0.0 ;
	}

	double N = 0.0 ;
	double T1 = 0.0 ;
	double T2 = 0.0 ;
	std::vector< double > e( row.number_of_samples() ) ;
	std::vector< double > e_squared( row.number_of_samples() ) ;

	GenRow::genotype_proportion_const_iterator
		i = row.begin_genotype_proportions(),
		end_i = row.end_genotype_proportions() ;

	for( std::size_t index = 0; i != end_i ; ++i, ++index ) {
		N += i->sum() ;
		double e = i->AB() + (2.0 * i->BB()) ;
		T1 += e*e ;
		T2 += e ;
	}
	
	if( N == 0.0 ) {
		return (T2 == 0.0) ? 1.0 : 0.0 ;
	}

	double theta = T2 / (2.0 * N) ;

	if( theta == 0 || theta == 1 ) {
		result = 1.0 ;
	}
	else {
		result = T1 - (T2*T2) ;
		result /= T2 * ( 1 - theta ) ;
	}
	
	return result ;
}

double PlainMachRsquaredStatistic::calculate_value( GenRowStatistics const& stats ) const {
	return MachRsquaredStatistic::calculate_value( stats.row() ) ;
}

double FillingMachRsquaredStatistic::calculate_value( GenRowStatistics const& stats ) const {
	InternalStorageGenRow row( stats.row() ) ;
	fill_genotype_probabilities( &row, 0.25, 0.5, 0.25 ) ;
	return MachRsquaredStatistic::calculate_value( row ) ;
}

double ScalingMachRsquaredStatistic::calculate_value( GenRowStatistics const& stats ) const {
	InternalStorageGenRow row( stats.row() ) ;
	rescale_genotype_probabilities( &row, 0.1 ) ;
	return MachRsquaredStatistic::calculate_value( row ) ;
}

