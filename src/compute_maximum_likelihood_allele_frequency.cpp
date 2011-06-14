#include <limits>
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "compute_maximum_likelihood_allele_frequency.hpp"

double compute_maximum_likelihood_allele_frequency( genfile::SingleSNPGenotypeProbabilities const& genotypes ) {
	// Under a model in which alleles are drawn randomly from a population of haplotypes,
	// with given allele frequency, return the maximum likelihood allele frequency.
	// Note that this deals with null genotype calls (= missing data) by dividing by the total amount
	// of data (not the total number of samples).
	double allele_count = 0.0 ;
	double data_count = 0.0 ;
	for( std::size_t i = 0; i < genotypes.size(); ++i ) {
		allele_count += genotypes( i, 1 ) + 2.0 * genotypes( i, 2 ) ; 
		data_count += genotypes.sum( i ) ;
	}

	if( data_count > 0.0 ) {
		return allele_count / (2.0 * data_count ) ;
	}
	else {
		return std::numeric_limits< double >::quiet_NaN() ;
	}
}

