#include <cmath>
#include "GenotypeAssayStatistics.hpp"
#include "LikelihoodRatioTestStatistic.hpp"

double MaximumLikelihoodForIndependentGenotypes::calculate_value( GenotypeAssayStatistics const& assay_statistics ) const {
	double sum = assay_statistics.get_genotype_amounts().sum() ;
	
	return
		  std::pow( assay_statistics.get_genotype_amounts().AA() / sum, assay_statistics.get_genotype_amounts().AA() )
		* std::pow( assay_statistics.get_genotype_amounts().AB() / sum, assay_statistics.get_genotype_amounts().AB() )
		* std::pow( assay_statistics.get_genotype_amounts().BB() / sum, assay_statistics.get_genotype_amounts().BB() ) ;
}


double MaximumLikelihoodForIndependentGenotypesInHardyWeinberg::calculate_value( GenotypeAssayStatistics const& assay_statistics ) const {
	double sum = 2.0 * assay_statistics.get_genotype_amounts().sum() ;
	AlleleProportions allele_proportions = assay_statistics.get_allele_amounts() / sum ;
	
	return
		std::pow( allele_proportions.A(), 2.0 * assay_statistics.get_genotype_amounts().AA() )
		* std::pow( 2.0 * allele_proportions.A() * allele_proportions.B(), assay_statistics.get_genotype_amounts().AB() )
		* std::pow( allele_proportions.B(), 2.0 * assay_statistics.get_genotype_amounts().BB() ) ;
}

