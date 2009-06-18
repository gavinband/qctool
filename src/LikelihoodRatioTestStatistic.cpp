#include <cmath>
#include "GenotypeAssayStatistics.hpp"
#include "LikelihoodRatioTestStatistic.hpp"

double LogMaximumLikelihoodForIndependentGenotypes::calculate_value( GenotypeAssayStatistics const& assay_statistics ) const {
	double sum = assay_statistics.get_genotype_amounts().sum() ;
	double AA = assay_statistics.get_genotype_amounts().AA() ;
	double AB = assay_statistics.get_genotype_amounts().AB() ;
	double BB = assay_statistics.get_genotype_amounts().BB() ;
	double result = 0.0 ;

	// Guard against log(0).  In case one proportion is zero,
	// the correct interpretation is that the term is 0, corresponding to a factor of 1 in the probability.
	if( AA != 0.0 )
		result += AA * std::log( AA / sum ) ;
	if( AB != 0.0 )
		result += AB * std::log( AB / sum ) ;
	if( BB != 0.0 )
		result += BB * std::log( BB / sum ) ;

	return result ;
}


double LogMaximumLikelihoodForIndependentGenotypesInHardyWeinberg::calculate_value( GenotypeAssayStatistics const& assay_statistics ) const {
	double sum = 2.0 * assay_statistics.get_genotype_amounts().sum() ;
	AlleleProportions allele_proportions = assay_statistics.get_allele_amounts() / sum ;
	double A = allele_proportions.A() ;
	double B = allele_proportions.B() ;
	double result = 0.0 ;

	// Guard against log(0).  In case one proportion is zero,
	// the correct interpretation is that the term is 0, corresponding to a factor of 1 in the probability.
	if( A != 0.0 ) {
		result += ( std::log( A ) * 2.0 * assay_statistics.get_genotype_amounts().AA() ) ;
	}
	if( A != 0.0 && B != 0.0 ) {
		result += ( std::log( 2.0 * A * B ) * assay_statistics.get_genotype_amounts().AB() ) ;
	}
	if( B != 0.0 ) {
		result += ( std::log( B ) * 2.0 * assay_statistics.get_genotype_amounts().BB() ) ;
	}
	return result ;
}

