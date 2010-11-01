#include <cassert>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <numeric>

#include "GenotypeAssayStatistics.hpp"
#include "SimpleGenotypeAssayStatistics.hpp"
#include "floating_point_utils.hpp"
#include "GenRow.hpp"

double MissingDataProportionStatistic::calculate_value( GenotypeAssayStatistics const& statistics ) const {
	if( statistics.number_of_samples() > 0 ) {
		return (statistics.number_of_samples() - statistics.get_genotype_amounts().sum()) / statistics.number_of_samples() ;
	}
	else {
		return std::numeric_limits< double >::quiet_NaN() ;
	}
}

double MinorAlleleProportionStatistic::calculate_value( GenotypeAssayStatistics const& statistics ) const {
	return statistics.get_mean_allele_proportions().minor() ;
}

double HeterozygosityStatistic::calculate_value( GenotypeAssayStatistics const& statistics ) const {
	return statistics.get_genotype_amounts().AB() / statistics.get_genotype_amounts().sum() ;
}

double AAStatistic::calculate_value( GenotypeAssayStatistics const& statistics ) const {
	return statistics.get_genotype_amounts().AA() ;
}

double ABStatistic::calculate_value( GenotypeAssayStatistics const& statistics ) const {
	return statistics.get_genotype_amounts().AB() ;
}

double BBStatistic::calculate_value( GenotypeAssayStatistics const& statistics ) const {
	return statistics.get_genotype_amounts().BB() ;
}


