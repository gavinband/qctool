#include <cassert>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <numeric>

#include "GenotypeAssayStatistics.hpp"
#include "SimpleGenotypeAssayStatistics.hpp"
#include "floating_point_utils.hpp"

double MissingDataStatistic::calculate_value( GenotypeAssayStatistics const& statistics ) const {
	return statistics.number_of_samples() - statistics.get_genotype_amounts().sum() ;
}

double MissingDataProportionStatistic::calculate_value( GenotypeAssayStatistics const& statistics ) const {
	return (statistics.number_of_samples() - statistics.get_genotype_amounts().sum()) / statistics.get_genotype_amounts().sum() ;
}

double NumberOfSamplesStatistic::calculate_value( GenotypeAssayStatistics const& statistics ) const {
	return statistics.number_of_samples() ;
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

double MinorAlleleProportionStatistic::calculate_value( GenotypeAssayStatistics const& statistics ) const {
	return statistics.get_mean_allele_proportions().minor() ;
}

double HeterozygosityStatistic::calculate_value( GenotypeAssayStatistics const& statistics ) const {
	return statistics.get_genotype_amounts().AB() / statistics.get_genotype_amounts().sum() ;
}
