#include "GenotypeAssayStatistics.hpp"
#include "AlleleProportions.hpp"
#include "SimpleGenotypeAssayStatistics.hpp"
#include "floating_point_utils.hpp"
#include "GenRow.hpp"
#include <limits>

double MissingDataProportionStatistic::calculate_value( GenotypeAssayStatistics const& statistics ) const {
	if( statistics.number_of_samples() > 0 ) {
		return (statistics.number_of_samples() - statistics.get_genotype_amounts().sum()) / statistics.number_of_samples() ;
	}
	else {
		return std::numeric_limits< double >::quiet_NaN() ;
	}
}

AlleleProportionStatistic::AlleleProportionStatistic( Choice choice ): m_choice( choice ) {}

double AlleleProportionStatistic::calculate_value( GenotypeAssayStatistics const& statistics ) const {
	double result = std::numeric_limits< double >::quiet_NaN() ;
	switch( m_choice ) {
		case minor_allele:
			result = statistics.get_mean_allele_proportions().minor_allele_proportion() ;
			break ;
		case major_allele:
			result = statistics.get_mean_allele_proportions().major_allele_proportion() ;
			break ;
		case first_allele:
			result = statistics.get_mean_allele_proportions().proportion_of_A() ;
			break ;
		case second_allele:
			result = statistics.get_mean_allele_proportions().proportion_of_B() ;
			break ;
	}
	return result ;
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


