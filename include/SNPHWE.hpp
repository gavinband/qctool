
#ifndef __GTOOL__SNPHWE_HPP__
#define __GTOOL__SNPHWE_HPP__

#include "GenotypeAssayStatistics.hpp"

double SNPHWE(int obs_hets, int obs_hom1, int obs_hom2) ;

struct SNPHWEStatistic: public GenotypeAssayStatistic
{
	double calculate_value( GenotypeAssayStatistics const& statistics ) const {
		return SNPHWE( statistics.get_genotype_amounts().AB(), statistics.get_genotype_amounts().AA(), statistics.get_genotype_amounts().BB() ) ;
	}
} ;

#endif
