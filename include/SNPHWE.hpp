
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


#ifndef __GTOOL__SNPHWE_HPP__
#define __GTOOL__SNPHWE_HPP__

#include <cmath>
#include "GenotypeAssayStatistics.hpp"

double SNPHWE(int obs_hets, int obs_hom1, int obs_hom2) ;

struct SNPHWEStatistic: public GenotypeAssayStatistic
{
	double calculate_value( GenotypeAssayStatistics const& statistics ) const {
		if( statistics.get_rounded_genotype_amounts().sum() > 0.0 ) {
			return SNPHWE( statistics.get_rounded_genotype_amounts().AB(), statistics.get_rounded_genotype_amounts().AA(), statistics.get_rounded_genotype_amounts().BB() ) ;
		}
		else {
			return 1.0 ;
		}
	}
} ;

struct MinusLog10SNPHWEStatistic: public GenotypeAssayStatistic
{
	double calculate_value( GenotypeAssayStatistics const& statistics ) const {
		return -std::log10( m_snp_hwe_statistic.calculate_value( statistics )) ;
	}
	
private:
	
	SNPHWEStatistic m_snp_hwe_statistic ;
	
} ;

#endif
