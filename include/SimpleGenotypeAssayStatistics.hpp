#ifndef __GTOOL_SIMPLEGENOTYPEASSAYSTATISTICS__
#define __GTOOL_SIMPLEGENOTYPEASSAYSTATISTICS__

#include <vector>
#include <iostream>
#include <map>
#include "GenotypeProportions.hpp"
#include "AlleleProportions.hpp"
#include "GToolException.hpp"
#include "GenotypeAssayStatistics.hpp"

// Statistic which just returns 0.0
struct NullStatistic: public GenotypeAssayStatistic
{
	double calculate_value( GenotypeAssayStatistics const& ) const {
		return 0.0 ;
	}
} ;

struct MissingDataProportionStatistic: public GenotypeAssayStatistic
{
	double calculate_value( GenotypeAssayStatistics const& ) const ;
} ;

struct MinorAlleleProportionStatistic: public GenotypeAssayStatistic
{
	double calculate_value( GenotypeAssayStatistics const& statistics ) const ;
} ;

struct HeterozygosityStatistic: public GenotypeAssayStatistic
{
	double calculate_value( GenotypeAssayStatistics const& statistics ) const ;
} ;



#endif

