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

// Statistic representing the fraction of data that's missing
// in an assay.
struct MissingDataStatistic: public GenotypeAssayStatistic
{
	double calculate_value( GenotypeAssayStatistics const& ) const ;
} ;

struct MissingDataProportionStatistic: public GenotypeAssayStatistic
{
	double calculate_value( GenotypeAssayStatistics const& ) const ;
} ;

struct NumberOfSamplesStatistic: public GenotypeAssayStatistic
{
	double calculate_value( GenotypeAssayStatistics const& ) const ;
} ;

struct AAStatistic: public GenotypeAssayStatistic
{
	double calculate_value( GenotypeAssayStatistics const& ) const ;
} ;

struct ABStatistic: public GenotypeAssayStatistic
{
	double calculate_value( GenotypeAssayStatistics const& ) const ;
} ;

struct BBStatistic: public GenotypeAssayStatistic
{
	double calculate_value( GenotypeAssayStatistics const& ) const ;
} ;

struct ProportionOfAAStatistic: public GenotypeAssayStatistic
{
	double calculate_value( GenotypeAssayStatistics const& ) const ;
} ;


#endif

