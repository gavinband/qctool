#ifndef SIMPLEGENOTYPEASSAYSTATISTICS_HPP
#define SIMPLEGENOTYPEASSAYSTATISTICS_HPP

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

struct AAStatistic: public GenotypeAssayStatistic
{
	double calculate_value( GenotypeAssayStatistics const& statistics ) const ;
} ;

struct ABStatistic: public GenotypeAssayStatistic
{
	double calculate_value( GenotypeAssayStatistics const& statistics ) const ;
} ;

struct BBStatistic: public GenotypeAssayStatistic
{
	double calculate_value( GenotypeAssayStatistics const& statistics ) const ;
} ;

#endif

