
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

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

struct AlleleProportionStatistic: public GenotypeAssayStatistic
{
	enum Choice { minor_allele = 0, major_allele = 1, first_allele = 2, second_allele = 3 } ;
	AlleleProportionStatistic( Choice choice ) ;
	double calculate_value( GenotypeAssayStatistics const& statistics ) const ;
private:
	Choice const m_choice ;
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

