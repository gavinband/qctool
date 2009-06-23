#ifndef __GENROWSTATISTIC_FACTORY_HPP
#define __GENROWSTATISTIC_FACTORY_HPP

#include <string>
#include <utility>
#include "GenotypeAssayStatistics.hpp"

struct GenRowStatisticFactory
{
	static std::auto_ptr< GenotypeAssayStatistic > create_statistic( std::string statistic_spec ) ;
	static void add_statistics( std::vector< std::string > statistic_specs	, GenotypeAssayStatistics& statistics ) ;
} ;

struct SampleRowStatisticFactory
{
	static std::auto_ptr< GenotypeAssayStatistic > create_statistic( std::string statistic_spec ) ;
	static void add_statistics( std::vector< std::string > statistic_specs	, GenotypeAssayStatistics& statistics ) ;
} ;

struct GenotypeAssayStatisticFactory
{
	static std::auto_ptr< GenotypeAssayStatistic > create_statistic( std::string statistic_spec ) ;
	static void add_statistics( std::vector< std::string > statistic_specs	, GenotypeAssayStatistics& statistics ) ;
} ;

#endif

