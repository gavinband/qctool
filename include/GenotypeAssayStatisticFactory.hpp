#ifndef __GENROWSTATISTIC_FACTORY_HPP
#define __GENROWSTATISTIC_FACTORY_HPP

#include <string>
#include <utility>
#include "GenotypeAssayStatistics.hpp"

struct GenotypeAssayStatisticFactory
{
	static std::auto_ptr< GenotypeAssayStatistic > create_statistic( std::string statistic_spec ) ;
} ;

#endif

