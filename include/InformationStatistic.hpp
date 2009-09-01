#ifndef INFORMATION_STATISTIC_HPP
#define INFORMATION_STATISTIC_HPP

#include "GenRowStatistics.hpp"

// Base class for individual statistics
struct InformationStatistic: public GenRowSpecificStatistic
{
	double calculate_value( GenRowStatistics const& ) const ;
} ;

#endif

