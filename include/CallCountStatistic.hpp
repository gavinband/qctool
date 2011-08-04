#ifndef QCTOOL_CALL_COUNT_STATISTIC_HPP
#define QCTOOL_CALL_COUNT_STATISTIC_HPP

#include "GenRowStatistics.hpp"

struct CallCountStatistic: public GenRowSpecificStatistic
{
	CallCountStatistic( std::size_t genotype, double threshhold = 0.9 ) ;
	double calculate_value( GenRowStatistics const& ) const ;
private:
	double const m_threshhold ;
	std::size_t m_genotype ;
} ;

#endif
