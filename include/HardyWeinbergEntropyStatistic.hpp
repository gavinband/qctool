#ifndef HARDY_WEINBERG_ENTROPY_INFORMATION_MEASURE_HPP
#define HARDY_WEINBERG_ENTROPY_INFORMATION_MEASURE_HPP

#include "GenRowStatistics.hpp"

// Implementation of information statistic for GenRows with no missing data.
struct HardyWeinbergEntropyStatistic: public GenRowSpecificStatistic
{
public:
	virtual ~HardyWeinbergEntropyStatistic() {} ;
	virtual double calculate_value( GenRow const& ) const ;
} ;

struct PlainHardyWeinbergEntropyStatistic: public HardyWeinbergEntropyStatistic
{
	double calculate_value( GenRowStatistics const& ) const ;
} ;

struct FillingHardyWeinbergEntropyStatistic: public HardyWeinbergEntropyStatistic
{
	double calculate_value( GenRowStatistics const& ) const ;
} ;

struct ScalingHardyWeinbergEntropyStatistic: public HardyWeinbergEntropyStatistic
{
	double calculate_value( GenRowStatistics const& ) const ;
} ;


#endif
