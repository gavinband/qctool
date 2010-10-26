#ifndef ENTROPY_INFORMATION_MEASURE_HPP
#define ENTROPY_INFORMATION_MEASURE_HPP

#include "GenRowStatistics.hpp"

// Implementation of information statistic for GenRows with no missing data.
struct EntropyStatistic: public GenRowSpecificStatistic
{
public:
	virtual ~EntropyStatistic() {} ;
	virtual double calculate_value( GenRow const& ) const ;
} ;

struct PlainEntropyStatistic: public EntropyStatistic
{
	double calculate_value( GenRowStatistics const& ) const ;
} ;

struct FillingEntropyStatistic: public EntropyStatistic
{
	double calculate_value( GenRowStatistics const& ) const ;
} ;

struct ScalingEntropyStatistic: public EntropyStatistic
{
	double calculate_value( GenRowStatistics const& ) const ;
} ;


#endif
