#ifndef INFORMATION_STATISTIC_HPP
#define INFORMATION_STATISTIC_HPP

#include "GenRowStatistics.hpp"

// Implementation of information statistic for GenRows with no missing data.
struct InformationStatistic: public GenRowSpecificStatistic
{
public:
	virtual ~InformationStatistic() {} ;
	virtual double calculate_value( GenRow const& ) const ;
} ;

struct PlainInformationStatistic: public InformationStatistic
{
	double calculate_value( GenRowStatistics const& ) const ;
} ;

struct FillingInformationStatistic: public InformationStatistic
{
	double calculate_value( GenRowStatistics const& ) const ;
} ;

struct ScalingInformationStatistic: public InformationStatistic
{
	double calculate_value( GenRowStatistics const& ) const ;
} ;

#endif

