#ifndef MACH_R_SQUARED_STATISTIC_HPP
#define MACH_R_SQUARED_STATISTIC_HPP

#include "GenRowStatistics.hpp"

// Implementation of information statistic for GenRows with no missing data.
struct MachRsquaredStatistic: public GenRowSpecificStatistic
{
public:
	virtual ~MachRsquaredStatistic() {} ;
	virtual double calculate_value( GenRow const& ) const ;
} ;

struct PlainMachRsquaredStatistic: public MachRsquaredStatistic
{
	double calculate_value( GenRowStatistics const& ) const ;
} ;

struct FillingMachRsquaredStatistic: public MachRsquaredStatistic
{
	double calculate_value( GenRowStatistics const& ) const ;
} ;

struct ScalingMachRsquaredStatistic: public MachRsquaredStatistic
{
	double calculate_value( GenRowStatistics const& ) const ;
} ;

#endif

