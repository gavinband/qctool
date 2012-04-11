
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

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

