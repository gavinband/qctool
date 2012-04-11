
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef INFORMATION_STATISTIC_HPP
#define INFORMATION_STATISTIC_HPP

#include "GenRowStatistics.hpp"

// Implementation of information statistic for GenRows with no missing data.
struct OldInformationStatistic: public GenRowSpecificStatistic
{
public:
	virtual ~OldInformationStatistic() {} ;
	virtual double calculate_value( GenRow const& ) const ;
	double calculate_value( GenRowStatistics const& ) const ;
	
} ;

struct InformationStatistic: public GenRowSpecificStatistic
{
public:
	virtual ~InformationStatistic() {} ;
	virtual double calculate_value( GenRow const& ) const ;
	double calculate_value( GenRowStatistics const& ) const ;
} ;

struct GavinsInformationStatistic: public GenRowSpecificStatistic
{
public:
	virtual ~GavinsInformationStatistic() {} ;
	virtual double calculate_value( GenRow const& ) const ;
	double calculate_value( GenRowStatistics const& ) const ;
} ;

#endif

