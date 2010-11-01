#ifndef MISSINGCALLPROPORTIONSTATISTIC_HPP
#define MISSINGCALLPROPORTIONSTATISTIC_HPP

#include "GenRowStatistics.hpp"

struct MissingCallProportionStatistic: public GenRowSpecificStatistic
{
	MissingCallProportionStatistic( double threshhold = 0.9 ) ;
	double calculate_value( GenRowStatistics const& ) const ;
private:
	double const m_threshhold ;
} ;


#endif
