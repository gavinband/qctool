
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "GenotypeAssayStatisticArithmetic.hpp"

StatisticBinaryOp::StatisticBinaryOp( std::auto_ptr< GenotypeAssayStatistic > statistic1, std::auto_ptr< GenotypeAssayStatistic > statistic2 )
	: m_1st_statistic( statistic1 ),
	m_2nd_statistic( statistic2 )
{}

StatisticRatio::StatisticRatio( std::auto_ptr< GenotypeAssayStatistic > statistic1, std::auto_ptr< GenotypeAssayStatistic > statistic2 )
	: StatisticBinaryOp( statistic1, statistic2 )
	{}

double StatisticRatio::calculate_value( GenotypeAssayStatistics const& stats ) const {
	return first_statistic().get_value<double>( stats ) / second_statistic().get_value<double>( stats ) ;
}

StatisticProduct::StatisticProduct( std::auto_ptr< GenotypeAssayStatistic > statistic1, std::auto_ptr< GenotypeAssayStatistic > statistic2 )
	: StatisticBinaryOp( statistic1, statistic2 )
	{}

double StatisticProduct::calculate_value( GenotypeAssayStatistics const& stats ) const {
	return first_statistic().get_value<double>( stats ) * second_statistic().get_value<double>( stats ) ;
}

StatisticSum::StatisticSum( std::auto_ptr< GenotypeAssayStatistic > statistic1, std::auto_ptr< GenotypeAssayStatistic > statistic2 )
	: StatisticBinaryOp( statistic1, statistic2 )
	{}

double StatisticSum::calculate_value( GenotypeAssayStatistics const& stats ) const {
	return first_statistic().get_value<double>( stats ) + second_statistic().get_value<double>( stats ) ;
}

StatisticDifference::StatisticDifference( std::auto_ptr< GenotypeAssayStatistic > statistic1, std::auto_ptr< GenotypeAssayStatistic > statistic2 )
	: StatisticBinaryOp( statistic1, statistic2 )
	{}

double StatisticDifference::calculate_value( GenotypeAssayStatistics const& stats ) const {
	return first_statistic().get_value<double>( stats ) - second_statistic().get_value<double>( stats ) ;
}
