#ifndef __GTOOL_GENOTYPEASSAYSTATISTICARITHMETIC__
#define __GTOOL_GENOTYPEASSAYSTATISTICARITHMETIC__

#include "GenotypeAssayStatistics.hpp"



struct StatisticBinaryOp: public GenotypeAssayStatistic 
{	
	StatisticBinaryOp( std::auto_ptr< GenotypeAssayStatistic > statistic1, std::auto_ptr< GenotypeAssayStatistic > statistic2 ) ;
	
	GenotypeAssayStatistic const& first_statistic() const { return *m_1st_statistic ;}
	GenotypeAssayStatistic const& second_statistic() const { return *m_2nd_statistic ;}

private:
	std::auto_ptr< GenotypeAssayStatistic > m_1st_statistic, m_2nd_statistic ;
} ;

struct StatisticRatio: public StatisticBinaryOp 
{	
	StatisticRatio( std::auto_ptr< GenotypeAssayStatistic > statistic1, std::auto_ptr< GenotypeAssayStatistic > statistic2 ) ;
	double calculate_value( GenotypeAssayStatistics const& ) const ;
} ;

struct StatisticProduct: public StatisticBinaryOp 
{	
	StatisticProduct( std::auto_ptr< GenotypeAssayStatistic > statistic1, std::auto_ptr< GenotypeAssayStatistic > statistic2 ) ;
	double calculate_value( GenotypeAssayStatistics const& ) const ;
} ;

struct StatisticSum: public StatisticBinaryOp 
{	
	StatisticSum( std::auto_ptr< GenotypeAssayStatistic > statistic1, std::auto_ptr< GenotypeAssayStatistic > statistic2 ) ;
	double calculate_value( GenotypeAssayStatistics const& ) const ;
} ;

struct StatisticDifference: public StatisticBinaryOp 
{	
	StatisticDifference( std::auto_ptr< GenotypeAssayStatistic > statistic1, std::auto_ptr< GenotypeAssayStatistic > statistic2 ) ;
	double calculate_value( GenotypeAssayStatistics const& ) const ;
} ;



#endif
