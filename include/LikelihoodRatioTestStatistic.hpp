#ifndef __LIKELIHOOD_RATIO_TEST_STATISTIC__
#define __LIKELIHOOD_RATIO_TEST_STATISTIC__

#include <cmath>
#include "GenotypeAssayStatistics.hpp"
#include "distributions.hpp"

//
// Implementation of the likelihood ratio test for Hardy-Weinberg equilibrium
//

class MaximumLikelihoodForIndependentGenotypes: public GenotypeAssayStatistic
{
	public:
		double calculate_value( GenotypeAssayStatistics const& assay_statistics ) const ;
} ;

class MaximumLikelihoodForIndependentGenotypesInHardyWeinberg: public GenotypeAssayStatistic
{
	public:
		double calculate_value( GenotypeAssayStatistics const& assay_statistics ) const ;
} ;

template< typename T1, typename T2 >
class StatisticRatio: public GenotypeAssayStatistic
{
	public:
		double calculate_value( GenotypeAssayStatistics const& assay_statistics ) const {
			double denominator = m_statistic2.calculate_value( assay_statistics ) ;
			if( denominator == 0.0 )
				return std::numeric_limits< double >::max() ;
			else
				return m_statistic1.calculate_value( assay_statistics ) / m_statistic2.calculate_value( assay_statistics ) ;
		}

	private:

		T1 m_statistic1 ;
		T2 m_statistic2 ;
} ;

template< typename T1, typename T2 >
class LikelihoodRatioTestStatistic: public GenotypeAssayStatistic
{
	public:
		double calculate_value( GenotypeAssayStatistics const& assay_statistics ) const {
			chi_squared_distribution distribution( 1 ) ;
			return -2.0 * std::log( m_statistic_ratio.calculate_value( assay_statistics ) ) ;
		}
		
	private:
		StatisticRatio< T1, T2 > m_statistic_ratio ;
} ;

typedef StatisticRatio< MaximumLikelihoodForIndependentGenotypesInHardyWeinberg, MaximumLikelihoodForIndependentGenotypes > HardyWeinbergLikelihoodRatioStatistic ;
typedef LikelihoodRatioTestStatistic< MaximumLikelihoodForIndependentGenotypesInHardyWeinberg, MaximumLikelihoodForIndependentGenotypes > HardyWeinbergLikelihoodRatioTestStatistic ;

#endif
