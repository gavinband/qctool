#ifndef __LIKELIHOOD_RATIO_TEST_STATISTIC__
#define __LIKELIHOOD_RATIO_TEST_STATISTIC__

#include <cmath>
#include <cassert>
#include "../config.hpp"
#if HAVE_BOOST_MATH
	#include <boost/math/distributions/chi_squared.hpp>
	namespace math = boost::math ;
#else
	namespace math {		
	}
#endif
#include "GenotypeAssayStatistics.hpp"
#include "distributions.hpp"

//
// Implementation of the likelihood ratio test for Hardy-Weinberg equilibrium
//

class LogMaximumLikelihoodForIndependentGenotypes: public GenotypeAssayStatistic
{
	public:
		double calculate_value( GenotypeAssayStatistics const& assay_statistics ) const ;
} ;

class LogMaximumLikelihoodForIndependentGenotypesInHardyWeinberg: public GenotypeAssayStatistic
{
	public:
		double calculate_value( GenotypeAssayStatistics const& assay_statistics ) const ;
} ;

template< typename T1, typename T2 >
class StatisticDifference: public GenotypeAssayStatistic
{
	public:
		double calculate_value( GenotypeAssayStatistics const& assay_statistics ) const {
				return m_statistic1.calculate_value( assay_statistics ) - m_statistic2.calculate_value( assay_statistics ) ;
		}

	private:

		T1 m_statistic1 ;
		T2 m_statistic2 ;
} ;

template< typename T1, typename T2 >
class LikelihoodRatioTestStatistic: public GenotypeAssayStatistic
{
	public:
		LikelihoodRatioTestStatistic()
		#ifdef HAVE_BOOST_MATH
		 :	m_chi_squared_distribution(1)
		#endif
		{};
		
		double calculate_value( GenotypeAssayStatistics const& assay_statistics ) const {
		#ifndef HAVE_BOOST_MATH
			assert( 0 ) ; // Boost.math is required for this test.
			return 0.0 ;
		#else
			double difference = m_statistic_difference.calculate_value( assay_statistics ) ;
			return math::cdf( math::complement( m_chi_squared_distribution, -2.0 * difference )) ;
		#endif
		}			
		
	private:
		StatisticDifference< T1, T2 > m_statistic_difference ;
		#ifdef HAVE_BOOST_MATH
		math::chi_squared_distribution<double> m_chi_squared_distribution ;
		#endif
} ;

typedef LikelihoodRatioTestStatistic< LogMaximumLikelihoodForIndependentGenotypesInHardyWeinberg, LogMaximumLikelihoodForIndependentGenotypes > HardyWeinbergLikelihoodRatioTestStatistic ;

#endif
