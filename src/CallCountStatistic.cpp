#include <cassert>
#include "CallCountStatistic.hpp"

CallCountStatistic::CallCountStatistic( std::size_t genotype, double threshhold ):
	m_threshhold( threshhold ),
	m_genotype( genotype )
{
	assert( genotype < 3 ) ;
}

double CallCountStatistic::calculate_value( GenRowStatistics const& stats ) const {
	GenRow const& row = stats.row() ;
	GenRow::genotype_proportion_const_iterator
		begin = row.begin_genotype_proportions(),
		end = row.end_genotype_proportions() ;
	double count = 0 ;
	for( ; begin != end; ++begin ) {
		count += ( (*begin)[ m_genotype ] > m_threshhold ) ? 1 : 0 ;
	}
	return double( count ) ;
}

