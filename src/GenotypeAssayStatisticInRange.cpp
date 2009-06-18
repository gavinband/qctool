
#include <cassert>
#include "RowCondition.hpp"

GenotypeAssayStatisticInInclusiveRange::GenotypeAssayStatisticInInclusiveRange( std::string const& statistic_name, double lower_bound, double upper_bound, double epsilon )
: m_statistic_name( statistic_name ),
	m_lower_bound( lower_bound ),
	m_upper_bound( upper_bound ),
	m_epsilon( epsilon )
{} ;

bool GenotypeAssayStatisticInInclusiveRange::check_if_satisfied( GenRow const& genRow, GenotypeAssayStatistics const * row_genotype_statistics_ptr ) const {
	assert( row_genotype_statistics_ptr != 0 ) ;
	
	double statistic_value = row_genotype_statistics_ptr->get_statistic_value< double >( m_statistic_name ) ;
	return ( statistic_value >= ( m_lower_bound - m_epsilon ))
		&& ( statistic_value <= ( m_upper_bound - m_epsilon )) ;
}

GenotypeAssayStatisticInExclusiveRange::GenotypeAssayStatisticInExclusiveRange( std::string const& statistic_name, double lower_bound, double upper_bound, double epsilon )
: m_statistic_name( statistic_name ),
	m_lower_bound( lower_bound ),
	m_upper_bound( upper_bound ),
	m_epsilon( epsilon )
{} ;

bool GenotypeAssayStatisticInExclusiveRange::check_if_satisfied( GenRow const& genRow, GenotypeAssayStatistics const * row_genotype_statistics_ptr ) const {
	assert( row_genotype_statistics_ptr != 0 ) ;
	
	double statistic_value = row_genotype_statistics_ptr->get_statistic_value< double >( m_statistic_name ) ;
	return ( statistic_value > ( m_lower_bound - m_epsilon ))
		&& ( statistic_value < ( m_upper_bound - m_epsilon )) ;
}

