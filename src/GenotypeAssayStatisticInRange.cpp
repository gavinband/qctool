
#include <cassert>
#include "RowCondition.hpp"

GenotypeAssayStatisticInInclusiveRange::GenotypeAssayStatisticInInclusiveRange( std::string const& statistic_name, double lower_bound, double upper_bound, double epsilon )
: m_statistic_name( statistic_name ),
	m_lower_bound( lower_bound ),
	m_upper_bound( upper_bound ),
	m_epsilon( epsilon )
{} ;

bool GenotypeAssayStatisticInInclusiveRange::check_if_satisfied( GenotypeAssayStatistics const& row_genotype_statistics ) const {
	double statistic_value = row_genotype_statistics.get_statistic_value< double >( m_statistic_name ) ;
	return ( statistic_value >= ( m_lower_bound - m_epsilon ))
		&& ( statistic_value <= ( m_upper_bound - m_epsilon )) ;
}

void GenotypeAssayStatisticInInclusiveRange::format_to_stream( std::ostream& oStream ) const {
	oStream
		<< m_statistic_name
		<< " in ["
	 	<< m_lower_bound
		<< ","
		<< m_upper_bound
		<< "]" ;
}

GenotypeAssayStatisticInExclusiveRange::GenotypeAssayStatisticInExclusiveRange( std::string const& statistic_name, double lower_bound, double upper_bound, double epsilon )
: m_statistic_name( statistic_name ),
	m_lower_bound( lower_bound ),
	m_upper_bound( upper_bound ),
	m_epsilon( epsilon )
{} ;

bool GenotypeAssayStatisticInExclusiveRange::check_if_satisfied( GenotypeAssayStatistics const& row_genotype_statistics ) const {
	double statistic_value = row_genotype_statistics.get_statistic_value< double >( m_statistic_name ) ;
	return ( statistic_value > ( m_lower_bound - m_epsilon ))
		&& ( statistic_value < ( m_upper_bound - m_epsilon )) ;
}

void GenotypeAssayStatisticInExclusiveRange::format_to_stream( std::ostream& oStream ) const {
	oStream
		<< m_statistic_name
		<< " in ("
	 	<< m_lower_bound
		<< ","
		<< m_upper_bound
		<< ")" ;
}

GenotypeAssayStatisticGreaterThan::GenotypeAssayStatisticGreaterThan( std::string const& statistic_name, double lower_bound, double epsilon )
: m_statistic_name( statistic_name ),
	m_lower_bound( lower_bound ),
	m_epsilon( epsilon )
{} ;

bool GenotypeAssayStatisticGreaterThan::check_if_satisfied( GenotypeAssayStatistics const& row_genotype_statistics ) const {
	double statistic_value = row_genotype_statistics.get_statistic_value< double >( m_statistic_name ) ;
	return ( statistic_value > ( m_lower_bound - m_epsilon )) ;
}

void GenotypeAssayStatisticGreaterThan::format_to_stream( std::ostream& oStream ) const {
	oStream
		<< m_statistic_name
		<< " > "
		<< m_lower_bound ;
}

GenotypeAssayStatisticLessThan::GenotypeAssayStatisticLessThan( std::string const& statistic_name, double upper_bound, double epsilon )
: m_statistic_name( statistic_name ),
	m_upper_bound( upper_bound ),
	m_epsilon( epsilon )
{} ;

bool GenotypeAssayStatisticLessThan::check_if_satisfied( GenotypeAssayStatistics const& row_genotype_statistics ) const {
	double statistic_value = row_genotype_statistics.get_statistic_value< double >( m_statistic_name ) ;
	return ( statistic_value < ( m_upper_bound - m_epsilon )) ;
}

void GenotypeAssayStatisticLessThan::format_to_stream( std::ostream& oStream ) const {
	oStream
		<< m_statistic_name
		<< " < "
		<< m_upper_bound ;
}

