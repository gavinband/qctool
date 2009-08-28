#include <cassert>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <numeric>

#include "GenotypeAssayStatistics.hpp"
#include "floating_point_utils.hpp"

GenotypeAssayStatistic::GenotypeAssayStatistic()
: m_precision( 4 )
{} ;

template<>
double GenotypeAssayStatistic::get_value<double>( GenotypeAssayStatistics const& statistics ) const {
	return calculate_value( statistics ) ;
}

template<>
std::string GenotypeAssayStatistic::get_value<std::string>( GenotypeAssayStatistics const& statistics ) const {
	return calculate_string_value( statistics ) ;
}

std::string GenotypeAssayStatistic::calculate_string_value( GenotypeAssayStatistics const& statistics ) const {
	// default implementation: serialise double value
	double value = get_value<double>( statistics ) ;
	std::ostringstream oStream ;
	oStream << std::fixed << std::setprecision( m_precision ) << value ;
	return oStream.str() ;
}


GenotypeAssayStatistics::GenotypeAssayStatistics()
{
}

GenotypeAssayStatistics::~GenotypeAssayStatistics() {
	statistics_t::iterator
		statistic_i = m_statistics.begin(),
		end_statistic_i = m_statistics.end() ;
	for( ; statistic_i != end_statistic_i; ++statistic_i ) {
		if( statistic_i->second != 0 ) {
			delete statistic_i->second ;
		}
	}
}

void GenotypeAssayStatistics::add_statistic( std::string const& name, std::auto_ptr< GenotypeAssayStatistic > statistic_ptr ) {
	statistics_t::const_iterator i = m_statistics.find( name ) ;
	if( i != m_statistics.end() ) {
		throw GenotypeAssayStatisticException( "Unable to add statistic \"" + name + "\": a statistic with that name already exists." ) ;
	}
	m_statistics[name] = statistic_ptr.release() ;
	m_statistic_names.push_back( name ) ;
}

void GenotypeAssayStatistics::reset() {
	base_t::reset() ;

	statistics_t::iterator i = m_statistics.begin(),
		end_i = m_statistics.end() ;
	for( ; i != end_i ; ++i ) {
		i->second->reset() ;
	}
}

bool GenotypeAssayStatistics::has_value( std::string const& name ) const {
	statistics_t::const_iterator i = m_statistics.find( name ) ;
	return ( i != m_statistics.end() ) ;
}

double GenotypeAssayStatistics::get_double_value( std::string const& name ) const {
	statistics_t::const_iterator i = m_statistics.find( name ) ;
	assert( i != m_statistics.end() ) ;
	return i->second->get_value<double>( *this ) ;
}

std::string GenotypeAssayStatistics::get_string_value( std::string const& name ) const {
	statistics_t::const_iterator i = m_statistics.find( name ) ;
	assert( i != m_statistics.end() ) ;
	return i->second->get_value<std::string>( *this ) ;
}

std::ostream& GenotypeAssayStatistics::format_column_headers( std::ostream& aStream ) const {
	base_t::format_column_headers( aStream ) ;
	for( std::vector< std::string >::const_iterator i = m_statistic_names.begin(); i != m_statistic_names.end(); ++i ) {
		aStream << std::setw( std::max( std::size_t(8), i->size() )) << std::left << (*i) << "  " ;
	}
	return aStream ;
}

std::ostream& GenotypeAssayStatistics::format_statistic_values( std::ostream& aStream ) const {
	base_t::format_statistic_values( aStream ) ;
	for( std::vector< std::string >::const_iterator i = m_statistic_names.begin(); i != m_statistic_names.end(); ++i ) {
		aStream << std::setw( std::max( std::size_t(8), i->size() )) << std::left << get_value< std::string >( *i ) << "  ";
	}
		
	return aStream ;
}
