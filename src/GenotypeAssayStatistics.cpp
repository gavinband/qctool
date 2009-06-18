#include <cassert>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <numeric>

#include "GenotypeAssayStatistics.hpp"
#include "floating_point_utils.hpp"

template<>
double const& GenotypeAssayStatistic::get_value<double>( GenotypeAssayStatistics const& statistics ) const {
	if( !m_value_is_calculated ) {
		m_value = calculate_value( statistics ) ;
		m_value_is_calculated = true ;
	}
	return m_value ;
}

template<>
std::string const& GenotypeAssayStatistic::get_value<std::string>( GenotypeAssayStatistics const& statistics ) const {
	if( !m_string_value_is_calculated ) {
		m_string_value = calculate_string_value( statistics ) ;
		m_string_value_is_calculated = true ;
	}
	return m_string_value ;
}

std::string GenotypeAssayStatistic::calculate_string_value( GenotypeAssayStatistics const& statistics ) const {
	// default implementation: serialise double value
	double value = get_value<double>( statistics ) ;
	std::ostringstream oStream ;
	oStream << std::fixed << std::setprecision( 2 ) << value ;
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
}

void GenotypeAssayStatistics::reset_statistics() {
	statistics_t::iterator i = m_statistics.begin(),
		end_i = m_statistics.end() ;
	for( ; i != end_i ; ++i ) {
		i->second->reset() ;
	}
}

void GenotypeAssayStatistic::reset() const {
	m_value_is_calculated = false ;
	m_string_value_is_calculated = false ;
}
		
template<>
double const& GenotypeAssayStatistics::get_statistic_value< double >( std::string const& name ) const {
	statistics_t::const_iterator i = m_statistics.find( name ) ;
	assert( i != m_statistics.end() ) ;
	return i->second->get_value<double>( *this ) ;
}

template<>
std::string const& GenotypeAssayStatistics::get_statistic_value< std::string >( std::string const& name ) const {
	statistics_t::const_iterator i = m_statistics.find( name ) ;
	assert( i != m_statistics.end() ) ;
	return i->second->get_value<std::string>( *this ) ;
}

std::ostream& GenotypeAssayStatistics::format_column_headers( std::ostream& aStream ) {
	base_t::format_column_headers( aStream ) ;
	for( statistics_t::const_iterator i = m_statistics.begin(); i != m_statistics.end(); ++i ) {
		aStream << std::setw(8) << std::left << i->first.substr( 0, 8 ) << "  " ;
	}
	return aStream ;
}

std::ostream& GenotypeAssayStatistics::format_statistic_values( std::ostream& aStream ) const {
	base_t::format_statistic_values( aStream ) ;
		for( statistics_t::const_iterator i = m_statistics.begin(); i != m_statistics.end(); ++i ) {
			aStream << std::setw(8) << std::left << get_statistic_value< std::string >( i->first ) << "  ";
		}
		
		return aStream ;
}

double MissingDataStatistic:: calculate_value( GenotypeAssayStatistics const& statistics ) const {
	return statistics.amount_of_missing_data() / statistics.number_of_samples();
}

