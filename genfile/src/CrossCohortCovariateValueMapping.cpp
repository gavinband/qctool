#include <set>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/CrossCohortCovariateValueMapping.hpp"
#include "genfile/string_utils.hpp"

namespace genfile {
	CrossCohortCovariateValueMapping::UniquePtr CrossCohortCovariateValueMapping::create( CohortIndividualSource::SingleColumnSpec const& column_spec ) {
		if( column_spec.is_continuous() ) {
			if( column_spec.is_phenotype() ) {
				return UniquePtr( new NormalisingCrossCohortCovariateValueMapping( column_spec.name() )) ;
			}
			else {
				return UniquePtr( new ContinuousVariableCrossCohortCovariateValueMapping( column_spec.name() )) ;
			}
		}
		else {
			return UniquePtr( new CategoricalCrossCohortCovariateValueMapping( column_spec.name() )) ;
		}
	}
	
	LevelCountingCrossCohortCovariateValueMapping::LevelCountingCrossCohortCovariateValueMapping(
		std::string const& column_name
	):
		m_column_name( column_name ),
		m_number_of_missing_values( 0u )
	{
	}
	
	void LevelCountingCrossCohortCovariateValueMapping::add_source( CohortIndividualSource const& source ) {
		add_entries_from_source( source, m_column_name ) ;
	}
	
	void LevelCountingCrossCohortCovariateValueMapping::add_entries_from_source( CohortIndividualSource const& source, std::string const& column_name ) {
		for( std::size_t i = 0; i < source.get_number_of_individuals(); ++i ) {
			Entry entry = source.get_entry( i, column_name ) ;
			if( entry.is_missing() ) {
				++m_number_of_missing_values ;
			}
			else {
				++m_entries[ entry ] ;
			}
		}
	}
	
	std::size_t LevelCountingCrossCohortCovariateValueMapping::get_number_of_distinct_mapped_values() const {
		return m_entries.size() ;
	}

	std::size_t LevelCountingCrossCohortCovariateValueMapping::get_number_of_missing_values() const {
		return m_number_of_missing_values ;
	}
	
	CategoricalCrossCohortCovariateValueMapping::CategoricalCrossCohortCovariateValueMapping( std::string const& column_name ):
		LevelCountingCrossCohortCovariateValueMapping( column_name )
	{	
	}

	CohortIndividualSource::Entry CategoricalCrossCohortCovariateValueMapping::get_unmapped_value( Entry const& level ) const {
		int i = level.as< int >() ;
		assert( i > 0 ) ;
		// Adjust because levels are always positive integers
		--i ;
		assert( std::size_t( i ) < entries().size() ) ;
		Entries::const_iterator where = entries().begin() ;
		std::advance( where, std::size_t( i ) ) ;
		return where->first ;
	}
	
	CohortIndividualSource::Entry CategoricalCrossCohortCovariateValueMapping::get_mapped_value( Entry const& entry ) const {
		Entries::const_iterator where = entries().find( entry ) ;
		assert( where != entries().end() ) ;
		// We always return positive integers, so add one.
		return Entry( int( std::distance( entries().begin(), where ) + 1 ) ) ;
	}
	
	
	std::string CategoricalCrossCohortCovariateValueMapping::get_summary( std::string const& prefix ) const {
		std::ostringstream ostr ;
		ostr << "missing  levels\n"
			<< prefix
			<< std::setw(8) << std::setfill( ' ' ) << std::left << std::fixed << std::setprecision( 5 )
			<< get_number_of_missing_values() ;
		for( Entries::const_iterator i = entries().begin(); i != entries().end(); ++i ) {
			ostr << " " << i->first << "(" << i->second << ")";
		}
		return ostr.str() ;
	}

	ContinuousVariableCrossCohortCovariateValueMapping::ContinuousVariableCrossCohortCovariateValueMapping( std::string const& column_name ):
		LevelCountingCrossCohortCovariateValueMapping( column_name )
	{	
	}

	CrossCohortCovariateValueMapping::Entry ContinuousVariableCrossCohortCovariateValueMapping::get_unmapped_value( Entry const& level ) const {
		return level ;
	}

	CrossCohortCovariateValueMapping::Entry ContinuousVariableCrossCohortCovariateValueMapping::get_mapped_value( Entry const& entry ) const {
		return entry ;
	}

	void ContinuousVariableCrossCohortCovariateValueMapping::calculate_mean_and_variance(
		Entries const& entries,
		double* result_mean,
		double* result_variance
	) {
		assert( result_mean != 0 ) ;
		assert( result_variance != 0 ) ;
		// Calculate mean and variance
		double mean = 0.0 ;
		double variance = 0.0 ;
		unsigned int count = 0 ;
		for( Entries::const_iterator i = entries.begin(); i != entries.end(); ++i ) {
			count += i->second ;
			mean += i->second * i->first.as< double >() ;
			variance += i->second * ( i->first.as< double >() * i->first.as< double >() ) ;
		}
		if( count > 0 ) {
			mean /= count ;
		}
		else {
			mean = std::numeric_limits< double >::quiet_NaN() ;
		}

		if( count > 1 ) {
			// Sample variance including Bessel's (N-1) correction.
			variance /= count ;
			variance -= ( mean * mean ) ;
			variance *= double( count ) / double( count - 1 ) ;
		}
		else {
			variance = std::numeric_limits< double >::quiet_NaN() ;
		}
		
		*result_mean = mean ;
		*result_variance = variance ;
	}

	std::string ContinuousVariableCrossCohortCovariateValueMapping::get_summary( std::string const& prefix ) const {
		// Calculate mean and variance
		double mean = 0.0 ;
		double variance = 0.0 ;

		calculate_mean_and_variance( entries(), &mean, &variance ) ;
		
		std::ostringstream ostr ;
			ostr
			<< "missing  min      max      mean     variance\n"
			<< prefix
			<< std::setw(8) << std::setfill( ' ' ) << std::left << std::fixed << std::setprecision( 5 )
			<< get_number_of_missing_values() ;

		if( entries().size() > 0 ) {
			Entries::const_iterator first = entries().begin() ;
			Entries::const_iterator last = entries().end() ;
			--last ;
			ostr
			<< " " << std::setw(8) << std::setfill( ' ' ) << std::left << std::fixed << std::setprecision( 5 )
			<< first->first
			<< " " << std::setw(8) << std::setfill( ' ' ) << std::left << std::fixed << std::setprecision( 5 )
			<< last->first
			<< " " << std::setw(8) << std::setfill( ' ' ) << std::left << std::fixed << std::setprecision( 5 )
			<< mean
			<< " " << std::setw(8) << std::setfill( ' ' ) << std::left << std::fixed << std::setprecision( 5 )
			<< variance ;
		}

		return ostr.str() ;
	}
	
	NormalisingCrossCohortCovariateValueMapping::NormalisingCrossCohortCovariateValueMapping(
		std::string const& column_name
	):
		ContinuousVariableCrossCohortCovariateValueMapping( column_name )
	{}
	
	void NormalisingCrossCohortCovariateValueMapping::add_source( CohortIndividualSource const& source ) {
		ContinuousVariableCrossCohortCovariateValueMapping::add_source( source ) ;
		analyse_values( entries() ) ;
	}
	
	void NormalisingCrossCohortCovariateValueMapping::analyse_values( Entries const& entries ) {
		calculate_mean_and_variance( entries, &m_mean, &m_variance ) ;
		if( m_variance == m_variance ) {
			m_standard_deviation = std::sqrt( m_variance ) ;
		}
		else {
			m_standard_deviation = std::numeric_limits< double >::quiet_NaN() ;
		}
	}
	
	CrossCohortCovariateValueMapping::Entry NormalisingCrossCohortCovariateValueMapping::get_unmapped_value( Entry const& level ) const {
		if( m_variance == m_variance ) { // test for NaN
			return Entry(( level.as< double >() * m_standard_deviation ) + m_mean ) ;
		}
		else {
			return Entry( level.as< double >() + m_mean ) ;
		}
	}
	
	CrossCohortCovariateValueMapping::Entry NormalisingCrossCohortCovariateValueMapping::get_mapped_value( Entry const& entry ) const {
		if( m_variance == m_variance ) { // test for NaN
			return Entry(( entry.as< double >() - m_mean ) / m_standard_deviation ) ;
		}
		else {
			return Entry( entry.as< double >() - m_mean ) ;
		}
	}
	
	std::string NormalisingCrossCohortCovariateValueMapping::get_summary( std::string const& prefix ) const {
		std::string result = "                 "
			+ ContinuousVariableCrossCohortCovariateValueMapping::get_summary( prefix + " (unnormalised): " )
			+ "\n" ;

		Entries::const_iterator first = entries().begin() ;
		Entries::const_iterator last = entries().end() ;
		--last ;
		std::ostringstream ostr ;
		ostr
			<< std::setw(8) << std::setfill( ' ' ) << std::left << std::fixed << std::setprecision( 5 )
			<< get_number_of_missing_values()
			<< " " << std::setw(8) << std::setfill( ' ' ) << std::left << std::fixed << std::setprecision( 5 )
			<< get_mapped_value( first->first )
			<< " " << std::setw(8) << std::setfill( ' ' ) << std::left << std::fixed << std::setprecision( 5 )
			<< get_mapped_value( last->first )
			<< " " << std::setw(8) << std::setfill( ' ' ) << std::left << std::fixed << std::setprecision( 5 )
			<< 0.0
			<< " " << std::setw(8) << std::setfill( ' ' ) << std::left << std::fixed << std::setprecision( 5 )
			<< 1.0 ;
		result += prefix + "   (normalised): " + ostr.str() ;
		return result ;
	}	
}
