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
	
	std::size_t LevelCountingCrossCohortCovariateValueMapping::get_number_of_unmapped_values() const {
		Entries::const_iterator i = m_entries.begin(), end_i = m_entries.end() ;
		std::size_t count = 0 ;
		for( ; i != end_i; ++i ) count += i->second ;
		return count ;
	}

	std::size_t LevelCountingCrossCohortCovariateValueMapping::get_number_of_distinct_unmapped_values() const {
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

	std::size_t CategoricalCrossCohortCovariateValueMapping::get_number_of_distinct_mapped_values() const {
		// mapping is one to one
		return get_number_of_distinct_unmapped_values() ;
	}

	ContinuousVariableCrossCohortCovariateValueMapping::ContinuousVariableCrossCohortCovariateValueMapping( std::string const& column_name ):
		LevelCountingCrossCohortCovariateValueMapping( column_name )
	{	
	}

	std::size_t ContinuousVariableCrossCohortCovariateValueMapping::get_number_of_distinct_mapped_values() const {
		// mapping is one to one
		return get_number_of_distinct_unmapped_values() ;
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
		
		std::size_t max_name_width = std::max( std::string( "unnormalised" ).size(), get_mapping_name().size() ) + 3 ;
		
		std::ostringstream ostr ;
			ostr
			<< std::string( max_name_width, ' ' )
			<< " "
			<< "missing  min      max      mean     variance\n"
			<< prefix
			<< std::setw( max_name_width ) << std::right << "(unnormalised):"
			<< " "
			<< std::setw(8) << std::setfill( ' ' ) << std::left << std::fixed << std::setprecision( 4 )
			<< get_number_of_missing_values() ;

		if( entries().size() > 0 ) {
			{
				double mean = 0.0 ;
				double variance = 0.0 ;
				calculate_mean_and_variance( entries(), &mean, &variance ) ;

				Entries::const_iterator first = entries().begin() ;
				Entries::const_iterator last = entries().end() ;
				--last ;
				ostr
				<< " " << std::setw(8) << std::setfill( ' ' ) << std::left << std::fixed << std::setprecision( 4 )
				<< first->first
				<< " " << std::setw(8) << std::setfill( ' ' ) << std::left << std::fixed << std::setprecision( 4 )
				<< last->first
				<< " " << std::setw(8) << std::setfill( ' ' ) << std::left << std::fixed << std::setprecision( 4 )
				<< mean
				<< " " << std::setw(8) << std::setfill( ' ' ) << std::left << std::fixed << std::setprecision( 4 )
				<< variance ;
			}
			{
				// calculate mean and variance of mapped values
				double mean, variance ;
				Entries mapped_values ;
				for( Entries::const_iterator i = entries().begin(); i != entries().end(); ++i ) {
					mapped_values[ get_mapped_value( i->first ) ] = i->second ;
				}
				if( mapped_values != entries() ) {
					calculate_mean_and_variance( mapped_values, &mean, &variance ) ;

					Entries::const_iterator first = entries().begin() ;
					Entries::const_iterator last = entries().end() ;
					--last ;
					ostr
						<< "\n"
						<< prefix
						<< std::setw( max_name_width ) << std::right << ( "(" + get_mapping_name() + "):" ) << " "
						<< std::setw(8) << std::setfill( ' ' ) << std::left << std::fixed << std::setprecision( 4 )
						<< get_number_of_missing_values()
						<< " " << std::setw(8) << std::setfill( ' ' ) << std::left << std::fixed << std::setprecision( 4 )
						<< get_mapped_value( first->first )
						<< " " << std::setw(8) << std::setfill( ' ' ) << std::left << std::fixed << std::setprecision( 4 )
						<< get_mapped_value( last->first )
						<< " " << std::setw(8) << std::setfill( ' ' ) << std::left << std::fixed << std::setprecision( 4 )
						<< mean
						<< " " << std::setw(8) << std::setfill( ' ' ) << std::left << std::fixed << std::setprecision( 4 )
						<< variance ;
				}

				Histogram histogram = get_histogram( mapped_values, 31 ) ;
				if( histogram.size() > 0 ) {
					ostr << "\n"
						<< prefix << std::setw( max_name_width ) << std::right << "(histogram):" << " "
						<< print_histogram( histogram, prefix + std::string( max_name_width + 1, ' ' ), 11, 1 ) ;
				}
			}
		}
		return ostr.str() ;
	}
	
	std::string ContinuousVariableCrossCohortCovariateValueMapping::print_histogram(
		Histogram const& histogram,
		std::string const& prefix,
		std::size_t const height,
		std::size_t const bin_width
	) const {
		// Find highest point of histogram
		Histogram::const_iterator
			highest_i = histogram.begin(),
			end_i = histogram.end() ;
		
		for( Histogram::const_iterator i = histogram.begin(); i != end_i; ++i ) {
			if( highest_i->second < i->second ) {
				highest_i = i ;
			}
		}
		
		// leave
		std::size_t x_axis_size = 2;
		std::size_t graph_height = height - x_axis_size ;
		double vertical_scale = highest_i->second / double( graph_height ) ;
		
		std::ostringstream o ;
		o << "\n" ;

		for( std::size_t row = 0; row < graph_height; ++row ) {
			o << prefix ;
			if( row % 4 == 0 ) {
				o << std::setw( 6 ) << std::right << std::fixed << std::setprecision( 0 ) << ( double( graph_height - row - 0.5 ) * vertical_scale ) << "-|" ;
			}
			else {
				o << std::string( 6, ' ' ) << " |" ;
			}
			for( Histogram::const_iterator i = histogram.begin(); i != end_i; ++i ) {
				if( i->second > ( double( graph_height - row - 0.5 ) * vertical_scale )) {
					o << std::string( bin_width, '*' ) ;
				}
				else {
					o << std::string( bin_width, ' ' ) ; ;
				}
			}
			o << "\n" ;
		}
		
		o << prefix << "       "
			<< "+" + std::string( histogram.size() * bin_width, '-' ) << "\n" ;

		o << prefix << "       "
			<< " "
			<< std::fixed << std::setw( 8 ) << std::left << std::setprecision( 2 ) << histogram.begin()->first.first
			<< std::fixed << std::setw( ( histogram.size() * bin_width ) - 8 )
				<< std::right << std::setprecision( 2 ) << (--histogram.end())->first.second ;

		return o.str() ;
	}
	
	ContinuousVariableCrossCohortCovariateValueMapping::Histogram
	ContinuousVariableCrossCohortCovariateValueMapping::get_histogram(
			Entries const& entries,
			std::size_t const number_of_bins
	) const {
		// divide range into
		Histogram result ;
		if( entries.size() < 10 ) {
			return result;
		}

		double low = entries.begin()->first.as< double >() ;
		double high = (--entries.end())->first.as< double >() ;

		low -= ( high - low ) / ( number_of_bins * 2.0 ) ;
		high += ( high - low ) / ( number_of_bins * 2.0 ) ;
			
		double step = ( high - low ) / double( number_of_bins ) ;
		for( std::size_t i = 0; i < number_of_bins; ++i ) {
			std::size_t bin_count = 0 ;
			double bin_low = low + i * step ;
			double bin_high = ( low + ( i + 1 ) * step ) ;
			for(
				Entries::const_iterator lower_i = entries.lower_bound( bin_low ) ;
				lower_i != entries.end() && lower_i->first.as< double >() < bin_high ;
			 	++lower_i
			) {
				++bin_count ;
			}
			result[ std::make_pair( bin_low, bin_high ) ] = bin_count ;
		}
		return result ;
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
	
	std::size_t NormalisingCrossCohortCovariateValueMapping::get_number_of_distinct_mapped_values() const {
		// mapping is one to one
		return get_number_of_distinct_unmapped_values() ;
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
}
