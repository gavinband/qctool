
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <set>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/ContinuousVariableCrossCohortCovariateValueMapping.hpp"
#include "genfile/string_utils.hpp"

namespace genfile {
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
		Histogram const& histogram,
		double* result_mean,
		double* result_variance
	) {
		assert( result_mean != 0 ) ;
		assert( result_variance != 0 ) ;
		// Calculate mean and variance
		double mean = 0.0 ;
		double variance = 0.0 ;
		unsigned int count = 0 ;
		for( Histogram::const_iterator i = histogram.begin(); i != histogram.end(); ++i ) {
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

		if( histogram().size() > 0 ) {
			{
				double mean = 0.0 ;
				double variance = 0.0 ;
				calculate_mean_and_variance( histogram(), &mean, &variance ) ;

				Histogram::const_iterator first = histogram().begin() ;
				Histogram::const_iterator last = histogram().end() ;
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
				Histogram mapped_values ;
				for( Histogram::const_iterator i = histogram().begin(); i != histogram().end(); ++i ) {
					mapped_values[ get_mapped_value( i->first ) ] = i->second ;
				}
				if( mapped_values != histogram() ) {
					calculate_mean_and_variance( mapped_values, &mean, &variance ) ;

					Histogram::const_iterator first = histogram().begin() ;
					Histogram::const_iterator last = histogram().end() ;
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

				NumericalHistogram histogram = get_numerical_histogram( mapped_values, 31 ) ;
				if( histogram.size() > 0 ) {
					ostr << "\n"
						<< prefix << std::setw( max_name_width ) << std::right << "(histogram):" << " "
						<< print_numerical_histogram( histogram, prefix + std::string( max_name_width + 1, ' ' ), 11, 1 ) ;
				}
			}
		}
		return ostr.str() ;
	}
	
	std::string ContinuousVariableCrossCohortCovariateValueMapping::print_numerical_histogram(
		NumericalHistogram const& histogram,
		std::string const& prefix,
		std::size_t const height,
		std::size_t const bin_width
	) const {
		// Find highest point of histogram
		NumericalHistogram::const_iterator
			highest_i = histogram.begin(),
			end_i = histogram.end() ;
		
		for( NumericalHistogram::const_iterator i = histogram.begin(); i != end_i; ++i ) {
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
			for( NumericalHistogram::const_iterator i = histogram.begin(); i != end_i; ++i ) {
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
	
	ContinuousVariableCrossCohortCovariateValueMapping::NumericalHistogram
	ContinuousVariableCrossCohortCovariateValueMapping::get_numerical_histogram(
		Histogram const& histogram,
		std::size_t const number_of_bins
	) const {
		// divide range into
		NumericalHistogram result ;
		if( histogram.size() < 10 ) {
			return result;
		}

		double low = histogram.begin()->first.as< double >() ;
		double high = (--histogram.end())->first.as< double >() ;

		low -= ( high - low ) / ( number_of_bins * 2.0 ) ;
		high += ( high - low ) / ( number_of_bins * 2.0 ) ;
			
		double step = ( high - low ) / double( number_of_bins ) ;
		for( std::size_t i = 0; i < number_of_bins; ++i ) {
			std::size_t bin_count = 0 ;
			double bin_low = low + i * step ;
			double bin_high = ( low + ( i + 1 ) * step ) ;
			for(
				Histogram::const_iterator lower_i = histogram.lower_bound( bin_low ) ;
				lower_i != histogram.end() && lower_i->first.as< double >() < bin_high ;
			 	++lower_i
			) {
				++bin_count ;
			}
			result[ std::make_pair( bin_low, bin_high ) ] = bin_count ;
		}
		return result ;
	}
}
