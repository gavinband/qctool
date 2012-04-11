
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <numeric>
#include <boost/math/distributions/normal.hpp>
#include "genfile/CrossCohortCovariateValueMapping.hpp"
#include "genfile/QuantileNormalisingCrossCohortCovariateValueMapping.hpp"

namespace genfile {
	QuantileNormalisingCrossCohortCovariateValueMapping::QuantileNormalisingCrossCohortCovariateValueMapping(
		std::string const& column_name
	):
		ContinuousVariableCrossCohortCovariateValueMapping( column_name )
	{}
	
	void QuantileNormalisingCrossCohortCovariateValueMapping::add_source( CohortIndividualSource const& source ) {
		ContinuousVariableCrossCohortCovariateValueMapping::add_source( source ) ;
		analyse_values( entries() ) ;
	}
	
	void QuantileNormalisingCrossCohortCovariateValueMapping::analyse_values( Entries const& entries ) {
		// Generate Normal quantiles, one per entry.
		std::size_t const N = get_number_of_unmapped_values() ;
		std::vector< double > normal_quantiles( N ) ;
		for( std::size_t i = 0; i < N; ++i ) {
			normal_quantiles[i] = boost::math::quantile( m_normal_distribution, double( i+1 ) / double( N + 1 ) ) ;
		}

		Entries mapped_entries ;

		Mapping mapping, reverse_mapping ;
		std::size_t quantile_i = 0 ;
		for(
			Entries::const_iterator i = entries.begin();
			i != entries.end();
			quantile_i += (i++)->second
		) {
			assert( ( quantile_i + i->second ) <= N ) ;
			// value for this entry is mean of the next m, where m is the multiplicity.
			double mean = std::accumulate( normal_quantiles.begin() + quantile_i, normal_quantiles.begin() + quantile_i + i->second, 0.0 ) ;
			mean /= double( i->second ) ;
			mapping[ i->first.as< double >() ] = mean ;
			reverse_mapping[ mean ] = i->first.as< double >() ;
			mapped_entries[ Entry( mean ) ] = i->second ;
		}
		assert( mapping.size() == entries.size() ) ;
		assert( reverse_mapping.size() == mapping.size() ) ;
		
		m_mapping = mapping ;
		m_reverse_mapping = reverse_mapping ;
	}
	
	std::size_t QuantileNormalisingCrossCohortCovariateValueMapping::get_number_of_distinct_mapped_values() const {
		return get_number_of_distinct_unmapped_values() ;
	}
	
	CrossCohortCovariateValueMapping::Entry QuantileNormalisingCrossCohortCovariateValueMapping::get_unmapped_value( Entry const& level ) const {
		Mapping::const_iterator where = m_reverse_mapping.find( level.as< double >() ) ;
		return Entry( where->second ) ;
	}
	
	CrossCohortCovariateValueMapping::Entry QuantileNormalisingCrossCohortCovariateValueMapping::get_mapped_value( Entry const& entry ) const {
		Mapping::const_iterator where = m_mapping.find( entry.as< double >() ) ;
		return Entry( where->second ) ;
	}
}

