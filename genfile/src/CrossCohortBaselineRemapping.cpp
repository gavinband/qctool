
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <map>
#include <string>
#include <cassert>
#include <boost/bimap.hpp>
#include "genfile/string_utils/string_utils.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/CrossCohortBaselineRemapping.hpp"
#include "genfile/Error.hpp"

namespace genfile {
	CrossCohortBaselineRemapping::UniquePtr CrossCohortBaselineRemapping::create( CrossCohortCovariateValueMapping::UniquePtr source, CrossCohortBaselineRemapping::Entry const baseline ) {
		return CrossCohortBaselineRemapping::UniquePtr( new CrossCohortBaselineRemapping( source, baseline )) ;
	}
	
	CrossCohortBaselineRemapping::CrossCohortBaselineRemapping( CrossCohortCovariateValueMapping::UniquePtr source, CrossCohortBaselineRemapping::Entry const baseline ):
		m_source( source ),
		m_baseline( baseline )
	{
		assert( m_source.get() ) ;
		recompute_map() ;
	}
	
	void CrossCohortBaselineRemapping::add_source( CohortIndividualSource const& source ) {
		m_source->add_source( source ) ;
		recompute_map() ;
	}
	
	CrossCohortBaselineRemapping::Histogram const& CrossCohortBaselineRemapping::histogram() const {
		return m_source->histogram() ;
	}
	
	void CrossCohortBaselineRemapping::recompute_map() {
		m_map.clear() ;
		Histogram const& histogram = m_source->histogram() ;
		if( histogram.find( m_baseline ) == histogram.end() ) {
			throw BadArgumentError(
				"genfile::CrossCohortBaselineRemapping::recompute_map()",
				"m_baseline=\"" + genfile::string_utils::to_string( m_baseline ) + "\"",
				"Baseline level is not found in the source."
			) ;
		} else {
			m_map.insert( Map::value_type( m_baseline, 1 ) ) ;
			for( int i = 0; i < get_number_of_distinct_mapped_values(); ++i ) {
				Entry const unmapped_value = m_source->get_unmapped_value( i+1 ) ;
				if( unmapped_value == m_baseline ) {
					// baseline already added
					std::cerr << "1:" << unmapped_value << "...\n" << std::flush ;
				} else if( unmapped_value < m_baseline ) {
					m_map.insert( Map::value_type( unmapped_value, i+2 )) ;
					std::cerr << i+2 << ":" << unmapped_value << "...\n" << std::flush ;
				} else {
					std::cerr << i+1 << ":" << unmapped_value << "...\n" << std::flush ;
					m_map.insert( Map::value_type( unmapped_value, i+1 )) ;
				}
			}
		}
	}

	std::size_t CrossCohortBaselineRemapping::get_number_of_unmapped_values() const {
		return m_source->get_number_of_unmapped_values() ;
	}
	std::size_t CrossCohortBaselineRemapping::get_number_of_distinct_unmapped_values() const {
		return m_source->get_number_of_distinct_unmapped_values() ;
	}
	std::size_t CrossCohortBaselineRemapping::get_number_of_distinct_mapped_values() const {
		return m_source->get_number_of_distinct_mapped_values() ;
	}
	std::size_t CrossCohortBaselineRemapping::get_number_of_missing_values() const {
		return m_source->get_number_of_missing_values() ;
	}
	

	CrossCohortBaselineRemapping::Entry CrossCohortBaselineRemapping::get_mapped_value( Entry const& level ) const {
		return m_map.left.find( level )->second ;
	}

	CrossCohortBaselineRemapping::Entry CrossCohortBaselineRemapping::get_unmapped_value( Entry const& level ) const {
		return m_map.right.find( level )->second ;
	}

	// Return a human-readable summary of the mapping and/or its values.
	std::string CrossCohortBaselineRemapping::get_summary( std::string const& prefix ) const {
		return m_source->get_summary( prefix ) + "; baseline=" + genfile::string_utils::to_string( m_baseline ) ;
	}
}
