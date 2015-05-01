
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
#include "genfile/LevelCountingCrossCohortCovariateValueMapping.hpp"
#include "genfile/string_utils.hpp"

namespace genfile {
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
				++m_histogram[ entry ] ;
			}
		}
	}
	
	LevelCountingCrossCohortCovariateValueMapping::Histogram const& LevelCountingCrossCohortCovariateValueMapping::histogram() const {
		return m_histogram ;
	}
	
	std::size_t LevelCountingCrossCohortCovariateValueMapping::get_number_of_unmapped_values() const {
		Histogram::const_iterator i = m_histogram.begin(), end_i = m_histogram.end() ;
		std::size_t count = 0 ;
		for( ; i != end_i; ++i ) count += i->second ;
		return count ;
	}

	std::size_t LevelCountingCrossCohortCovariateValueMapping::get_number_of_distinct_unmapped_values() const {
		return m_histogram.size() ;
	}

	std::size_t LevelCountingCrossCohortCovariateValueMapping::get_number_of_missing_values() const {
		return m_number_of_missing_values ;
	}
}
