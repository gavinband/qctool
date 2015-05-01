
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
#include "genfile/NormalisingCrossCohortCovariateValueMapping.hpp"
#include "genfile/string_utils.hpp"

namespace genfile {
	NormalisingCrossCohortCovariateValueMapping::NormalisingCrossCohortCovariateValueMapping(
		std::string const& column_name
	):
		ContinuousVariableCrossCohortCovariateValueMapping( column_name )
	{}
	
	void NormalisingCrossCohortCovariateValueMapping::add_source( CohortIndividualSource const& source ) {
		ContinuousVariableCrossCohortCovariateValueMapping::add_source( source ) ;
		analyse_values( histogram() ) ;
	}
	
	void NormalisingCrossCohortCovariateValueMapping::analyse_values( Histogram const& histogram ) {
		calculate_mean_and_variance( histogram, &m_mean, &m_variance ) ;
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
