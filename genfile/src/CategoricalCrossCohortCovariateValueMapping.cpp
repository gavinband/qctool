
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
#include "genfile/CategoricalCrossCohortCovariateValueMapping.hpp"
#include "genfile/string_utils.hpp"

namespace genfile {
	CategoricalCrossCohortCovariateValueMapping::CategoricalCrossCohortCovariateValueMapping( std::string const& column_name ):
		LevelCountingCrossCohortCovariateValueMapping( column_name )
	{	
	}

	CohortIndividualSource::Entry CategoricalCrossCohortCovariateValueMapping::get_unmapped_value( Entry const& level ) const {
		int i = level.as< int >() ;
		assert( i > 0 ) ;
		// Adjust because levels are always positive integers
		--i ;
		assert( std::size_t( i ) < histogram().size() ) ;
		Histogram::const_iterator where = histogram().begin() ;
		std::advance( where, std::size_t( i ) ) ;
		return where->first ;
	}
	
	CohortIndividualSource::Entry CategoricalCrossCohortCovariateValueMapping::get_mapped_value( Entry const& entry ) const {
		Histogram::const_iterator where = histogram().find( entry ) ;
		assert( where != histogram().end() ) ;
		// We always return positive integers, so add one.
		return Entry( int( std::distance( histogram().begin(), where ) + 1 ) ) ;
	}
	
	
	std::string CategoricalCrossCohortCovariateValueMapping::get_summary( std::string const& prefix ) const {
		std::ostringstream ostr ;
		ostr << "missing  levels\n"
			<< prefix
			<< std::setw(8) << std::setfill( ' ' ) << std::left << std::fixed << std::setprecision( 5 )
			<< get_number_of_missing_values() ;
		std::size_t count = 0 ;
		for( Histogram::const_iterator i = histogram().begin(); i != histogram().end() && count < 10; ++i, ++count ) {
			ostr << " " << i->first << "(" << i->second << ")";
		}
		if( count < histogram().size() ) {
			ostr << "... (+ " << ( histogram().size() - count ) << " other levels)" ;
		}
		return ostr.str() ;
	}

	std::size_t CategoricalCrossCohortCovariateValueMapping::get_number_of_distinct_mapped_values() const {
		// mapping is one to one
		return get_number_of_distinct_unmapped_values() ;
	}
}
