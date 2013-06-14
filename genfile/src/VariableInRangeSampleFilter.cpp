
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


#include <memory>
#include <vector>
#include <set>
#include <boost/ptr_container/ptr_vector.hpp>
#include "genfile/VariantEntry.hpp"
#include "genfile/SampleFilter.hpp"
#include "genfile/impl/cast_types_for_comparison.hpp"
#include "genfile/VariableInRangeSampleFilter.hpp"

namespace genfile {
	VariableInRangeSampleFilter::VariableInRangeSampleFilter( std::string const& variable, genfile::VariantEntry lower_bound, genfile::VariantEntry upper_bound ):
		m_variable( variable ),
		m_lower_bound( lower_bound ),
		m_upper_bound( upper_bound )
	{
		// Always compare numbers as doubles.
		if( m_lower_bound.is_int() ) {
			m_lower_bound = double( m_lower_bound.as< int >() ) ;
		}
		if( m_upper_bound.is_int() ) {
			m_upper_bound = double( m_upper_bound.as< int >() ) ;
		}
	
		if(
			m_lower_bound.is_string() && m_upper_bound.is_string() 
			|| m_lower_bound.is_double() && m_upper_bound.is_double()
		) {
			// ok
		} else {
			throw genfile::BadArgumentError(
				"VariableInRangeSampleFilter::VariableInRangeSampleFilter()",
				"upper_bound=\"" + genfile::string_utils::to_string( upper_bound ) + "\"",
				"Lower and upper bound must both be strings or numbers."
			) ;
		}	
	}

	bool VariableInRangeSampleFilter::test( genfile::CohortIndividualSource const& source, std::size_t i, DetailBlock* detail ) const {
		genfile::VariantEntry value = source.get_entry( i, variable() ) ;
		bool in_range = false ;
		if( !value.is_missing() ) {
			in_range = compare_impl( value, m_lower_bound, m_upper_bound ) ;
		}
		if( detail ) {
			(*detail)( 0, 0 ) = in_range ;
		}
		return in_range ;
	}

	VariableInInclusiveRangeSampleFilter::VariableInInclusiveRangeSampleFilter( std::string const& variable, genfile::VariantEntry lower_bound, genfile::VariantEntry upper_bound ):
		VariableInRangeSampleFilter( variable, lower_bound, upper_bound )
	{}

	void VariableInInclusiveRangeSampleFilter::summarise( std::ostream& oStream ) const {
		oStream << "VariableInInclusiveRangeSampleFilter( " << lower_bound() << ", " << upper_bound() << " )" ;
	}

	bool VariableInInclusiveRangeSampleFilter::compare_impl( genfile::VariantEntry const& value, genfile::VariantEntry const& lower, genfile::VariantEntry const& upper ) const {
		bool result = true ;
		genfile::VariantEntry v = value, a = lower, b = upper ;
		impl::cast_types_for_comparison( v, a ) ;
		result &= ( v >= a ) ;
		v = value ;
		impl::cast_types_for_comparison( v, b ) ;
		result &= ( v <= b ) ;
		return result ;
	}

	VariableInExclusiveRangeSampleFilter::VariableInExclusiveRangeSampleFilter( std::string const& variable, genfile::VariantEntry lower_bound, genfile::VariantEntry upper_bound ):
		VariableInRangeSampleFilter( variable, lower_bound, upper_bound )
	{}

	void VariableInExclusiveRangeSampleFilter::summarise( std::ostream& oStream ) const {
		oStream << "VariableInExclusiveRangeSampleFilter( " << lower_bound() << ", " << upper_bound() << " )" ;
	}

	bool VariableInExclusiveRangeSampleFilter::compare_impl( genfile::VariantEntry const& value, genfile::VariantEntry const& lower, genfile::VariantEntry const& upper ) const {
		bool result = true ;
		genfile::VariantEntry v = value, a = lower, b = upper ;
		impl::cast_types_for_comparison( v, a ) ;
		result &= ( v > a ) ;
		v = value ;
		impl::cast_types_for_comparison( v, b ) ;
		result &= ( v < b ) ;
		return result ;
	}
}
