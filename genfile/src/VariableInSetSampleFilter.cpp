
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
#include "genfile/VariableInSetSampleFilter.hpp"

namespace genfile {
	VariableInSetSampleFilter::VariableInSetSampleFilter( std::string const& variable ):
		m_variable( variable )
	{}

	void VariableInSetSampleFilter::summarise( std::ostream& o ) const {
		o << m_variable << " in (" ;
		for( std::size_t i = 0; i < std::min( std::size_t( 5 ), m_levels.size() ); ++i ) {
			o << ( i > 0 ? ", " : "" ) << m_levels[i] ;
		}
		if( m_levels.size() > 5 ) {
			o << ", ..." ;
		}
		o << ")" ;
	}

	void VariableInSetSampleFilter::add_level( genfile::VariantEntry value ) {
		if( value.is_int() ) {
			value = double( value.as< int >() ) ;
		}
		if( !( value.is_double() || value.is_string() ) ) {
			throw genfile::BadArgumentError(
				"VariableInSetSampleFilter::add_level()",
				"value=\"" + genfile::string_utils::to_string( value ) + "\"",
				"values must be strings or numbers."
			) ;
		}
		m_levels.push_back( value ) ;
	}

	bool VariableInSetSampleFilter::test(
		genfile::CohortIndividualSource const& source,
		std::size_t i,
		DetailBlock* detail
	) const {
		genfile::VariantEntry value = source.get_entry( i, m_variable ) ;
		bool in_set = false ;
		if( !value.is_missing() ) {
			for( std::size_t i = 0; i < m_levels.size() && !in_set; ++i ) {
				genfile::VariantEntry v = value ;
				genfile::VariantEntry level = m_levels[i] ;
				impl::cast_types_for_comparison( v, level ) ;
				in_set = ( v == level ) ;
			}
		}
		if( detail ) {
			(*detail)( 0, 0 ) = in_set ;
		}
		return in_set ;
	}
}

