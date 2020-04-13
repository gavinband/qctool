
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
		o << m_variable << " IN (" ;
		Levels::const_iterator i = m_levels.begin() ;
		for( std::size_t count = 0; count < std::min( std::size_t( 5 ), m_levels.size() ); ++count, ++i ) {
			o << ( count > 0 ? ", " : "" ) << *i ;
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
		if( !( value.is_double() || value.is_string() || value.is_missing() ) ) {
			throw genfile::BadArgumentError(
				"VariableInSetSampleFilter::add_level()",
				"value=\"" + genfile::string_utils::to_string( value ) + "\"",
				"values must be strings or numbers."
			) ;
		}
		//m_levels.push_back( value ) ;
		m_levels.insert( value ) ;
		if( value.is_string() ) {
			try {
				value = genfile::string_utils::to_repr< double >( value.as< std::string >() ) ;
				m_levels.insert( value ) ;
			}
			catch( genfile::string_utils::StringConversionError const& e ) {
				// nothing to do
			}
		}
	}

	bool VariableInSetSampleFilter::test(
		genfile::CohortIndividualSource const& source,
		std::size_t i,
		DetailBlock* detail
	) const {
		genfile::VariantEntry value = source.get_entry( i, m_variable ) ;
		// integers are compared as doubles
		if( value.is_int() ) {
			value = double( value.as< int >() ) ;
		}
		return m_levels.find( value ) != m_levels.end() ;
#if 0
		bool in_set = false ;
		if( !value.is_missing() ) {
			for( std::size_t i = 0; i < m_levels.size() && !in_set; ++i ) {
				genfile::VariantEntry v = value ;
				genfile::VariantEntry level = m_levels[i] ;
				impl::cast_types_for_comparison( v, level ) ;
				in_set = ( v == level ) ;
			}
		} else {
			for( std::size_t i = 0; i < m_levels.size() && !in_set; ++i ) {
				in_set = m_levels[i].is_missing() ;
			}
		}
		if( detail ) {
			(*detail)( 0, 0 ) = in_set ;
		}
		return in_set ;
#endif
	}

	VariableNotInSetSampleFilter::VariableNotInSetSampleFilter( std::string const& variable ):
		m_inverse_filter( variable )
	{}

	void VariableNotInSetSampleFilter::summarise( std::ostream& o ) const {
		o << m_inverse_filter.m_variable << " NOT IN (" ;
		VariableInSetSampleFilter::Levels::const_iterator i = m_inverse_filter.m_levels.begin() ;
		VariableInSetSampleFilter::Levels::const_iterator end_i = m_inverse_filter.m_levels.end() ;
		for( std::size_t count = 0; count < std::min( std::size_t( 5 ), m_inverse_filter.m_levels.size() ); ++count, ++i ) {
			o << ( count > 0 ? ", " : "" ) << *i ;
		}
		if( m_inverse_filter.m_levels.size() > 5 ) {
			o << ", ..." ;
		}
		o << ")" ;
	}

	void VariableNotInSetSampleFilter::add_level( genfile::VariantEntry value ) {
		m_inverse_filter.add_level( value ) ;
	}

	bool VariableNotInSetSampleFilter::test(
		genfile::CohortIndividualSource const& source,
		std::size_t i,
		DetailBlock* detail
	) const {
		bool result = !m_inverse_filter.test( source, i, 0 ) ;
		if( detail ) {
			(*detail)( 0, 0 ) = result ;
		}
		return result ;
	}
}

