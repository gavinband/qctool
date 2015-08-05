
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <boost/bimap.hpp>
#include "statfile/BuiltInTypeStatSource.hpp"
#include "statfile/BuiltInTypeStatSourceChain.hpp"
#include "statfile/FilteringStatSource.hpp"
#include "genfile/VariantEntry.hpp"

namespace statfile {
	namespace impl {
		// 
		// Homogenise types for value comparison.
		// We do this:
		// * 1. convert any ints to doubles
		// * 2. of one operand is int or double and the other operand is string then we cast the string to a double if possible.
		//
		void cast_types_for_comparison( genfile::VariantEntry& value1, genfile::VariantEntry& value2 ) {
			using genfile::string_utils::to_repr ;
			using genfile::string_utils::StringConversionError ;
			if( value1.is_int() ) {
				value1 = double( value1.as< int >() ) ;
			}
			if( value2.is_int() ) {
				value2 = double( value2.as< int >() ) ;
			}

			if( value1.is_double() && value2.is_string() ) {
				try {
					value2 = to_repr< double >( value2.as< std::string >() ) ;
				} catch( StringConversionError const& ) {
				}
			}
			else if( value2.is_double() && value1.is_string() ) {
				try {
					value1 = to_repr< double >( value1.as< std::string >() ) ;
				} catch( StringConversionError const& ) {
				}
			}
		}
	}
	
	Constraint Constraint::equal_to( genfile::VariantEntry const& value ) {
		return Constraint( value, genfile::VariantEntry(), eEqualTo ) ;
	}

	Constraint Constraint::not_equal_to( genfile::VariantEntry const& value ) {
		return Constraint( value, genfile::VariantEntry(), eNotEqualTo ) ;
	}

	Constraint Constraint::less_than( genfile::VariantEntry const& value ) {
		return Constraint( value, genfile::VariantEntry(), eLessThan ) ;
	}

	Constraint Constraint::less_than_or_equal_to( genfile::VariantEntry const& value ) {
		return Constraint( value, genfile::VariantEntry(), eLessThanOrEqualTo ) ;
	}

	Constraint Constraint::greater_than( genfile::VariantEntry const& value ) {
		return Constraint( value, genfile::VariantEntry(), eGreaterThan ) ;
	}

	Constraint Constraint::greater_than_or_equal_to( genfile::VariantEntry const& value ) {
		return Constraint( value, genfile::VariantEntry(), eGreaterThanOrEqualTo ) ;
	}

	Constraint Constraint::between( genfile::VariantEntry const& value1, genfile::VariantEntry const& value2 ) {
		return Constraint( value1, value2, eBetween ) ;
	}
	
	Constraint::Constraint( Constraint const& other ):
		m_type( other.m_type ),
		m_value1( other.m_value1 ),
		m_value2( other.m_value2 )
	{}
	
	Constraint& Constraint::operator=( Constraint const& other ) {
		m_type = other.m_type ;
		m_value1 = other.m_value1 ;
		m_value2 = other.m_value2 ;
		return *this ;
	}

	bool Constraint::test( genfile::VariantEntry value ) const {
		genfile::VariantEntry refValue = m_value1 ;
		genfile::VariantEntry refValue2 = m_value2 ;
		impl::cast_types_for_comparison( value, refValue ) ;
		switch( m_type ) {
			case eEqualTo:
				return (value == refValue) ;
				break ;
			case eNotEqualTo:
				return (value != refValue) ;
				break ;
			case eLessThan:
				return (value < refValue) ;
				break ;
			case eLessThanOrEqualTo:
				return (value <= refValue) ;
				break ;
			case eGreaterThan:
				return (value > refValue) ;
				break ;
			case eGreaterThanOrEqualTo:
				return (value >= refValue) ;
				break ;
			case eBetween:
				impl::cast_types_for_comparison( value, refValue2 ) ;
				return (value >= refValue) && (value <= refValue2) ;
				break ;
			default:
				assert(0) ;
		}
	}

	Constraint::Constraint( genfile::VariantEntry const& value1, genfile::VariantEntry const& value2, Type type ):
		m_type( type ),
		m_value1( value1 ),
		m_value2( value2 )
	{}

	FilteringStatSource::FilteringStatSource(
		BuiltInTypeStatSource::UniquePtr source,
		std::string const column,
		Constraint const& constraint
	):
		m_source( source ),
		m_constraint( constraint )
	{
		m_columns.insert( ColumnMap::value_type( column, m_source->index_of_column( column )) ) ;
		move_to_next_matching_row() ;
	}

	FilteringStatSource::operator bool() const {
		return (*m_source) ;
	}

	void FilteringStatSource::read_value( int32_t& value ) {
		(*m_source) >> value ;
	}
	void FilteringStatSource::read_value( uint32_t& value ) {
		(*m_source) >> value ;
	}
	void FilteringStatSource::read_value( std::string& value ) {
		(*m_source) >> value ;
	}
	void FilteringStatSource::read_value( double& value ) {
		(*m_source) >> value ;
	}
	void FilteringStatSource::ignore_value() {
		(*m_source) >> ignore() ;
	}
	void FilteringStatSource::ignore_all() {
		(*m_source) >> statfile::ignore_all() ;
	}
	void FilteringStatSource::end_row() {
		(*m_source) >> statfile::end_row() ;
	} ;
	void FilteringStatSource::restart_row() {
		(*m_source) >> statfile::restart_row() ;
	}
	void FilteringStatSource::move_to_next_row_impl() {
		move_to_next_matching_row() ;
	}

	void FilteringStatSource::move_to_next_matching_row() {
		while( (*m_source) && !row_matches() ) {
			(*m_source) >> statfile::ignore_all() ;
		}
		if( *m_source ) {
			(*m_source) >> statfile::restart_row() ;
		}
	}

	bool FilteringStatSource::row_matches() {
		assert( m_source->current_column() == 0 ) ;
		boost::bimap< std::string, std::size_t >::left_const_iterator i = m_columns.left.begin() ;
		boost::bimap< std::string, std::size_t >::left_const_iterator const end_i = m_columns.left.end() ;
		std::string value ;
		for( ; i != end_i; ++i ) {
			(*m_source) >> ignore( i->second - m_source->current_column() ) >> value ;
			if( !m_constraint.test( value ) ) {
				return false ;
			}
		}
		return true ;
	}
}
