
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <boost/variant.hpp>
#include <boost/noncopyable.hpp>
#include "genfile/VariantEntry.hpp"
#include "genfile/endianness_utils.hpp"

namespace genfile {
	VariantEntry::VariantEntry():
		m_entrydata( MissingValue() )
	{}
	
	VariantEntry::VariantEntry( int const value ):
			m_entrydata( Integer( value ))
	{}
	
	VariantEntry& VariantEntry::operator=( int const value ) {
		this->m_entrydata = Integer( value ) ;
		return *this ;
	}
	
	bool VariantEntry::is_missing() const {
		return m_entrydata.which() == eMissing ;
	}

	bool VariantEntry::is_string() const {
		return m_entrydata.which() == eString ;
	}

	bool VariantEntry::is_int() const {
		return m_entrydata.which() == eInteger ;
	}

	bool VariantEntry::is_double() const {
		return m_entrydata.which() == eDouble ;
	}

	bool VariantEntry::is_chromosome() const {
		return m_entrydata.which() == eChromosome ;
	}

	bool VariantEntry::is_position() const {
		return m_entrydata.which() == eGenomePosition ;
	}
	
	template<> double VariantEntry::as() const {
		if( Integer const* v = boost::get< Integer >( &m_entrydata )) {
			return double( *v ) ;
		}
		else {
			return boost::get< double >( m_entrydata ) ;
		} 
	}

	template<> int VariantEntry::as() const {
		Integer integer = boost::get< Integer >( m_entrydata ) ;
		int result( integer ) ;
		assert( result == integer ) ;
		return result ;
	}
	
	bool VariantEntry::operator==( VariantEntry const& rhs ) const {
		return m_entrydata == rhs.m_entrydata ;
	}

	bool VariantEntry::operator!=( VariantEntry const& rhs ) const {
		return !( m_entrydata == rhs.m_entrydata ) ;
	}
	
	bool VariantEntry::operator<( VariantEntry const& rhs ) const {
		return m_entrydata < rhs.m_entrydata ;
	}

	bool VariantEntry::operator<=( VariantEntry const& rhs ) const {
		return m_entrydata < rhs.m_entrydata || m_entrydata == rhs.m_entrydata ;
	}

	bool VariantEntry::operator>( VariantEntry const& rhs ) const {
		return !( m_entrydata < rhs.m_entrydata ) && !( m_entrydata == rhs.m_entrydata ) ;
	}

	bool VariantEntry::operator>=( VariantEntry const& rhs ) const {
		return !( m_entrydata < rhs.m_entrydata ) ;
	}
	
	std::ostream& operator<<( std::ostream& ostr, VariantEntry const& entry ) {
		return ostr << entry.m_entrydata ;
	}
}
