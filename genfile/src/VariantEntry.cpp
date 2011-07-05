#include <boost/variant.hpp>
#include "genfile/VariantEntry.hpp"

namespace genfile {
	VariantEntry::VariantEntry():
		m_entrydata( MissingValue() )
	{}
	
	bool VariantEntry::is_missing() const {
		return m_entrydata.which() == 0 ;
	}

	bool VariantEntry::is_string() const {
		return m_entrydata.which() == 1 ;
	}

	bool VariantEntry::is_int() const {
		return m_entrydata.which() == 2 ;
	}

	bool VariantEntry::is_double() const {
		return m_entrydata.which() == 3 ;
	}

	bool VariantEntry::is_chromosome() const {
		return m_entrydata.which() == 4 ;
	}

	bool VariantEntry::is_position() const {
		return m_entrydata.which() == 5 ;
	}
	
	template<> double VariantEntry::as() const {
		if( int const* v = boost::get< int >( &m_entrydata )) {
			return double( *v ) ;
		}
		else {
			return boost::get< double >( m_entrydata ) ;
		} 
	}
	
	bool VariantEntry::operator==( VariantEntry const& rhs ) const {
		return m_entrydata == rhs.m_entrydata ;
	}
	
	bool VariantEntry::operator<( VariantEntry const& rhs ) const {
		return m_entrydata < rhs.m_entrydata ;
	}
	
	std::ostream& operator<<( std::ostream& ostr, VariantEntry const& entry ) {
		return ostr << entry.m_entrydata ;
	}
}
