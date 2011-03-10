#include <boost/variant.hpp>
#include "genfile/VariantEntry.hpp"

namespace genfile {
	VariantEntry::VariantEntry():
		m_entrydata( MissingValue() )
	{}
	
	bool VariantEntry::is_missing() const {
		return boost::get< MissingValue >( &m_entrydata ) ;
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
