#ifndef GENFILE_VARIANTENTRY_HPP
#define GENFILE_VARIANTENTRY_HPP

#include <string>
#include <boost/variant.hpp>
#include "genfile/MissingValue.hpp"

namespace genfile {
	class VariantEntry {
	public:
		template< typename T > VariantEntry( T const& value ) ;
		template< typename T > VariantEntry& operator=( T const& value ) ;
		VariantEntry() ; // Initialise with MissingValue.
	public:
		bool is_missing() const ;
		bool is_string() const ;
		bool is_int() const ;
		bool is_double() const ;
		template< typename T > T as() const ;
	public:
		bool operator==( VariantEntry const& rhs ) const ;
		bool operator<( VariantEntry const& rhs ) const ;
		friend std::ostream& operator<<( std::ostream&, VariantEntry const& ) ;
	private:
		typedef boost::variant< MissingValue, std::string, int, double > EntryData ;
		EntryData m_entrydata ;
	private:
		// prevent implicit conversion to VariantEntry in comparisons.
		template< typename T > bool operator==( T const& rhs ) const ; 
		template< typename T > bool operator<( T const& rhs ) const ; // this prevents implicit conversion to VariantEntry.
	} ;
	
	template< typename T >
	VariantEntry::VariantEntry( T const& value ):
		m_entrydata( value )
	{}

	template< typename T >
	VariantEntry& VariantEntry::operator=( T const& value ) {
		this->m_entrydata = value ;
		return *this ;
	}
	
	template< typename T >
	T VariantEntry::as() const {
		return boost::get< T >( m_entrydata ) ;
	}
	
	// The following specialisation allows the user to get either a double or integer entry as a double.
	// This simplifies some uses.
	template<> double VariantEntry::as() const ;
}

#endif
