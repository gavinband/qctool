#ifndef GENFILE_VARIANTENTRY_HPP
#define GENFILE_VARIANTENTRY_HPP

#include <string>
#include <boost/variant.hpp>
#include "genfile/MissingValue.hpp"
#include "genfile/GenomePosition.hpp"

namespace genfile {
	namespace impl {
		class serializer ;
		class deserializer ;
	}
	class VariantEntry {
	public:
		template< typename T > VariantEntry( T const& value ) ;
		VariantEntry( int const value ) ;
		template< typename T > VariantEntry& operator=( T const& value ) ;
		VariantEntry& operator=( int const value ) ;
		VariantEntry() ; // Initialise with MissingValue.
	public:
		bool is_missing() const ;
		bool is_string() const ;
		bool is_int() const ;
		bool is_double() const ;
		bool is_chromosome() const ;
		bool is_position() const ;
		template< typename T > T as() const ;

		std::size_t get_serialized_size() const ;
		char* serialize( char* buffer, char* const end ) const ;
		char const* deserialize( char const* buffer, char const* const end ) ;
	public:
		bool operator==( VariantEntry const& rhs ) const ;
		bool operator<( VariantEntry const& rhs ) const ;
		friend std::ostream& operator<<( std::ostream&, VariantEntry const& ) ;
		typedef int64_t Integer ;
	private:
		enum Types { eMissing = 0, eString = 1, eInteger = 2, eDouble = 3, eChromosome = 4, eGenomePosition = 5 } ;
		typedef boost::variant< MissingValue, std::string, Integer, double, Chromosome, GenomePosition > EntryData ;
		EntryData m_entrydata ;
		friend class impl::serializer ;
		friend class impl::deserializer ;
		
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
	// The following specialisation allows the user to get either a value as an int.
	template<> int VariantEntry::as() const ;
}

#endif
