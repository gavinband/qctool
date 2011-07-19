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
	
	bool VariantEntry::operator<( VariantEntry const& rhs ) const {
		return m_entrydata < rhs.m_entrydata ;
	}
	
	std::ostream& operator<<( std::ostream& ostr, VariantEntry const& entry ) {
		return ostr << entry.m_entrydata ;
	}
	
	namespace impl {
		class serialization_sizer: public boost::static_visitor< std::size_t >
		{
			typedef VariantEntry::Integer Integer ;
		public:
			serialization_sizer() {}
			serialization_sizer( serialization_sizer const& other ) {}

			std::size_t operator()( MissingValue const& operand ) const ;
			std::size_t operator()( Integer const& operand ) const ;
			std::size_t operator()( std::string const& operand ) const ;
			std::size_t operator()( double const& operand ) const ;
			std::size_t operator()( Chromosome const& operand ) const ;
			std::size_t operator()( GenomePosition const& operand ) const ;
		} ;

		class serializer: public boost::static_visitor< char* >
		{
			typedef VariantEntry::Integer Integer ;
		public:
			serializer( char* begin, char* const end ): m_begin( begin ), m_end( end ) {}
			serializer( serializer const& other ): m_begin( other.m_begin ), m_end( other.m_end ) {}

			char* operator()( MissingValue const& operand ) ;
			char* operator()( Integer const& operand ) ;
			char* operator()( std::string const& operand ) ;
			char* operator()( double const& operand ) ;
			char* operator()( Chromosome const& operand ) ;
			char* operator()( GenomePosition const& operand ) ;

		private:
			char* m_begin ;
			char* const m_end ;
		} ;

		class deserializer: public boost::static_visitor< char const* >
		{
			typedef VariantEntry::Integer Integer ;
		public:
			deserializer( char const* begin, char const* const end ): m_begin( begin ), m_end( end ) {}
			deserializer( deserializer const& other ): m_begin( other.m_begin ), m_end( other.m_end ) {}
			
			char const* operator()( MissingValue& operand ) ;
			char const* operator()( Integer& operand ) ;
			char const* operator()( std::string& operand ) ;
			char const* operator()( double& operand ) ;
			char const* operator()( Chromosome& operand ) ;
			char const* operator()( GenomePosition& operand ) ;

		private:
			char const* m_begin ;
			char const* const m_end ;
		};

		std::size_t serialization_sizer::operator()( MissingValue const& operand ) const {
			return 0 ;
		}

		char* serializer::operator()( MissingValue const& operand )
		{
			assert( m_end >= m_begin ) ;
			return m_begin ;
		}

		char const* deserializer::operator()( MissingValue& operand )
		{
			assert( m_end >= m_begin ) ;
			return m_begin ;
		}

		std::size_t serialization_sizer::operator()( Integer const& operand ) const {
			return sizeof( Integer ) ;
		}

		char* serializer::operator()( Integer const& operand )
		{
			m_begin = genfile::write_little_endian_integer( m_begin, m_end, operand ) ;
			return m_begin ;
		}

		char const* deserializer::operator()( Integer& operand )
		{
			m_begin = read_little_endian_integer( m_begin, m_end, &operand ) ;
			return m_begin ;
		}
		
		std::size_t serialization_sizer::operator()( std::string const& operand ) const {
			return get_small_integer_size( operand.size() ) + operand.size() ;
		}

		char* serializer::operator()( std::string const& operand )
		{
			m_begin = write_small_integer( m_begin, m_end, operand.size() ) ;
			assert( m_end >= m_begin + operand.size() ) ;
			m_begin = std::copy( operand.begin(), operand.end(), m_begin ) ;
			return m_begin ;
		}

		char const* deserializer::operator()( std::string& operand )
		{
			std::size_t size ;
			m_begin = read_small_integer( m_begin, m_end, &size ) ;
			assert( m_end >= m_begin + size ) ;
			operand.assign( m_begin, m_begin + size ) ;
			m_begin += size ;
			return m_begin ;
		}
		
		std::size_t serialization_sizer::operator()( double const& operand ) const {
			return sizeof( double ) ;
		}

		char* serializer::operator()( double const& operand )
		{
			assert( m_end >= m_begin + sizeof( double ) ) ;
			char const* ptr = reinterpret_cast< char const* >( &operand ) ;
			m_begin = std::copy( ptr, ptr + sizeof( double ), m_begin ) ;
			return m_begin ;
		}

		char const* deserializer::operator()( double& operand )
		{
			assert( m_end >= m_begin + sizeof( double ) ) ;
			char* ptr = reinterpret_cast< char* >( &operand ) ;
			m_begin = std::copy( m_begin, m_begin + sizeof( double ), ptr ) ;
			return m_begin ;
		}

		std::size_t serialization_sizer::operator()( Chromosome const& operand ) const {
			return 1 ;
		}

		char* serializer::operator()( Chromosome const& operand )
		{
			assert( m_end >= m_begin + 1 ) ;
			*m_begin++ = char( ChromosomeEnum( operand ) ) ;
			return m_begin ;
		}

		char const* deserializer::operator()( Chromosome& operand )
		{
			assert( m_end >= m_begin + 1 ) ;
			operand = Chromosome( *m_begin++ ) ;
			return m_begin ;
		}

		std::size_t serialization_sizer::operator()( GenomePosition const& operand ) const {
			return 1 + sizeof( Position ) ;
		}

		char* serializer::operator()( GenomePosition const& operand )
		{
			assert( m_end >= m_begin + 1 + sizeof( Position ) ) ;
			*m_begin++ = char( ChromosomeEnum( operand.chromosome() ) ) ;
			m_begin = write_little_endian_integer( m_begin, m_end, operand.position() ) ;
			return m_begin ;
		}

		char const* deserializer::operator()( GenomePosition& operand )
		{
			assert( m_end >= m_begin + 1 + sizeof( Position ) ) ;
			operand.chromosome() = Chromosome( *m_begin++ ) ;
			m_begin = read_little_endian_integer( m_begin, m_end, &operand.position() ) ;
			return m_begin ;
		}
	}

	std::size_t VariantEntry::get_serialized_size() const {
		return 1 + boost::apply_visitor( impl::serialization_sizer(), m_entrydata ) ;
	}

	char* VariantEntry::serialize( char* buffer, char* const end ) const {
		assert( end >= buffer + 1 ) ;
		*buffer++ = char( m_entrydata.which() ) ;
		impl::serializer serializer( buffer, end ) ;
		return boost::apply_visitor( serializer, m_entrydata ) ;
	}

	const char* VariantEntry::deserialize( char const* buffer, const char* const end ) {
		assert( end >= buffer + 1 ) ;
		char type = *buffer++ ;
		switch( type ) {
			case eMissing: m_entrydata = MissingValue() ; break ;
			case eInteger: m_entrydata = Integer( 0 ) ; break ;
			case eString: m_entrydata = std::string() ; break ;
			case eDouble: m_entrydata = double( 0.0 ) ; break ;
			case eChromosome: m_entrydata = Chromosome() ; break ;
			case eGenomePosition: m_entrydata = GenomePosition() ; break ;
		}
		impl::deserializer deserializer( buffer, end ) ;
		return boost::apply_visitor( deserializer, m_entrydata ) ;
	}
}
