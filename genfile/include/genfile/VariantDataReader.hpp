#ifndef GENFILE_VARIANT_DATA_READER_HPP
#define GENFILE_VARIANT_DATA_READER_HPP

#include <memory>
#include <vector>
#include <string>
#include <set>
#include <boost/noncopyable.hpp>
#include <boost/function.hpp>
#include "genfile/VariantEntry.hpp"
#include "genfile/get_set.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "vcf/Types.hpp"

namespace genfile {
	class VariantDataReader: public boost::noncopyable
	{
	public:
		typedef std::auto_ptr< VariantDataReader > UniquePtr ;
		typedef genfile::VariantEntry Entry ;
	public:
		struct PerSampleSetter: public vcf::EntriesSetter, public boost::noncopyable {
			virtual ~PerSampleSetter() throw() {}
			virtual void set_number_of_samples( std::size_t n ) = 0 ;
			virtual void set_sample( std::size_t i ) = 0 ;
		} ;
		typedef boost::function< void ( std::string ) > SpecSetter ;
	public:
		virtual ~VariantDataReader() {} ;
		virtual VariantDataReader& get( std::string const& spec, PerSampleSetter& setter ) = 0 ;

		virtual bool supports( std::string const& spec ) const = 0 ;
		virtual void get_supported_specs( SpecSetter ) const = 0 ;
		
		struct VectorSetter: public PerSampleSetter {
		public:
			VectorSetter( std::vector< std::vector< Entry > >& data ):
				m_data( data )
			{}

			void set_number_of_samples( std::size_t n ) { m_data.resize( n ) ; }
			void set_sample( std::size_t n ) { assert( n < m_data.size() ) ; m_sample = n ; }
			void set_number_of_entries( std::size_t n ) { m_data[ m_sample ].resize( n ) ; m_entry_i = 0 ; }

		private:
			template< typename T >
			void set( T value ) {
				assert( m_entry_i < m_data[ m_sample ].size() ) ;
				m_data[ m_sample ][ m_entry_i++ ] = value ;
			}
		public:
			void operator()( MissingValue const value ) { set( value ) ; }
			void operator()( std::string& value ) { set( value ) ; }
			void operator()( Integer const value ) { set( value ) ; }
			void operator()( double const value ) { set( value ) ; }

		private:
			std::vector< std::vector< Entry > >& m_data ;
			std::size_t m_number_of_samples ;
			std::size_t m_sample ;
			std::size_t m_entry_i ;
		} ;
		
		// Convenience method.
		VariantDataReader& get( std::string const& spec, std::vector< std::vector< Entry > >& data ) {
			VectorSetter setter( data ) ;
			get( spec, setter ) ;
			return *this ;
		}

		// Convenience method setting SingleSNPGenotypeProbabilities.
		VariantDataReader& get( std::string const& spec, SingleSNPGenotypeProbabilities& data ) ;

		template< typename Setter >
		static GenotypeSetter< Setter > set( Setter setter ) {
			return GenotypeSetter< Setter >( setter ) ;
		}
	} ;
}

#endif
