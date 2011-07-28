#ifndef GENFILE_VARIANT_DATA_READER_HPP
#define GENFILE_VARIANT_DATA_READER_HPP

#include <memory>
#include <vector>
#include <string>
#include <set>
#include "boost/noncopyable.hpp"
#include "boost/function.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/get_set.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"

namespace genfile {
	class VariantDataReader: public boost::noncopyable
	{
	public:
		typedef std::auto_ptr< VariantDataReader > UniquePtr ;
		typedef genfile::VariantEntry Entry ;
	public:
		typedef boost::function< void ( std::size_t i, std::vector< Entry >& ) > PerSampleSetter ;
		typedef boost::function< void ( std::string ) > SpecSetter ;
	public:
		virtual ~VariantDataReader() {} ;
		virtual VariantDataReader& get( std::string const& spec, PerSampleSetter setter ) = 0 ;
		virtual bool supports( std::string const& spec ) const = 0 ;
		virtual void get_supported_specs( SpecSetter ) const = 0 ;
		
		struct VectorSetter {
		public:
			VectorSetter( std::vector< std::vector< Entry > >& values ):
				m_values( values )
			{}

			VectorSetter( VectorSetter const& other ):
				m_values( other.m_values )
			{}
			
			void operator()( std::size_t i, std::vector< VariantDataReader::Entry >& values ) {
				if( m_values.size() < (i+1) ) {
					m_values.resize( i+1 ) ;
				}
				m_values[i].swap( values ) ;
			}
		private:
			std::vector< std::vector< Entry > >& m_values ;
		} ;
		
		static VectorSetter set( std::vector< std::vector< Entry > >& values ) {
			return VectorSetter( values ) ;
		}

		template< typename Setter >
		struct GenotypeSetter {
		public:
			typedef VariantDataReader::PerSampleSetter PerSampleSetter ;

			GenotypeSetter( Setter setter ): m_setter( setter ) {}

			void operator()( std::size_t i, std::vector< VariantDataReader::Entry > const& values ) {
				double AA = 0.0, AB = 0.0, BB = 0.0 ; // zero genotypes to represent missing call.
				if( !values.empty() ) {
					if(
						values[0].is_int()
					) {
						if(
							values.size() == 2
						) {
							int A = values[0].as< int >(),
								B = values[1].as< int >() ;

							if( A >= 0 && A < 2 && B >= 0 && B < 2 ) {
								AA = ( A == 0 && B == 0 ) ;
								AB = ( A + B == 1 ) ;
								BB = ( A == 1 && B == 1 ) ;
							}
						}
					}
					else if( values[0].is_double() && ( values.size() == 3 || values.size() == 4 )) {
						AA = values[0].is_missing() ? 0.0 : values[0].as< double >() ;
						AB = values[1].is_missing() ? 0.0 : values[1].as< double >() ;
						BB = values[2].is_missing() ? 0.0 : values[2].as< double >() ;
					}
				}
				m_setter( i, AA, AB, BB ) ;
			}
			
		private:
			Setter m_setter ;
		} ;

		template< typename Setter >
		static GenotypeSetter< Setter > set( Setter setter ) {
			return GenotypeSetter< Setter >( setter ) ;
		}

		static GenotypeSetter< ::genfile::GenotypeSetter< SingleSNPGenotypeProbabilities > >
			set( SingleSNPGenotypeProbabilities& genotypes ) {
			return set( set_genotypes( genotypes ) ) ;
		}
	} ;
}

#endif
