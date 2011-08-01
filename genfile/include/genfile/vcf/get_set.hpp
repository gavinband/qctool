#ifndef GENFILE_VCF_GET_SET_HPP
#define GENFILE_VCF_GET_SET_HPP

#include "genfile/VariantEntry.hpp"
#include "genfile/Error.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "genfile/VariantDataReader.hpp"

namespace genfile {
	namespace vcf {
		template< typename Setter >
		struct GenotypeProbabilitySetter: public VariantDataReader::PerSampleSetter
		{
			GenotypeProbabilitySetter( Setter const& setter ): m_setter( setter ) {}
			GenotypeProbabilitySetter( SingleSNPGenotypeProbabilities& probabilities ):
				m_setter( genfile::set_genotypes( probabilities ) )
			{}

			void set_number_of_samples( std::size_t n ) {
				// destination is supposed to know its size.
				m_number_of_samples = n ;
			}

			void set_sample( std::size_t n ) {
				assert( n < m_number_of_samples ) ;
				m_sample = n ;
				m_missing = false ;
			}
			void set_number_of_entries( std::size_t n ) {
				m_number_of_entries = n ;
				m_entry_i = 0 ;
			}

			void operator()( MissingValue const value ) {
				// if any prob is missing, all are.
				m_missing = true ;
				if( ++m_entry_i == m_number_of_entries ) {
					set() ;
				}
			}

			template< typename T >
			void store( T const value ) {
				if( m_number_of_entries == 2 ) {
					// Treat as two calls.
					assert( value == 0 || value == 1 ) ;
					m_A += ( value == 0 ) ? 1 : 0 ;
					m_B += ( value == 0 ) ? 0 : 1 ;
				}
				else if( m_number_of_entries == 3 || m_number_of_entries == 4 ) {
					// treat as probabilities.  Ignore the fourth probability, which we interpret as NULL call.
					if( m_entry_i < 3 ) {
						m_store[ m_entry_i ] = value ;
					}
				}
				else {
					assert(0) ;
				}
				++m_entry_i ;
				if( m_entry_i == m_number_of_entries ) {
					set() ;
				}
			}

			void set() {
				if( m_missing ) {
					m_setter( m_sample, 0.0, 0.0, 0.0 ) ;
				}
				else if( m_number_of_entries == 2 ) {
					if( m_A == 0 && m_B == 2 ) {
						m_setter( m_sample, 0.0, 0.0, 1.0 ) ;
					}
					else if( m_A == 1 && m_B == 1 ) {
						m_setter( m_sample, 1.0, 0.0, 0.0 ) ;
					}
					else {
						m_setter( m_sample, 0.0, 1.0, 0.0 ) ;
					}
				}
				else {
					m_setter( m_sample, m_store[0], m_store[1], m_store[2] ) ;
				}
			}

			void operator()( Integer const value ) {
				store( value ) ;
			}

			void operator()( double const value ) {
				store( value ) ;
			}

		private:
			Setter const& m_setter ;
			std::size_t m_number_of_samples ;
			std::size_t m_sample ;
			std::size_t m_number_of_entries ;
			std::size_t m_entry_i ;
			double m_store[3] ;
			double m_A, m_B ;
			bool m_missing ;
		} ;
		
		template< typename Setter >
		GenotypeProbabilitySetter< Setter > make_genotype_probability_setter( Setter setter ) {
			return GenotypeProbabilitySetter< Setter >( setter ) ;
		}
	}
}

#endif
