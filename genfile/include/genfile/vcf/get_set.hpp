
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_VCF_GET_SET_HPP
#define GENFILE_VCF_GET_SET_HPP

#include <cassert>
#include <limits>
#include "genfile/VariantEntry.hpp"
#include "genfile/Error.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "genfile/VariantDataReader.hpp"

namespace genfile {
	namespace vcf {
		struct VectorSetter: public VariantDataReader::PerSampleSetter {
		public:
			VectorSetter( std::vector< std::vector< Entry > >& data ):
				m_data( data )
			{}

			void initialise( std::size_t nSamples, std::size_t nAlleles ) { assert( nAlleles == 2 ); m_data.resize( nSamples ) ; }
			bool set_sample( std::size_t n ) { assert( n < m_data.size() ) ; m_sample = n ; return true ; }
			void set_number_of_entries( uint32_t, std::size_t n, OrderType const order_type, ValueType const value_type ) {
				m_data[ m_sample ].resize( n ) ; m_entry_i = 0 ;
			}

		private:
			template< typename T >
			void set( T value ) {
				assert( m_entry_i < m_data[ m_sample ].size() ) ;
				m_data[ m_sample ][ m_entry_i++ ] = value ;
			}
		public:
			void set_value( MissingValue const value ) { set( value ) ; }
			void set_value( std::string& value ) { set( value ) ; }
			void set_value( Integer const value ) { set( value ) ; }
			void set_value( double const value ) { set( value ) ; }

			void finalise() {}
		private:
			std::vector< std::vector< Entry > >& m_data ;
			std::size_t m_sample ;
			std::size_t m_entry_i ;
		} ;
		
		struct GenotypeSetterBase: public VariantDataReader::PerSampleSetter
		{
			GenotypeSetterBase( std::string const& scale = "identity" ) ;
			~GenotypeSetterBase() throw() ;

			void initialise( std::size_t nSamples, std::size_t nAlleles ) ;
			bool set_sample( std::size_t n ) ;
			void set_number_of_entries( uint32_t, std::size_t n, OrderType const, ValueType const ) ;
			void set_value( MissingValue const value ) ;
			void set_value( Integer const value ) ;
			void set_value( double const value ) ;
			void finalise() {}

		protected:
			virtual void set( std::size_t, double, double, double ) = 0 ;

		private:
			enum ProbabilityScale { eIdentityScale = 0, ePhredScale = 1 } ;
			ProbabilityScale m_probability_scale ;
			std::size_t m_number_of_samples ;
			std::size_t m_sample ;
			OrderType m_order_type ;
			ValueType m_value_type ;
			std::size_t m_number_of_entries ;
			std::size_t m_entry_i ;
			double m_store[3] ;
			double m_A, m_B ;
			bool m_missing ;

			void set() ;

			template< typename T >
			void store( T const value ) {
				if( m_value_type == eDosage && m_order_type == ePerSample ) {
					// A single dosage value
					assert( value == 0 || value == 1 || value == 2 ) ;
					m_A = 2 - value ;
					m_B = value ;
				}
				else if( m_value_type == eAlleleIndex && ( m_order_type == ePerOrderedHaplotype || m_order_type == ePerUnorderedHaplotype )) {
					// GT-style genotype calls
					assert( value == 0 || value == 1 ) ;
					m_A += ( value == 0 ) ? 1 : 0 ;
					m_B += ( value == 0 ) ? 0 : 1 ;
				}
				else if( m_value_type == eProbability && m_order_type == ePerUnorderedGenotype ) {
					// genotype probabilities.
					if( m_entry_i < 3 ) {
						switch( m_probability_scale ) {
							case ePhredScale:
								m_store[ m_entry_i ] = std::pow( -value / 10.0, 10 ) ;
								break ;
							case eIdentityScale:
								m_store[ m_entry_i ] = value ;
								break ;
						}
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
		} ;

		template< typename Setter > struct GenotypeSetter ;

		template<>
		struct GenotypeSetter< boost::function< void ( std::size_t, double, double, double ) > >: public GenotypeSetterBase
		{
			typedef boost::function< void ( std::size_t, double, double, double ) > Setter ;
			GenotypeSetter( Setter const& setter ): GenotypeSetterBase( "identity" ), m_setter( setter ) {}
			void set( std::size_t sample_i, double AA, double AB, double BB ) {
				m_setter( sample_i, AA, AB, BB ) ;
			}
		private:
			Setter const& m_setter ;
		} ;

		// Genotype setter which stores genotype probabilities in a SingleSNPGenotypeProbabilities object.
		template<>
		struct GenotypeSetter< SingleSNPGenotypeProbabilities >: public GenotypeSetterBase
		{
			GenotypeSetter( SingleSNPGenotypeProbabilities& result ) ;
			void initialise( std::size_t nSamples, std::size_t nAlleles ) ;
			void set( std::size_t sample_i, double AA, double AB, double BB ) ;
		private:
			SingleSNPGenotypeProbabilities& m_result ;
		} ;

		// Genotype setter which stores genotype probabilities as a single vector of doubles.
		// Genotype call probabilities for individual i are at indices (3*i), (3*i)+1, (3*i)+2.
		// A missing call is indicated with three zeroes.
		template<>
		struct GenotypeSetter< std::vector< double > >: public GenotypeSetterBase
		{
			GenotypeSetter( std::vector< double >& result ) ;
			void initialise( std::size_t nSamples, std::size_t nAlleles ) ;
			void set( std::size_t sample_i, double AA, double AB, double BB ) ;
		private:
			std::vector< double >& m_result ;
		} ;

		template< typename T > struct ThreshholdingGenotypeSetter {} ;

		// Genotype setter which stores hard genotype calls as VariantEntries.
		// They are either integers 0, 1, or 2, or MissingValue.
		template<>
		struct ThreshholdingGenotypeSetter< std::vector< VariantEntry > >: public GenotypeSetterBase
		{
			ThreshholdingGenotypeSetter( std::vector< VariantEntry >& result, double threshhold ) ;
			void initialise( std::size_t nSamples, std::size_t nAlleles ) ;
			void set( std::size_t sample_i, double AA, double AB, double BB ) ;
		private:
			std::vector< VariantEntry >& m_result ;
			double const m_threshhold ;
		} ;

		// Genotype setter which stores hard genotype calls as integers with specified values for each genotype.
		template< typename T >
		struct ThreshholdingGenotypeSetter< std::vector< T > >: public GenotypeSetterBase
		{
			ThreshholdingGenotypeSetter( std::vector< T >& result, double threshhold, T missing_value = -1, T AA_value = 0, T AB_value = 1, T BB_value = 2 ):
				GenotypeSetterBase( "identity" ),
				m_result( result ),
				m_missing_value( missing_value ),
				m_AA_value( AA_value ),
				m_AB_value( AB_value ),
				m_BB_value( BB_value ),
				m_threshhold( threshhold )
			{}
				
			void initialise( std::size_t nSamples, std::size_t nAlleles ) {
				assert( nAlleles == 2 ) ;
				m_result.clear() ;
				m_result.resize( nSamples, -1 ) ;
				GenotypeSetterBase::initialise( nSamples, nAlleles ) ;
			}
			
			void set( std::size_t sample_i, double AA, double AB, double BB ) {
				if( AA > m_threshhold ) {
					m_result[ sample_i ] = m_AA_value ;
				}
				else if( AB > m_threshhold ) {
					m_result[ sample_i ] = m_AB_value ;
				}
				else if( BB > m_threshhold ) {
					m_result[ sample_i ] = m_BB_value ;
				}
				else {
					m_result[ sample_i ] = m_missing_value ;
				}
			}
		private:
			std::vector< T >& m_result ;
			int const m_missing_value ;
			int const m_AA_value ;
			int const m_AB_value ;
			int const m_BB_value ;
			double const m_threshhold ;
		} ;

		template< typename Matrix >
		struct MatrixSetter: public VariantDataReader::PerSampleSetter
		{
			MatrixSetter(
				Matrix& result,
				double const missing_value = std::numeric_limits< double >::quiet_NaN()
			):
				m_result( result ),
				m_nonmissingness( 0 ),
				m_missing_value( missing_value )
			{}
			MatrixSetter(
				Matrix& result,
				Matrix& nonmissingness,
				double const missing_value = 0
			):
				m_result( result ),
				m_nonmissingness( &nonmissingness ),
				m_missing_value( missing_value )
			{}

			void initialise( std::size_t nSamples, std::size_t nAlleles ) { assert( nAlleles == 2 ) ; m_number_of_samples = nSamples ; }
			bool set_sample( std::size_t n ) {
				assert( n < m_number_of_samples ) ; 
				m_sample = n ;
				m_entry_i = 0 ;
				return true ;
			}
			void set_number_of_entries( uint32_t, std::size_t n, OrderType const, ValueType const ) {
				if( m_sample == 0 ) {
					m_result.resize( m_number_of_samples, n ) ;
					if( m_nonmissingness ) {
						m_nonmissingness->resize( m_number_of_samples, n ) ;
						m_nonmissingness->setZero() ;
					}
				}
				else if( n != m_result.cols() ) {
					throw BadArgumentError( "genfile::vcf::MatrixSetter::set_number_of_entries()", "n" ) ;
				}
			}
			void set_value( MissingValue const value ) {
				m_result( m_sample, m_entry_i++ ) = m_missing_value ;
				if( m_nonmissingness ) {
					(*m_nonmissingness )( m_sample, m_entry_i ) = 0 ;
				}
			}
			void set_value( double const value ) {
				if( m_nonmissingness ) {
					(*m_nonmissingness )( m_sample, m_entry_i ) = 1 ;
				}
				m_result( m_sample, m_entry_i++ ) = value ;
			}

			void finalise() {} ;
			
		private:
			Matrix& m_result ;
			Matrix* m_nonmissingness ;
			double const m_missing_value ;
			std::size_t m_number_of_samples ;
			std::size_t m_sample ;
			std::size_t m_number_of_entries ;
			std::size_t m_entry_i ;
		} ;
	}
}

#endif
