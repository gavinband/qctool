
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_VCF_GET_SET_HPP
#define GENFILE_VCF_GET_SET_HPP

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
			std::size_t m_sample ;
			std::size_t m_entry_i ;
		} ;
		
		struct GenotypeSetterBase: public VariantDataReader::PerSampleSetter
		{
			GenotypeSetterBase( std::string const& scale = "identity" ) ;
			virtual ~GenotypeSetterBase() throw() ;

			virtual void set_number_of_samples( std::size_t n ) ;
			virtual void set_sample( std::size_t n ) ;
			virtual void set_number_of_entries( std::size_t n ) ;
			virtual void operator()( MissingValue const value ) ;
			virtual void operator()( Integer const value ) ;
			virtual void operator()( double const value ) ;

		protected:
			virtual void set( std::size_t, double, double, double ) = 0 ;

		private:
			enum ProbabilityScale { eIdentityScale = 0, ePhredScale = 1 } ;
			ProbabilityScale m_probability_scale ;
			std::size_t m_number_of_samples ;
			std::size_t m_sample ;
			std::size_t m_number_of_entries ;
			std::size_t m_entry_i ;
			double m_store[3] ;
			double m_A, m_B ;
			bool m_missing ;

			void set() ;

			template< typename T >
			void store( T const value ) {
				if( m_number_of_entries == 1 ) {
					// Treat as a dosage
					assert( value == 0 || value == 1 || value == 2 ) ;
					m_A = 2 - value ;
					m_B = value ;
				} else if( m_number_of_entries == 2 ) {
					// Treat as two calls.
					assert( value == 0 || value == 1 ) ;
					m_A += ( value == 0 ) ? 1 : 0 ;
					m_B += ( value == 0 ) ? 0 : 1 ;
				}
				else if( m_number_of_entries == 3 || m_number_of_entries == 4 ) {
					// treat as probabilities.  Ignore the fourth probability, which we interpret as NULL call.
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
/*
		template< typename Setter >
		struct GenotypeSetter: public GenotypeSetterBase
		{
			GenotypeSetter( Setter const& setter ): m_setter( setter ) {}
			void set( std::size_t sample_i, double AA, double AB, double BB ) {
				m_setter( sample_i, AA, AB, BB ) ;
			}
		private:
			Setter const& m_setter ;
		} ;
*/
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
			void set_number_of_samples( std::size_t n ) ;
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
			void set_number_of_samples( std::size_t n ) ;
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
			void set_number_of_samples( std::size_t n ) ;
			void set( std::size_t sample_i, double AA, double AB, double BB ) ;
		private:
			std::vector< VariantEntry >& m_result ;
			double const m_threshhold ;
		} ;

		// Genotype setter which stores hard genotype calls as integers with specified values for each genotype.
		template<>
		struct ThreshholdingGenotypeSetter< std::vector< int > >: public GenotypeSetterBase
		{
			ThreshholdingGenotypeSetter( std::vector< int >& result, double threshhold, int missing_value = -1, int AA_value = 0, int AB_value = 1, int BB_value = 2 ) ;
			void set_number_of_samples( std::size_t n ) ;
			void set( std::size_t sample_i, double AA, double AB, double BB ) ;
		private:
			std::vector< int >& m_result ;
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

			void set_number_of_samples( std::size_t n ) { m_number_of_samples = n ; }
			void set_sample( std::size_t n ) {
				assert( n < m_number_of_samples ) ; 
				m_sample = n ;
				m_entry_i = 0 ;
			}
			void set_number_of_entries( std::size_t n ) {
				if( m_sample == 0 ) {
					m_result.resize( m_number_of_samples, n ) ;
					if( m_nonmissingness ) {
						m_nonmissingness->resize( m_number_of_samples, n ) ;
						m_nonmissingness->setZero() ;
					}
				}
				else if( !n == m_result.rows() ) {
					throw BadArgumentError( "genfile::vcf::MatrixSetter::set_number_of_entries()", "n" ) ;
				}
			}
			void operator()( MissingValue const value ) {
				m_result( m_sample, m_entry_i++ ) = m_missing_value ;
			}
			void operator()( double const value ) {
				if( m_nonmissingness ) {
					(*m_nonmissingness )( m_sample, m_entry_i ) = 1 ;
				}
				m_result( m_sample, m_entry_i++ ) = value ;
			}

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
