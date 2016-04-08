
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <utility>
#include <map>
#include <iostream>
#include <boost/noncopyable.hpp>
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/shared_ptr.hpp>
#include "genfile/MissingValue.hpp"
#include "genfile/types.hpp"
#include "genfile/string_utils/string_utils.hpp"
#include "genfile/bgen/bgen.hpp"
#include "genfile/Error.hpp"
#include "genfile/VariantDataReader.hpp"

// #define DEBUG_TO_GP 1

namespace genfile {

	template< typename Setter >
	class ToGPImpl: public boost::noncopyable {
	public:
		typedef std::auto_ptr< ToGPImpl > UniquePtr ;
		typedef boost::shared_ptr< ToGPImpl > SharedPtr ;

		virtual ~ToGPImpl() {}
		virtual void initialise( Setter& setter, std::size_t number_of_alleles, uint32_t ploidy, std::size_t number_of_entries ) = 0 ;
		virtual bool initialised() const = 0 ;
		virtual void set_value( std::size_t value_i, int64_t const value ) = 0 ;
		virtual void set_value( std::size_t value_i, genfile::MissingValue const value ) = 0 ;
		virtual void set_value( std::size_t value_i, double const value ) = 0 ;
		virtual void finalise() = 0 ;
		
	} ;

	template< typename Setter >
	struct ToGP: public VariantDataReader::PerSampleSetter {
	public:
		typedef ToGPImpl< Setter > Impl ;
	private:
		typedef std::map< std::pair< OrderType, ValueType >, typename Impl::SharedPtr > ImplMap ;
		//typedef boost::ptr_map< std::pair< OrderType, ValueType >, Impl > ImplMap ;
	public:
		~ToGP() throw() {}
		
		ToGP( Setter& setter ):
			m_setter( setter )
		{}

		ToGP( ToGP const& other ):
			m_setter( other.m_setter ),
			m_impls( other.m_impls )
		{}

		void add_impl( OrderType const order_type, ValueType const value_type, typename Impl::SharedPtr impl ) {
			m_impls[ std::make_pair( order_type, value_type ) ] = impl ;
		}

		void initialise( std::size_t number_of_samples, std::size_t number_of_alleles ) {
			m_number_of_samples = number_of_samples ;
			m_number_of_alleles = number_of_alleles ;
			m_setter.initialise( number_of_samples, number_of_alleles ) ;
#if DEBUG_TO_GP
			std::cerr << "NEW SNP\n" ;
#endif
		}
		
		bool set_sample( std::size_t i ) {
			if( m_current_impl && m_current_impl->initialised() ) {
				m_current_impl->finalise() ;
			}
			m_current_impl.reset() ;
			return m_setter.set_sample( i ) ;
		}
		
		void set_number_of_entries( uint32_t ploidy, std::size_t n, OrderType const order_type, ValueType const value_type ) {
			typename ImplMap::iterator where = m_impls.find( std::make_pair( order_type, value_type )) ;
			if( where == m_impls.end() ) {
				throw genfile::BadArgumentError(
					"genfile::ToGP::set_number_of_entries()",
					"order_type=" + string_utils::to_string( order_type ) + ", value_type=" + string_utils::to_string(value_type),
					"Unsupported order_type/value_type combination."
				) ;
			} else {
				m_current_impl = where->second ;
				m_current_impl->initialise( m_setter, m_number_of_alleles, ploidy, n ) ;
			}
		}
		
		void set_value( std::size_t value_i, genfile::MissingValue const value ) {
			m_current_impl->set_value( value_i, value ) ;
		}

		void set_value( std::size_t value_i, int64_t const value ) {
			m_current_impl->set_value( value_i, value ) ;
		}
		
		void set_value( std::size_t value_i, double const value ) {
			m_current_impl->set_value( value_i, value ) ;
		}
		
		void finalise() {
			if( m_current_impl->initialised() ) {
				m_current_impl->finalise() ;
			}
			m_setter.finalise() ;
		}
		
	private:
		Setter& m_setter ;
		ImplMap m_impls ;
		std::size_t m_number_of_samples ;
		std::size_t m_number_of_alleles ;
		typename Impl::SharedPtr m_current_impl ;
	} ;

	namespace impl {
		// Generate all possible genotypes for a given ploidy and up to K alleles.
		// The order is the order such that genotypes carrying at least one allele k,...K
		// later in the order than those carrying no alleles k,....,K
		// This order is especially useful because (for the given ploidy) the genotypes at 
		// biallelic variants are ordered in the same way as those at a triallelic variant
		// but carrying none of the third allele - thus the orders are nested.
		//
		// Genotypes are stored as uint16_t with 4 bits per allele.  This gives
		// up to four alleles and up to ploidy of 15.
		std::pair< uint32_t, std::vector< uint16_t > > enumerate_unphased_genotypes( std::size_t ploidy ) ;

		struct GTToGPUnphasedBase {
		protected:
			static std::vector< std::pair< uint32_t, std::vector< uint16_t > > > m_tables ;
		} ;
		
		std::string format_call( uint16_t call, uint32_t bitsPerAllele ) ;
	}

	// This class receives unphased (or phased) GT-style genotypes
	// and outputs unphased genotype probabilities.
	template< typename Setter >
	struct GTToGPUnphased: public ToGPImpl< Setter >, impl::GTToGPUnphasedBase
	{
		typedef typename ToGPImpl< Setter >::SharedPtr SharedPtr ;

		static SharedPtr create() {
			if( m_tables.empty() ) {
				for( std::size_t ploidy = 0; ploidy < 16; ++ploidy ) {
					m_tables.push_back( impl::enumerate_unphased_genotypes( ploidy )) ;
				}
			}
			return SharedPtr( new GTToGPUnphased() ) ;
		}
		
		GTToGPUnphased():
			m_setter( 0 ),
			m_missing( false )
		{
		}

		void initialise( Setter& setter, std::size_t number_of_alleles, uint32_t ploidy, std::size_t number_of_entries ) {
			m_setter = &setter ;
			m_number_of_alleles = number_of_alleles ;
			m_ploidy = ploidy ;
			// Compute number of bits required to store values up to ploidy.
			assert( ploidy < 64 ) ;
			m_number_of_entries = number_of_entries ;
			m_missing = false ;
			m_encoded_call = 0 ;
#if DEBUG_TO_GP
			std::cerr << "Initialised, ploidy = " << m_ploidy << ".\n" ;
#endif
		}
		
		bool initialised() const { return (m_setter != 0)  ; }

		void set_value( std::size_t value_i, genfile::MissingValue const value ) {
			m_missing = true ;
		}

		void set_value( std::size_t value_i, double const value ) {
			throw BadArgumentError( "genfile::GTToGPUnphased::set_value()", "value", "Expected an integer or missing value" ) ;
		}

		void set_value( std::size_t value_i, int64_t const value ) {
			m_encoded_call += encode_call( value, m_tables[ m_ploidy ].first ) ;
#if DEBUG_TO_GP
			std::cerr << "Encoded value #" << value_i << "(" << value << ") encoded as:" << impl::format_call( m_encoded_call ) << "\n" ;
#endif
		}
		
		uint16_t encode_call( int64_t const value, uint32_t const bitsPerAllele ) {
			uint16_t const bitMask = uint16_t( 0xFFFF ) >> ( 16 - bitsPerAllele ) ;
			std::size_t const numberOfAlleles = 16 / bitsPerAllele ;
			assert( value <= numberOfAlleles ) ;
			return(( value == 0 ) ? uint16_t(0) : (1 << ((value-1)*bitsPerAllele))) ;
		}
		
		void finalise() {
			std::size_t const count = bgen::impl::n_choose_k( m_ploidy + m_number_of_alleles - 1, m_number_of_alleles - 1 ) ;
			m_setter->set_number_of_entries( m_ploidy, count, ePerUnorderedGenotype, eProbability ) ;
			
			if( m_missing ) {
#if DEBUG_TO_GP
				std::cerr << "count = " << count << ", call is ./././.\n" ;
#endif
				
				for( std::size_t i = 0 ; i < count; ++i ) {
					m_setter->set_value( i, 0.0 ) ;
				}
			} else {
				std::size_t index_of_nonzero_probability = m_tables[ m_ploidy ].second[ m_encoded_call ] ;
#if DEBUG_TO_GP
				std::cerr << "count = " << count << ", call is: " << impl::format_call( m_encoded_call, m_tables[ m_ploidy ].first ) ;
				std::cerr << ", index is: " << index_of_nonzero_probability << "\n" ;
#endif
				std::size_t i = 0 ;
				for( ; i < index_of_nonzero_probability; ++i ) {
					m_setter->set_value( i, 0.0 ) ;
				}
				m_setter->set_value( i++, 1.0 ) ;
				for( ; i < count; ++i ) {
					m_setter->set_value( i, 0.0 ) ;
				}
			}
			m_setter = 0 ;
		}
		
	private:
		Setter* m_setter ;
		std::size_t m_number_of_alleles ;
		std::size_t m_ploidy ;
		std::size_t m_number_of_entries ;
		uint16_t m_encoded_call ;
		bool m_missing ;
	} ;


	// This class receives unphased (or phased) GT-style genotypes
	// and outputs unphased genotype probabilities.
	template< typename Setter >
	struct GTToGPPhased: public ToGPImpl< Setter >, impl::GTToGPUnphasedBase
	{
		typedef typename ToGPImpl< Setter >::SharedPtr SharedPtr ;


		static SharedPtr create() {
			return SharedPtr( new GTToGPPhased() ) ;
		}
		
		GTToGPPhased():
			m_setter( 0 ),
			m_missing( false )
		{
		}

		void initialise( Setter& setter, std::size_t number_of_alleles, uint32_t ploidy, std::size_t number_of_entries ) {
			m_setter = &setter ;
			m_number_of_alleles = number_of_alleles ;
			m_ploidy = ploidy ;
			m_number_of_entries = number_of_entries ;
			assert( m_ploidy == m_number_of_entries ) ;
			m_missing = false ;
			m_values.resize( m_number_of_entries ) ;
//			std::cerr << "Initialised.\n" ;
		}
		
		bool initialised() const { return ( m_setter != 0 ) ; }

		void set_value( std::size_t value_i, genfile::MissingValue const value ) {
			m_missing = true ;
		}

		void set_value( std::size_t value_i, int64_t const value ) {
			m_values[value_i] = value ;
		}

		void set_value( std::size_t value_i, double const value ) {
			throw BadArgumentError( "genfile::GTToGPUnphased::set_value()", "value", "Expected an integer or missing value" ) ;
		}
		
		void finalise() {
			std::size_t const count = m_number_of_alleles * m_ploidy ;
			m_setter->set_number_of_entries( m_ploidy, count, ePerPhasedHaplotypePerAllele, eProbability ) ;
			
			if( m_missing ) {
//				std::cerr << "count = " << count << ", call is .|.\n" ;
				
				for( std::size_t i = 0 ; i < count; ++i ) {
					m_setter->set_value( i, 0.0 ) ;
				}
			} else {
#if DEBUG_TO_GP
				std::cerr << "ploidy = " << m_ploidy << ", count = " << count << "\n" ;
#endif
				for( std::size_t i = 0; i < m_ploidy; ++i ) {
					std::size_t const index_of_nonzero_probability = m_values[ i ] ;
#if DEBUG_TO_GP
					std::cerr << "allele " << i << ", index is: " << index_of_nonzero_probability << "\n" ;
#endif
					std::size_t allele = 0 ;
					for( ; allele < index_of_nonzero_probability; ++allele ) {
						m_setter->set_value( i*m_number_of_alleles+allele, 0.0 ) ;
					}
					m_setter->set_value( i * m_number_of_alleles + allele, 1.0 ) ;
					for( ++allele; allele < m_number_of_alleles; ++allele ) {
						m_setter->set_value( i*m_number_of_alleles+allele, 0.0 ) ;
					}
				}
			}
			m_setter = 0 ;
		}
		
	private:
		Setter* m_setter ;
		std::size_t m_number_of_alleles ;
		std::size_t m_ploidy ;
		std::size_t m_number_of_entries ;
		std::vector< std::size_t > m_values ;
		bool m_missing ;
	} ;
	
	// This class receives unphased (or phased) GT-style genotypes
	// and outputs unphased genotype probabilities.
	template< typename Setter >
	struct GPToGPUnphased: public ToGPImpl< Setter >, impl::GTToGPUnphasedBase
	{
		typedef typename ToGPImpl< Setter >::SharedPtr SharedPtr ;


		static SharedPtr create() {
			return SharedPtr( new GPToGPUnphased() ) ;
		}
		
		GPToGPUnphased():
			m_setter( 0 )
		{
		}

		void initialise( Setter& setter, std::size_t number_of_alleles, uint32_t ploidy, std::size_t number_of_entries ) {
			m_setter = &setter ;
			m_setter->set_number_of_entries( ploidy, number_of_entries, ePerUnorderedGenotype, eProbability ) ;
		}
		
		bool initialised() const { return ( m_setter != 0 ) ; }

		void set_value( std::size_t value_i, genfile::MissingValue const value ) {
			m_setter->set_value( value_i, value ) ;
		}

		void set_value( std::size_t value_i, int64_t const value ) {
			throw BadArgumentError( "genfile::GTToGPUnphased::set_value()", "value", "Expected a double or missing value" ) ;
		}

		void set_value( std::size_t value_i, double const value ) {
			m_setter->set_value( value_i, value ) ;
		}
		
		void finalise() {
			m_setter = 0 ; 
		}

	private:
		Setter* m_setter ;
	} ;
	
	template< typename Setter >
	ToGP< Setter > to_GP_unphased( Setter& setter ) {
		ToGP< Setter > result( setter ) ;
		result.add_impl( ePerUnorderedHaplotype, eAlleleIndex, GTToGPUnphased< Setter >::create() ) ;
		result.add_impl( ePerOrderedHaplotype, eAlleleIndex, GTToGPUnphased< Setter >::create() ) ;
		result.add_impl( ePerUnorderedGenotype, eProbability, GPToGPUnphased< Setter >::create() ) ;
		return result ;
	}

	template< typename Setter >
	ToGP< Setter > to_GP_phased( Setter& setter ) {
		ToGP< Setter > result( setter ) ;
		result.add_impl( ePerOrderedHaplotype, eAlleleIndex, GTToGPPhased< Setter >::create() ) ;
		return result ;
	}
	
	template< typename Setter >
	ToGP< Setter > to_GP( Setter& setter ) {
		ToGP< Setter > result( setter ) ;
		result.add_impl( ePerUnorderedHaplotype, eAlleleIndex, GTToGPUnphased< Setter >::create() ) ;
		result.add_impl( ePerOrderedHaplotype, eAlleleIndex, GTToGPPhased< Setter >::create() ) ;
		result.add_impl( ePerUnorderedGenotype, eProbability, GPToGPUnphased< Setter >::create() ) ;
		return result ;
	}
	
}
