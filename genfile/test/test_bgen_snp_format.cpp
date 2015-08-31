
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cassert>
#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <boost/function.hpp>
#include "test_case.hpp"
#include "genfile/bgen/bgen.hpp"
#include "genfile/string_utils/hex.hpp"
#include "stdint.h"

#define DEBUG 1

BOOST_AUTO_TEST_SUITE( test_bgen )

// The following section contains a simple snp block writer.
namespace data {
	std::string construct_snp_block_v10(
		uint32_t number_of_samples,
		unsigned char max_id_size,
		std::string SNPID,
		std::string RSID,
		genfile::Chromosome chromosome,
		uint32_t SNP_position,
		std::string a_allele,
		std::string b_allele,
		boost::function< double ( std::size_t i, std::size_t g ) > get_probs
	) {
		std::ostringstream oStream ;
		genfile::write_little_endian_integer( oStream, number_of_samples ) ;
		genfile::write_little_endian_integer( oStream, max_id_size ) ;
		genfile::write_little_endian_integer( oStream, static_cast< char >( SNPID.size() )) ;
		oStream.write( SNPID.data(), SNPID.size() ) ;
		oStream.write( "                ", max_id_size - SNPID.size()) ;
		genfile::write_little_endian_integer( oStream, static_cast< char >( RSID.size() )) ;
		oStream.write( RSID.data(), RSID.size() ) ;
		oStream.write( "                ", max_id_size - RSID.size()) ;
		unsigned char chr = chromosome ;
		genfile::write_little_endian_integer( oStream, chr ) ;
		genfile::write_little_endian_integer( oStream, SNP_position ) ;

		assert( a_allele.size() == 1 ) ;
		assert( b_allele.size() == 1 ) ;
		unsigned char const a_allele_size = a_allele.size() ;
		unsigned char const b_allele_size = a_allele.size() ;
		oStream.put( a_allele_size ) ;
		oStream.write( a_allele.data(), a_allele.size() ) ;
		oStream.put( b_allele_size ) ;
		oStream.write( b_allele.data(), b_allele.size() ) ;
		
		for( std::size_t i = 0; i < number_of_samples; ++i ) {
			// We store the probability
			// ((i*100)+o) / 10000.0
			uint16_t
				AA = std::floor( 0.5 + get_probs( i, 0 ) * 10000 ),
				AB = std::floor( 0.5 + get_probs( i, 1 ) * 10000 ),
				BB = std::floor( 0.5 + get_probs( i, 2 ) * 10000 )
			;
			genfile::write_little_endian_integer( oStream, AA ) ;
			genfile::write_little_endian_integer( oStream, AB ) ;
			genfile::write_little_endian_integer( oStream, BB ) ;
		}

		return oStream.str() ;
	}

	std::string construct_snp_block_v11(
		uint32_t number_of_samples,
		std::string SNPID,
		std::string RSID,
		genfile::Chromosome chromosome,
		uint32_t SNP_position,
		std::string a_allele,
		std::string b_allele,
		boost::function< double ( std::size_t i, std::size_t g ) > get_probs
	) {
		std::ostringstream oStream ;
		genfile::write_little_endian_integer( oStream, number_of_samples ) ;
		genfile::write_little_endian_integer( oStream, static_cast< uint16_t >( SNPID.size() )) ;
		oStream.write( SNPID.data(), SNPID.size() ) ;
		genfile::write_little_endian_integer( oStream, static_cast< uint16_t >( RSID.size() )) ;
		oStream.write( RSID.data(), RSID.size() ) ;
		std::string chr( chromosome ) ;
		genfile::write_little_endian_integer( oStream, static_cast< uint16_t >( chr.size() ) ) ;
		oStream.write( chr.data(), chr.size() ) ;
		genfile::write_little_endian_integer( oStream, SNP_position ) ;

		genfile::write_little_endian_integer( oStream, static_cast< uint32_t >( a_allele.size() )) ;
		oStream.write( a_allele.data(), a_allele.size() ) ;
		genfile::write_little_endian_integer( oStream, static_cast< uint32_t >( b_allele.size() )) ;
		oStream.write( b_allele.data(), b_allele.size() ) ;
		
		for( std::size_t i = 0; i < number_of_samples; ++i ) {
			uint16_t
				AA = std::floor( 0.5 + get_probs( i, 0 ) * 32768.0 ),
				AB = std::floor( 0.5 + get_probs( i, 1 ) * 32768.0 ),
				BB = std::floor( 0.5 + get_probs( i, 2 ) * 32768.0 )
			;
			genfile::write_little_endian_integer( oStream, AA ) ;
			genfile::write_little_endian_integer( oStream, AB ) ;
			genfile::write_little_endian_integer( oStream, BB ) ;
		}

		return oStream.str() ;
	}
	
	void round_probs_to_simplex( double *p, std::size_t n, std::size_t number_of_bits ) {
		double const scale = uint64_t( 0xFFFFFFFFFFFFFFFF ) >> ( 64 - number_of_bits ) ;
		std::multimap< double, double*, std::greater< double > > fractional_parts ;
		double total_fractional_part = 0.0 ;
		for( std::size_t i = 0; i < n; ++i ) {
			*(p+i) *= scale ;
			double const fractional_part = *(p+i) - std::floor(*(p+i)) ;
			fractional_parts.insert( std::make_pair( fractional_part, (p+i) ) ) ;
			total_fractional_part += fractional_part ;
		}
		std::size_t const upper = std::floor( 0.5 + total_fractional_part ) ;
#if DEBUG > 2
		std::cerr << "round_probs_to_simplex(): number_of_bits = " << number_of_bits << ", scale = " << scale << ", total_fractional_part = " << total_fractional_part << ", upper = " << upper << ".\n" ;
		std::cerr << "round_probs_to_simplex(): p1 = " << *p << ".\n" ;
#endif
		std::multimap< double, double*, std::greater< double > >::const_iterator
			i = fractional_parts.begin(),
			end_i = fractional_parts.end() ;
		for( std::size_t count = 0; i != end_i; ++i, ++count ) {
			if( count < upper ) {
				*(i->second) = std::ceil( *(i->second) ) / scale ;
			} else {
				*(i->second) = std::floor( *(i->second) ) / scale ;
			}
		}
	}

	void write_probs( uint64_t** p, std::size_t* offset, uint64_t* const end_p, double const* probs, std::size_t const n, std::size_t number_of_bits ) {
		double const scale = uint64_t( 0xFFFFFFFFFFFFFFFF ) >> ( 64 - number_of_bits ) ; // 2^bits - 1.
		for( std::size_t i = 0; i < (n-1); ++i ) {
#if DEBUG > 1
			std::cerr << "write_probs(): p = " << std::hex << (*p)
				<< std::dec << ", offset = " << (*offset) << ", end_p = " << (end_p)
				<< ", number_of_bits = " << number_of_bits
				<< ", value = " << probs[i]
				<< ", unscaled value = " << uint64_t( probs[i] * scale )
				<< ".\n";
#endif
			uint64_t const storedValue = uint64_t( probs[i] * scale ) ;
			**p |= storedValue << *offset ;
			(*offset) += number_of_bits ;
			if( (*offset) >= 64 ) {
				// move to next 64-bit word
				(*offset) -= 64 ;
				++(*p) ;
				// Make sure and store unwritten remainder of value...
				(**p) |= ( storedValue >> ( number_of_bits - *offset ) ) ;
			}
		}
	}
	
	std::string construct_snp_block_v12(
		uint32_t number_of_samples,
		std::string SNPID,
		std::string RSID,
		genfile::Chromosome chromosome,
		uint32_t SNP_position,
		std::string a_allele,
		std::string b_allele,
		std::size_t const bits_per_probability,
		boost::function< double ( std::size_t i, std::size_t g ) > get_probs,
		std::string const& type
	) {
		assert( bits_per_probability <= 64 ) ;

		std::ostringstream oStream ;
		genfile::write_little_endian_integer( oStream, static_cast< uint16_t >( SNPID.size() )) ;
		oStream.write( SNPID.data(), SNPID.size() ) ;
		genfile::write_little_endian_integer( oStream, static_cast< uint16_t >( RSID.size() )) ;
		oStream.write( RSID.data(), RSID.size() ) ;
		std::string chr( chromosome ) ;
		genfile::write_little_endian_integer( oStream, static_cast< uint16_t >( chr.size() ) ) ;
		oStream.write( chr.data(), chr.size() ) ;
		genfile::write_little_endian_integer( oStream, SNP_position ) ;

		uint16_t const number_of_alleles = 2 ;
		genfile::write_little_endian_integer( oStream, number_of_alleles ) ;
		genfile::write_little_endian_integer( oStream, static_cast< uint32_t >( a_allele.size() )) ;
		oStream.write( a_allele.data(), a_allele.size() ) ;
		genfile::write_little_endian_integer( oStream, static_cast< uint32_t >( b_allele.size() )) ;
		oStream.write( b_allele.data(), b_allele.size() ) ;
		
		// Total size of probability data is:
		uint32_t const stored_probability_size = (((number_of_samples*2*bits_per_probability)+7)/8) ;
		// To this we must add space for the header block
		uint32_t const buffer_size = 10 + number_of_samples + stored_probability_size ;
		// And four bytes indicating the total length.
		genfile::write_little_endian_integer( oStream, buffer_size ) ;

		// Write number of samples and ploidies.
		genfile::write_little_endian_integer( oStream, number_of_samples ) ;
		genfile::write_little_endian_integer( oStream, number_of_alleles ) ;
		genfile::write_little_endian_integer( oStream, uint8_t( 2 ) ) ;
		genfile::write_little_endian_integer( oStream, uint8_t( 2 ) ) ;
		for( std::size_t i = 0; i < number_of_samples; ++i ) {
			uint8_t ploidy = 2 ;
			genfile::write_little_endian_integer( oStream, ploidy ) ;
		}

		uint8_t const phased = ( type == "phased" ) ? 1 : 0 ;
		genfile::write_little_endian_integer( oStream, phased ) ;
		genfile::write_little_endian_integer( oStream, uint8_t( bits_per_probability ) ) ;

		uint64_t const two_to_the_bits = ( uint64_t( 1 ) << bits_per_probability ) ;
		double scale = two_to_the_bits - 1 ;
		std::vector< char > probability_data( std::ceil( 2.0 * number_of_samples * bits_per_probability / 64.0 ) * 8, 0 ) ;
		uint64_t* const buffer = reinterpret_cast< uint64_t* >( &probability_data[0] ) ;
		uint64_t* const end = reinterpret_cast< uint64_t* const >( &probability_data[0] + probability_data.size() ) ;
		uint64_t* p = buffer ;
		std::size_t offset = 0 ;
		if( type == "unphased" ) {
			// Construct and write probability data.
			{
				double probs[3] ;
				for( std::size_t i = 0; i < number_of_samples; ++i ) {
					probs[0] = get_probs( i, 0 ) ;
					probs[1] = get_probs( i, 1 ) ;
					probs[2] = get_probs( i, 2 ) ;
					
					round_probs_to_simplex( &probs[0], 3, bits_per_probability ) ;
#if DEBUG
					std::cerr << ( boost::format( "sample %d of %d, bits_per_probability = %d, two_to_the_bits=%d, scale = %f, AA=%f, AB=%f, sum = %f\n" )
						% i % number_of_samples % bits_per_probability % two_to_the_bits % scale % probs[0] % probs[1] % (probs[0]+probs[1]) ).str() ;
#endif
					write_probs( &p, &offset, end, probs, 3, bits_per_probability ) ;
				}
			}
		} else if( type == "phased" ) {
			// Construct and write probability data.
			double probs[2] ;
			for( std::size_t i = 0; i < number_of_samples; ++i ) {
				for( std::size_t hap = 0; hap < 2; ++hap ) {
					probs[0] = get_probs( i, 0+(2*hap) ) ;
					probs[1] = get_probs( i, 1+(2*hap) ) ;
					round_probs_to_simplex( &probs[0], 2, bits_per_probability ) ;
#if DEBUG
					std::cerr << ( boost::format( "sample %d of %d, hap %d, bits_per_probability = %d, two_to_the_bits=%d, scale = %f, AA=%f, AB=%f, sum = %f\n" )
						% i % number_of_samples % hap % bits_per_probability % two_to_the_bits % scale % probs[0] % probs[1] % (probs[0]+probs[1]) ).str() ;
#endif
					write_probs( &p, &offset, end, &probs[0], 2, bits_per_probability ) ;
				}
			}
		} else {
			assert(0) ;
		}
		std::size_t const numBytes = ((p-buffer)*8) + ((offset+7)/8) ;
		BOOST_CHECK_EQUAL( numBytes, stored_probability_size ) ;
		oStream.write( &probability_data[0], numBytes ) ;
		
		return oStream.str() ;
	}
	
	std::string construct_snp_block(
		std::string const& version,
		uint32_t number_of_samples,
		std::string SNPID,
		std::string RSID,
		genfile::Chromosome chromosome,
		uint32_t SNP_position,
		std::string a_allele,
		std::string b_allele,
		std::size_t const bits_per_probability,
		boost::function< double ( std::size_t i, std::size_t g ) > get_probs,
		std::string const& type
	) {
#if DEBUG
		std::cerr << "construct_snp_block(): version=" << version << ".\n" ;
#endif
		if( version != "v12" ) {
			assert( bits_per_probability == 16 ) ;
		}

		if( version == "v10" ) {
			return construct_snp_block_v10(
				number_of_samples,
				std::max( SNPID.size(), RSID.size() ) + 1,
				SNPID,
				RSID,
				chromosome,
				SNP_position,
				a_allele,
				b_allele,
				get_probs
			) ;
		} else if( version == "v11" ) {
			return construct_snp_block_v11(
				number_of_samples,
				SNPID,
				RSID,
				chromosome,
				SNP_position,
				a_allele,
				b_allele,
				get_probs
			) ;
		} else if( version == "v12" ) {
			return construct_snp_block_v12(
				number_of_samples,
				SNPID,
				RSID,
				chromosome,
				SNP_position,
				a_allele,
				b_allele,
				bits_per_probability,
				get_probs,
				type
			) ;
		} else {
			assert(0) ;
		}
	}
}

namespace {
	enum { e_v10Layout = 0, e_v11Layout = 0x4, e_v12Layout = 0x8 } ;
	enum { e_CompressedSNPBlocks = 1 } ;
}
// The following section defines the needed objects for use with the bgen.hpp implementation.
template< typename T >
struct Setter
{
	Setter( T& field ): m_field( field ) {} ;
	template< typename T2 >
	void operator()( T2 const& value ) { m_field = T(value) ; }
private:
	T& m_field ;
} ;

template< typename T >
Setter< T > make_setter( T& field ) { return Setter<T>( field ) ; }


struct probabilities {
	double AA, AB, BB ;
} ;

struct ProbabilitySetter: public genfile::VariantDataReader::PerSampleSetter {
	typedef boost::function< double( std::size_t i, std::size_t g ) > GetExpectedProbs ;
	enum State { eNone, eSetNumberOfSamples, eSetSample, eSetNumberOfEntries, eSetValue } ;

	ProbabilitySetter(
		std::size_t n,
		GetExpectedProbs get_expected_probs
	):
		m_number_of_samples( n ),
		m_get_expected_probs( get_expected_probs ),
		m_sample_i( std::numeric_limits< std::size_t >::max() ),
		m_number_of_entries( std::numeric_limits< std::size_t >::max() ),
		m_entry_i( std::numeric_limits< std::size_t >::max() ),
		m_state( eNone )
	{}

	~ProbabilitySetter() throw() {
		BOOST_CHECK_EQUAL( m_sample_i + 1, m_number_of_samples ) ;
		BOOST_CHECK_EQUAL( m_entry_i, m_number_of_entries ) ;
	}
	void set_number_of_samples( std::size_t nSamples, std::size_t nAlleles ) {
		BOOST_CHECK_EQUAL( m_state, eNone ) ;
		BOOST_CHECK_EQUAL( nSamples, m_number_of_samples ) ;
		BOOST_CHECK_EQUAL( nAlleles, 2 ) ;
		m_state = eSetNumberOfSamples ;
	}
	bool set_sample( std::size_t i ) {
		TEST_ASSERT( i < m_number_of_samples ) ;
		TEST_ASSERT( m_state == eSetNumberOfSamples || m_state == eSetValue ) ;
		BOOST_CHECK_EQUAL( m_entry_i, m_number_of_entries ) ;
		m_sample_i = i ;
		m_state = eSetSample ;
		return true ;
	}
	void set_number_of_entries( std::size_t n, OrderType const order_type, ValueType const value_type ) {
#if DEBUG > 2
		std::cerr << "ProbabilitySetter::set_number_of_entries(): n = " << n
			<< ", order_type = " << order_type << ".\n" ;
#endif
		BOOST_CHECK_EQUAL( m_state, eSetSample ) ;
		TEST_ASSERT(
			( order_type == genfile::ePerUnorderedGenotype && m_number_of_entries == 3)
			||
			( order_type == genfile::ePerPhasedHaplotypePerAllele && m_number_of_entries == 4)
		) ;
		BOOST_CHECK_EQUAL( value_type, genfile::eProbability ) ;
		m_order_type = order_type ;
		m_state = eSetNumberOfEntries ;
		m_number_of_entries = n ;
		m_entry_i = 0 ;
	}
	void operator()( genfile::MissingValue const value ) {
		TEST_ASSERT(0) ;
	}

	void operator()( std::string& value ) {
		TEST_ASSERT(0) ;
	}

	void operator()( Integer const value ) {
		TEST_ASSERT(0) ;
	}

	void operator()( double const value ) {
		TEST_ASSERT(
			(m_entry_i == 0 && m_state == eSetNumberOfEntries)
			 || (m_entry_i > 0 && m_state == eSetValue)
		) ;
		m_state = eSetValue ;

		TEST_ASSERT( m_entry_i < m_number_of_entries ) ;
#if DEBUG > 2
		std::cerr << ( boost::format( "ProbabilitySetter: sample %d, entry %d of %d.\n" ) % m_sample_i % m_entry_i % m_number_of_entries ).str() ;
#endif
		BOOST_CHECK_CLOSE( value, m_get_expected_probs( m_sample_i, m_entry_i ), 0.00000001 ) ;
		++m_entry_i ;
	}

private:
	std::size_t m_number_of_samples ;
	GetExpectedProbs m_get_expected_probs ;
	std::size_t m_sample_i ;
	std::size_t m_number_of_entries ;
	std::size_t m_entry_i ;
	OrderType m_order_type ;
	std::set< std::pair< std::size_t, std::size_t > > m_set_values ;

	State m_state ;
} ;

double get_input_probability(
	std::size_t const number_of_samples,
	std::size_t i,
	std::size_t g,
	std::string const& type = "unphased"
) {
	double x = 0 ;
	if( type == "phased" ) {
		assert( g < 4 ) ;
		// two haplotypes, each of whose probs sum to one
		if( g == 0 ) {
			x = double(i) / double(number_of_samples-1) ;
		}
		else if( g == 1 ) {
			x = 1.0 - (double(i) / double(number_of_samples-1)) ; ;
		} else if( g == 2 ) {
			x = double((i+1) % number_of_samples) / double(number_of_samples-1) ;
		} else if( g == 3 ) {
			x = 1.0 - double((i+1) % number_of_samples) / double(number_of_samples-1) ;
		}
	} else {
		x = double(i) / double(number_of_samples-1) ;
		assert( g < 3 ) ;
		if( g == 1 ) {
			// values of second prob are 1 minus 1/3rd of first prob.
			x = 0.25 * ( 1.0 - x ) ;
		} else if( g == 2 ) {
			// values of second prob are 1 minus 2/3rd of first prob.
			x = 0.75 * ( 1.0 - x ) ;
		}
	}
	return x ;
}

double get_expected_stored_probability(
	std::size_t const number_of_samples,
	std::size_t i,
	std::size_t g,
	std::string const& bgen_version,
	std::size_t bits_per_probability,
	std::string const& type = "unphased"
) {
	if( bgen_version == "v10" ) {
		return std::floor( 0.5 + get_input_probability( number_of_samples, i, g ) * 10000.0 ) / 10000.0 ;
	} else if( bgen_version == "v11" ) {
		return std::floor( 0.5 + get_input_probability( number_of_samples, i, g ) * 32768.0 ) / 32768.0 ;
	} else if( bgen_version == "v12" ){
		double v[4] ;
		if( type == "phased" ) {
			for( std::size_t l = 0; l < 4; ++l ) {
				v[l] = get_input_probability( number_of_samples, i, l, type ) ;
			}
			data::round_probs_to_simplex( &v[0]+(2*(g/2)), 2, bits_per_probability ) ;
#if DEBUG > 1
			std::cerr << ( boost::format( "get_expected_stored_probability(): expected probs are: %f, %f\n" ) % v[(2*(g/2))] % v[1+(2*(g/2))] ).str() ;
#endif
			return *(v+g) ;
		} else {
			double v[3] ;
			for( std::size_t l = 0; l < 3; ++l ) {
				v[l] = get_input_probability( number_of_samples, i, l, type ) ;
			}
			data::round_probs_to_simplex( &v[0], 3, bits_per_probability ) ;
#if DEBUG > 1
			std::cerr << ( boost::format( "get_expected_stored_probability(): expected probs are: %f, %f, %f\n" ) % v[0] % v[1] % v[2] ).str() ;
#endif
			return *(v+g) ;
		}
	} else {
		assert(0) ;
	}
}

// The following section contains the main tests.
void do_snp_block_read_test( 
		std::string const& bgen_version,
		uint32_t number_of_individuals,
		std::string SNPID,
		std::string RSID,
		genfile::Chromosome chromosome,
		uint32_t SNP_position,
		std::string a,
		std::string b,
		std::size_t bits_per_probability = 16,
		std::string const& type = "unphased"
) {
	genfile::bgen::Context context ;
	context.number_of_samples = number_of_individuals ;
	boost::function< double ( std::size_t i, std::size_t g ) > get_probs ;
	if( bgen_version == "v11" ) {
		context.flags = e_v11Layout ;
	} else if( bgen_version == "v12" ) {
		context.flags = e_v12Layout ;
	} else {
		context.flags = e_v10Layout ;
	}
	

	std::istringstream inStream ;
	inStream.str(
		data::construct_snp_block(
			bgen_version,
			number_of_individuals,
			SNPID,
			RSID,
			chromosome,
			SNP_position,
			a,
			b,
			bits_per_probability,
			boost::bind( &get_input_probability, number_of_individuals, _1, _2, type ),
			type
		)
	) ;
#if DEBUG
	std::cerr << "do_snp_block_read_test(): bgen_version=" << bgen_version << ".\n" ;
	std::cerr << "do_snp_block_read_test():        data is: " << genfile::string_utils::to_hex( inStream.str() ) << "\n" ;
	std::cerr << "do_snp_block_read_test(): data (char) is: " << genfile::string_utils::to_hex_char( inStream.str() ) << "\n" ;
#endif
	
	std::string SNPID2 ;
	std::string RSID2 ;
	std::string chromosome_string ;
	genfile::Chromosome chromosome2 ;
	uint32_t SNP_position2 ;
	std::string a2 ;
	std::string b2 ;

	std::vector< char > buffer1, buffer2 ;

	genfile::bgen::read_snp_identifying_data(
		inStream,
		context,
		&SNPID2,
		&RSID2,
		&chromosome_string,
		&SNP_position2,
		&a2,
		&b2
	) ;
	chromosome2 = genfile::Chromosome( chromosome_string ) ;

	BOOST_CHECK_EQUAL( SNPID2, SNPID ) ;
	BOOST_CHECK_EQUAL( RSID2, RSID ) ;
	BOOST_CHECK_EQUAL( chromosome2, chromosome ) ;
	BOOST_CHECK_EQUAL( SNP_position2, SNP_position ) ;
	BOOST_CHECK_EQUAL( a2, a ) ;
	BOOST_CHECK_EQUAL( b2, b ) ;

	ProbabilitySetter setter(
		number_of_individuals,
		boost::bind(
			&get_expected_stored_probability,
			number_of_individuals,
			_1,
			_2,
			bgen_version,
			bits_per_probability,
			type
		)
	) ;

	genfile::bgen::read_and_parse_probability_data(
		inStream,
		context,
		setter,
		&buffer1,
		&buffer2
	) ;
}

void do_snp_block_write_test(
		std::string const& bgen_version,
		uint32_t number_of_individuals,
		std::string SNPID,
		std::string RSID,
		genfile::Chromosome chromosome,
		uint32_t SNP_position,
		std::string a,
		std::string b,
		std::size_t bits_per_probability = 16,
		std::string const& type = "unphased"
) {
	genfile::bgen::Context context ;
	context.number_of_samples = number_of_individuals ;
	if( bgen_version == "v11" ) {
		context.flags = e_v11Layout ;
	} else if( bgen_version == "v12" ) {
		context.flags = e_v12Layout ;
	} else {
		assert(0) ;	
	}
	
#if DEBUG
	std::cerr
		<< "do_snp_block_write_test(): bgen_version=" << bgen_version
		<< ", number_of_samples = " << number_of_individuals
		<< ", number_of_bits = " << bits_per_probability
		<< ".\n" ;
	std::cerr << "!! This test fails, phased output not implemented in bgen.hpp yet.\n" ;
#endif
	
	std::ostringstream outStream ;
	
	genfile::bgen::write_snp_identifying_data( 
		outStream,
		context,
		std::max( SNPID.size(), RSID.size() ) + 1,
		SNPID,
		RSID,
		chromosome,
		SNP_position,
		a,
		b
	) ;

	std::vector< char > buffer ;
	std::vector< char > buffer2 ;
	genfile::bgen::write_snp_probability_data( 
		outStream,
		context,
		boost::bind( &get_input_probability, number_of_individuals, _1, 0, type ),
		boost::bind( &get_input_probability, number_of_individuals, _1, 1, type ),
		boost::bind( &get_input_probability, number_of_individuals, _1, 2, type ),
		bits_per_probability,
		&buffer,
		&buffer2
	) ;	

	std::string expected = data::construct_snp_block(
		bgen_version,
		number_of_individuals,
		SNPID,
		RSID,
		chromosome,
		SNP_position,
		a,
		b,
		bits_per_probability,
		boost::bind( &get_input_probability, number_of_individuals, _1, _2, type ),
		type
	) ;

#if DEBUG	
	std::cerr << "          actual: " << genfile::string_utils::to_hex( outStream.str() ) << "\n" ;
	std::cerr << "   actual (char): " << genfile::string_utils::to_hex_char( outStream.str() ) << "\n" ;
	std::cerr << "        expected: " << genfile::string_utils::to_hex( expected ) << "\n" ;
	std::cerr << " expected (char): " << genfile::string_utils::to_hex_char( expected ) << "\n" ;
#endif
		
	TEST_ASSERT( outStream.str() == expected ) ;
}


// Test that reading a properly formatted bgen snp block works.
AUTO_TEST_CASE( test_snp_block_input_unphased ) {
	std::cout << "test_snp_block_input_unphased\n" ;
//	do_snp_block_read_test( "v10", 0, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C" ) ;
	do_snp_block_read_test( "v10", 6, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C" ) ;
	do_snp_block_read_test( "v10", 15, "01234567890123456789012345678901234567890123456789", "01234567890123456789012345678901234567890123456789", genfile::Chromosome22, 4294967295u, "G", "T" ) ;
	do_snp_block_read_test( "v10", 100, "01234567890123456789012345678901234567890123456789", "01234567890123456789012345678901234567890123456789", genfile::Chromosome22, 4294967295u, "G", "T" ) ;
	do_snp_block_read_test( "v10", 1001, "01234567890123456789012345678901234567890123456789", "01234567890123456789012345678901234567890123456789", genfile::Chromosome22, 4294967295u, "G", "T" ) ;

//	do_snp_block_read_test( "v11", 0, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C" ) ;
	do_snp_block_read_test( "v11", 6, "01234567890123456789012345678901234567890123456789", "01234567890123456789012345678901234567890123456789", genfile::Chromosome22, 4294967295u, "G", "T" ) ;
	do_snp_block_read_test( "v11", 15, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C" ) ;
	do_snp_block_read_test( "v11", 100, "01234567890123456789012345678901234567890123456789", "01234567890123456789012345678901234567890123456789", genfile::Chromosome22, 4294967295u, "G", "T" ) ;
	do_snp_block_read_test( "v11", 1001, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C" ) ;

	for( std::size_t number_of_bits = 1; number_of_bits <= 32; ++number_of_bits ) {
//		do_snp_block_read_test( "v12", 0, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits ) ;
		do_snp_block_read_test( "v12", 2, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits ) ;
		do_snp_block_read_test( "v12", 6, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits ) ;
		do_snp_block_read_test( "v12", 15, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits ) ;
		do_snp_block_read_test( "v12", 37, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits ) ;
		do_snp_block_read_test( "v12", 100, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits ) ;
		do_snp_block_read_test( "v12", 1001, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits ) ;
	}
}

AUTO_TEST_CASE( test_snp_block_input_phased ) {
	std::cout << "test_snp_block_input_phased\n" ;
	for( std::size_t number_of_bits = 1; number_of_bits <= 32; ++number_of_bits ) {
		do_snp_block_read_test( "v12", 0, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits, "phased" ) ;
		do_snp_block_read_test( "v12", 2, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits, "phased" ) ;
		do_snp_block_read_test( "v12", 6, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits, "phased" ) ;
		do_snp_block_read_test( "v12", 15, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits, "phased" ) ;
		do_snp_block_read_test( "v12", 37, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits, "phased" ) ;
		do_snp_block_read_test( "v12", 100, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits, "phased" ) ;
		do_snp_block_read_test( "v12", 1001, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits, "phased" ) ;
	}
}

// Test that writing a bgen snp block gives properly formatted output.
AUTO_TEST_CASE( test_snp_block_output_unphased ) {
	std::cout << "test_snp_block_output_unphased\n" ;
	do_snp_block_write_test( "v11", 6, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C" ) ;
	do_snp_block_write_test( "v11", 100, "01234567890123456789012345678901234567890123456789", "01234567890123456789012345678901234567890123456789", genfile::Chromosome22, 4294967295u, "G", "T" ) ;
	for( std::size_t number_of_bits = 1; number_of_bits <= 32; ++number_of_bits ) {
		do_snp_block_write_test( "v12", 0, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits ) ;
		do_snp_block_write_test( "v12", 2, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits ) ;
		do_snp_block_write_test( "v12", 6, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits ) ;
		do_snp_block_write_test( "v12", 15, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits ) ;
		do_snp_block_write_test( "v12", 37, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits ) ;
		do_snp_block_write_test( "v12", 100, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits ) ;
		do_snp_block_write_test( "v12", 1001, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits ) ;
	}
}

AUTO_TEST_CASE( test_snp_block_output_phased ) {
	std::cout << "test_snp_block_output_phased\n" ;
	for( std::size_t number_of_bits = 1; number_of_bits <= 32; ++number_of_bits ) {
		do_snp_block_write_test( "v12", 0, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits, "phased" ) ;
		do_snp_block_write_test( "v12", 2, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits, "phased" ) ;
		do_snp_block_write_test( "v12", 6, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits, "phased" ) ;
		do_snp_block_write_test( "v12", 15, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits, "phased" ) ;
		do_snp_block_write_test( "v12", 37, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits, "phased" ) ;
		do_snp_block_write_test( "v12", 100, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits, "phased" ) ;
		do_snp_block_write_test( "v12", 1001, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits, "phased" ) ;
	}
}

BOOST_AUTO_TEST_SUITE_END()
