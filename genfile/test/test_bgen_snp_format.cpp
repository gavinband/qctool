
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
		oStream.write( a_allele.data(), a_allele.size() ) ;
		oStream.write( b_allele.data(), b_allele.size() ) ;
		
		for( std::size_t i = 0; i < number_of_samples; ++i ) {
			// We store the probability
			// ((i*100)+o) / 10000.0
			uint16_t
				AA = get_probs( i, 0 ) * 10000,
				AB = get_probs( i, 1 ) * 10000,
				BB = get_probs( i, 2 ) * 10000
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
				AA = get_probs( i, 0 ) * 32768,
				AB = get_probs( i, 1 ) * 32768,
				BB = get_probs( i, 2 ) * 32768
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
		std::size_t const upper = std::round( total_fractional_part ) ;
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

	void write_probs( uint64_t** p, std::size_t* offset, uint64_t* const end_p, double p1, double p2, std::size_t number_of_bits ) {
#if DEBUG > 1
		std::cerr << "write_probs(): p1 = " << p1 << ", p2 = " << p2 << ".\n" ;
#endif
		assert( number_of_bits <= 64 ) ;
		assert( p1 >= 0.0 ) ;
		assert( p1 <= 1.0 ) ;
		assert( p2 >= 0.0 ) ;
		assert( p2 <= 1.0 ) ;
		assert( p1+p2 <= 1.0 ) ;
		double const scale = uint64_t( 0xFFFFFFFFFFFFFFFF ) >> ( 64 - number_of_bits ) ; // 2^bits - 1.
		double ps[3] ;
		ps[0] = p1 ;
		ps[1] = p2 ;
		ps[2] = ( 1 - p1 - p2 ) ;
		round_probs_to_simplex( &ps[0], 3, number_of_bits ) ;
		for( std::size_t i = 0; i < 2; ++i ) {
#if DEBUG > 1
			std::cerr << "write_probs(): p = " << std::hex << (*p)
				<< std::dec << ", offset = " << (*offset) << ", end_p = " << (end_p)
				<< ", number_of_bits = " << number_of_bits
				<< ", value = " << ps[i]
				<< ", unscaled value = " << uint64_t( ps[i] * scale )
				<< ".\n";
			std::cerr << "write_probs(): *p = " << genfile::string_utils::to_hex( reinterpret_cast< unsigned char* >( *p ), reinterpret_cast< unsigned char* >( *p ) + 8 ) << ".\n" ;
#endif
			uint64_t const storedValue = uint64_t( ps[i] * scale ) ;
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
		boost::function< double ( std::size_t i, std::size_t g ) > get_probs
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
		uint32_t const uncompressed_data_size = 10 + number_of_samples + ( number_of_samples * 2 * bits_per_probability ) ;
		genfile::write_little_endian_integer( oStream, uncompressed_data_size ) ;

		// Write number of samples and ploidies.
		genfile::write_little_endian_integer( oStream, number_of_samples ) ;
		genfile::write_little_endian_integer( oStream, number_of_alleles ) ;
		genfile::write_little_endian_integer( oStream, uint8_t( 2 ) ) ;
		genfile::write_little_endian_integer( oStream, uint8_t( 2 ) ) ;
		for( std::size_t i = 0; i < number_of_samples; ++i ) {
			uint8_t ploidy = 2 ;
			genfile::write_little_endian_integer( oStream, ploidy ) ;
		}

		genfile::write_little_endian_integer( oStream, uint8_t( 0 ) ) ;
		genfile::write_little_endian_integer( oStream, uint8_t( bits_per_probability ) ) ;

		{
			// Construct and write probability data.
			std::size_t const numberOfBytes = std::ceil( 3.0 * number_of_samples * bits_per_probability / 8.0 ) ;
			std::vector< char > probability_data( std::ceil( 3.0 * number_of_samples * bits_per_probability / 64.0 ) * 8, 0 ) ;
			{
				uint64_t const two_to_the_bits = ( uint64_t( 1 ) << bits_per_probability ) ;
				double scale = two_to_the_bits - 1 ;
				uint64_t* p = reinterpret_cast< uint64_t* >( &probability_data[0] ) ;
				uint64_t* const end_p = reinterpret_cast< uint64_t* const >( &probability_data[0] + probability_data.size() ) ;
				std::size_t offset = 0 ;
				for( std::size_t i = 0; i < number_of_samples; ++i ) {
					double const AA = get_probs( i, 0 ) ;
					double const AB = get_probs( i, 1 ) ;
#if DEBUG
					std::cerr << ( boost::format( "sample %d, bits_per_probability = %d, two_to_the_bits=%d, scale = %f, AA=%f, AB=%f, sum = %f\n" )
						% i % bits_per_probability % two_to_the_bits % scale % AA % AB % (AA+AB) ).str() ;
#endif
					write_probs( &p, &offset, end_p, AA, AB, bits_per_probability ) ;
				}
			}
			
			oStream.write( &probability_data[0], std::ceil( 3.0 * number_of_samples * bits_per_probability / 8.0 ) ) ;
		}
		
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
		boost::function< double ( std::size_t i, std::size_t g ) > get_probs
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
				get_probs
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
	enum State { eNone, eSetSample, eSetNumberOfEntries, eSetValue } ;

	ProbabilitySetter(
		std::size_t n,
		GetExpectedProbs get_expected_probs
	):
		m_number_of_samples( n ),
		m_get_expected_probs( get_expected_probs ),
		m_sample_i( std::numeric_limits< std::size_t >::max() ),
		m_entry_i( std::numeric_limits< std::size_t >::max() ),
		m_state( eNone )
	{}

	~ProbabilitySetter() throw() {
		TEST_ASSERT( m_sample_i < 1000 ) ;
		TEST_ASSERT( m_entry_i < 1000 ) ;
	}

	void set_number_of_samples( std::size_t n ) {
		BOOST_CHECK_EQUAL( n, m_number_of_samples ) ;
	}
	void set_number_of_alleles( std::size_t n ) {
		TEST_ASSERT( n == 2 ) ;
	}
	void set_sample( std::size_t i ) {
		TEST_ASSERT( i < m_number_of_samples ) ;
		TEST_ASSERT( m_state == eNone || m_state == eSetValue ) ;
		m_sample_i = i ;
		m_state = eSetSample ;
	}
	void set_order_type( OrderType const order_type, ValueType const value_type ) {
		TEST_ASSERT( order_type == ePerUnorderedGenotype ) ;
		TEST_ASSERT( value_type == eProbability ) ;
	}
	void set_number_of_entries( std::size_t n ) {
		TEST_ASSERT( n == 3 ) ;
		TEST_ASSERT( m_state == eSetSample ) ;
		m_state = eSetNumberOfEntries ;
		m_entry_i = 0 ;
	}
	void set_order_type( OrderType const type ) {
		TEST_ASSERT( type == eOrderedList ) ;
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

		TEST_ASSERT( m_entry_i < 3 ) ;
#if DEBUG > 2
		std::cerr << ( boost::format( "ProbabilitySetter: sample %d, entry %d.\n" ) % m_sample_i % m_entry_i ).str() ;
#endif
		BOOST_CHECK_CLOSE( value, m_get_expected_probs( m_sample_i, m_entry_i ), 0.00000001 ) ;
		++m_entry_i ;
	}

private:
	std::size_t m_number_of_samples ;
	GetExpectedProbs m_get_expected_probs ;
	std::size_t m_sample_i ;
	std::size_t m_entry_i ;
	std::set< std::pair< std::size_t, std::size_t > > m_set_values ;

	State m_state ;
} ;

double get_input_probability(
	std::size_t const number_of_samples,
	std::size_t i,
	std::size_t g
) {
	double x = double(i) / double(number_of_samples-1) ;
	if( g == 1 ) {
		// values of second prob are 1 minus 1/3rd of first prob.
		x = 0.25 * ( 1.0 - x ) ;
	} else if( g == 2 ) {
		// values of second prob are 1 minus 2/3rd of first prob.
		x = 0.75 * ( 1.0 - x ) ;
	}
	return x ;
}

double get_expected_stored_probability(
	std::size_t const number_of_samples,
	std::size_t i,
	std::size_t g,
	std::string const& bgen_version,
	std::size_t bits_per_probability
) {
	if( bgen_version == "v10" ) {
		return std::floor( get_input_probability( number_of_samples, i, g ) * 10000.0 ) / 10000.0 ;
	} else if( bgen_version == "v11" ) {
		return std::floor( get_input_probability( number_of_samples, i, g ) * 32768.0 ) / 32768.0 ;
	} else if( bgen_version == "v12" ){
		double v[3] ;
		for( std::size_t l = 0; l < 3; ++l ) {
			v[l] = get_input_probability( number_of_samples, i, l ) ;
		}
#if DEBUG > 1
		std::cerr << ( boost::format( "get_expected_stored_probability(): probs are: %f, %f, %f\n" ) % v[0] % v[1] % v[2] ).str() ;
#endif
		data::round_probs_to_simplex( &v[0], 3, bits_per_probability ) ;
#if DEBUG > 1
		std::cerr << ( boost::format( "get_expected_stored_probability(): expected probs are: %f, %f, %f\n" ) % v[0] % v[1] % v[2] ).str() ;
#endif
		return *(v+g) ;
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
		std::size_t bits_per_probability = 16
) {
	uint32_t flags = 0 ;
	boost::function< double ( std::size_t i, std::size_t g ) > get_probs ;
	if( bgen_version == "v11" ) {
		flags = e_v11Layout ;
	} else if( bgen_version == "v12" ) {
		flags = e_v12Layout ;
	} else {
		flags = e_v10Layout ;
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
			boost::bind( &get_input_probability, number_of_individuals, _1, _2 )
		)
	) ;
#if DEBUG
	std::cerr << "do_snp_block_read_test(): bgen_version=" << bgen_version << ".\n" ;
	std::cerr << "do_snp_block_read_test():        data is: " << genfile::string_utils::to_hex( inStream.str() ) << "\n" ;
//	std::cerr << "do_snp_block_read_test(): data (char) is: " << genfile::string_utils::to_hex_char( inStream.str() ) << "\n" ;
#endif
	
	std::string SNPID2 ;
	std::string RSID2 ;
	unsigned char chromosome_enum ;
	genfile::Chromosome chromosome2 ;
	uint32_t SNP_position2 ;
	std::string a2 ;
	std::string b2 ;

	std::vector< char > buffer1, buffer2 ;

	genfile::bgen::read_snp_identifying_data(
		inStream,
		flags,
		&SNPID2,
		&RSID2,
		&chromosome_enum,
		&SNP_position2,
		&a2,
		&b2
	) ;
	chromosome2 = genfile::Chromosome( chromosome_enum ) ;

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
			bits_per_probability
		)
	) ;

	genfile::bgen::read_snp_probability_data(
		inStream,
		flags,
		number_of_individuals,
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
		std::size_t bits_per_probability = 16
) {
	uint32_t flags = 0 ;
	if( bgen_version == "v11" ) {
		flags = e_v11Layout ;
	} else if( bgen_version == "v12" ) {
		flags = e_v12Layout ;
	} else {
		flags = e_v10Layout ;
	}
	
	std::ostringstream outStream ;
	
	genfile::bgen::write_snp_identifying_data( 
		outStream,
		flags,
		number_of_individuals,
		std::max( SNPID.size(), RSID.size() ) + 1,
		SNPID,
		RSID,
		chromosome,
		SNP_position,
		a,
		b
	) ;

	genfile::bgen::write_snp_probability_data( 
		outStream,
		flags,
		number_of_individuals,
		boost::bind( &get_input_probability, number_of_individuals, _1, 0 ),
		boost::bind( &get_input_probability, number_of_individuals, _1, 1 ),
		boost::bind( &get_input_probability, number_of_individuals, _1, 2 )
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
		boost::bind( &get_input_probability, number_of_individuals, _1, _2 )
	) ;

#if DEBUG	
	std::cerr << "          actual: " << genfile::string_utils::to_hex( outStream.str() ) << "\n" ;
//	std::cerr << "   actual (char): " << genfile::string_utils::to_hex_char( outStream.str() ) << "\n" ;
	std::cerr << "        expected: " << genfile::string_utils::to_hex( expected ) << "\n" ;
//	std::cerr << " expected (char): " << genfile::string_utils::to_hex_char( expected ) << "\n" ;
#endif
		
	TEST_ASSERT( outStream.str() == expected ) ;
}


// Test that reading a properly formatted bgen snp block works.
AUTO_TEST_CASE( test_snp_block_input ) {
	std::cout << "test_snp_block_input\n" ;
	do_snp_block_read_test( "v10", 6, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C" ) ;
	do_snp_block_read_test( "v10", 100, "01234567890123456789012345678901234567890123456789", "01234567890123456789012345678901234567890123456789", genfile::Chromosome22, 4294967295u, "G", "T" ) ;
	do_snp_block_read_test( "v11", 6, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C" ) ;
	do_snp_block_read_test( "v11", 100, "01234567890123456789012345678901234567890123456789", "01234567890123456789012345678901234567890123456789", genfile::Chromosome22, 4294967295u, "G", "T" ) ;
	for( std::size_t number_of_bits = 1; number_of_bits <= 32; ++number_of_bits ) {
		do_snp_block_read_test( "v12", 2, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits ) ;
		do_snp_block_read_test( "v12", 6, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits ) ;
		do_snp_block_read_test( "v12", 15, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits ) ;
		do_snp_block_read_test( "v12", 37, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits ) ;
		do_snp_block_read_test( "v12", 100, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits ) ;
		do_snp_block_read_test( "v12", 1000, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C", number_of_bits ) ;
	}
}

// Test that writing a bgen snp block gives properly formatted output.
AUTO_TEST_CASE( test_snp_block_output ) {
	std::cout << "test_snp_block_output\n" ;
	do_snp_block_write_test( "v10", 6, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C" ) ;
	do_snp_block_write_test( "v10", 100, "01234567890123456789012345678901234567890123456789", "01234567890123456789012345678901234567890123456789", genfile::Chromosome22, 4294967295u, "G", "T" ) ;
	do_snp_block_write_test( "v11", 6, "SNP01", "RS01", genfile::Chromosome1, 1000001, "A", "C" ) ;
	do_snp_block_write_test( "v11", 100, "01234567890123456789012345678901234567890123456789", "01234567890123456789012345678901234567890123456789", genfile::Chromosome22, 4294967295u, "G", "T" ) ;
}
