#include <iostream>
#include <iomanip>
#include <sstream>
#include <cassert>
#include "test_case.hpp"
#include "genfile/bgen.hpp"
#include "stdint.h"

// The following section contains a simple snp block writer.
namespace data {
	std::string construct_snp_block(
		uint32_t number_of_samples,
		unsigned char max_id_size,
		std::string SNPID,
		std::string RSID,
		genfile::Chromosome chromosome,
		uint32_t SNP_position,
		char a_allele,
		char b_allele
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
		oStream.put( a_allele ) ;
		oStream.put( b_allele ) ;
		for( std::size_t i = 0; i < number_of_samples; ++i ) {
			uint16_t AA = i, AB = i, BB = i ;
			genfile::write_little_endian_integer( oStream, AA ) ;
			genfile::write_little_endian_integer( oStream, AB ) ;
			genfile::write_little_endian_integer( oStream, BB ) ;
		}

		return oStream.str() ;
	}
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

struct ProbabilitySetter {
	ProbabilitySetter( std::vector< probabilities >& probabilities ): m_probabilities( probabilities ) {}
	void operator() ( std::size_t i, double aa, double ab, double bb ) {
		m_probabilities[i].AA = aa ;  
		m_probabilities[i].AB = ab ;  
		m_probabilities[i].BB = bb ;  
	}

private:
	std::vector<probabilities>& m_probabilities ;
} ;


// The following section contains the main tests.
void do_snp_block_read_test( 
		uint32_t number_of_individuals,
		unsigned char ID_field_length,
		std::string SNPID,
		std::string RSID,
		genfile::Chromosome chromosome,
		uint32_t SNP_position,
		char a,
		char b
	) {
	std::istringstream inStream ;
	inStream.str(
		data::construct_snp_block(
			number_of_individuals,
			ID_field_length,
			SNPID,
			RSID,
			chromosome,
			SNP_position,
			a,
			b
		)
	) ;

	uint32_t number_of_individuals2 ;
	std::string SNPID2 ;
	std::string RSID2 ;
	genfile::Chromosome chromosome2 ;
	uint32_t SNP_position2 ;
	char a2 ;
	char b2 ;
	std::vector< probabilities > genotype_probabilities ;
	genotype_probabilities.resize( 1000 ) ;

	genfile::bgen::read_snp_block(
		inStream,
		make_setter( number_of_individuals2 ),
		make_setter( SNPID2 ),
		make_setter( RSID2 ),
		make_setter( chromosome2 ),
		make_setter( SNP_position2 ),
		make_setter( a2 ),
		make_setter( b2 ),
		ProbabilitySetter( genotype_probabilities )
	) ;

	TEST_ASSERT( number_of_individuals2 == number_of_individuals ) ;
	TEST_ASSERT( SNPID2 == SNPID ) ;
	TEST_ASSERT( RSID2 == RSID ) ;
	TEST_ASSERT( chromosome2 == chromosome ) ;
	TEST_ASSERT( SNP_position2 == SNP_position ) ;
	TEST_ASSERT( a2 == a ) ;
	TEST_ASSERT( b2 == b ) ;
	for( std::size_t i = 0; i < number_of_individuals; ++i ) {
		TEST_ASSERT( genotype_probabilities[i].AA == (i / 10000.0) ) ;
		TEST_ASSERT( genotype_probabilities[i].AB == (i / 10000.0) ) ;
		TEST_ASSERT( genotype_probabilities[i].BB == (i / 10000.0) ) ;
	}
}

double get_test_probability( std::size_t i ) { return static_cast< double> (i) / 10000.0 ; }

void do_snp_block_write_test( 
		uint32_t number_of_individuals,
		unsigned char ID_field_length,
		std::string SNPID,
		std::string RSID,
		genfile::Chromosome chromosome,
		uint32_t SNP_position,
		char a,
		char b
	) {
	std::ostringstream outStream ;
	
	genfile::bgen::write_snp_block( 
		outStream,
		number_of_individuals,
		ID_field_length,
		SNPID,
		RSID,
		chromosome,
		SNP_position,
		a,
		b,
		&get_test_probability,
		&get_test_probability,
		&get_test_probability
	) ;

	std::string expected = data::construct_snp_block(
		number_of_individuals,
		ID_field_length,
		SNPID,
		RSID,
		chromosome,
		SNP_position,
		a,
		b
		) ;
	
	// std::cout << "\"" << genfile::string_utils::to_hex( outStream.str() ) << "\"\n" ;
	// std::cout << "\"" << genfile::string_utils::to_hex( expected ) << "\"\n" ;

	TEST_ASSERT( outStream.str() == expected ) ;
}


AUTO_TEST_CASE( test_snp_block_input ) {
	std::cout << "test_snp_block_input\n" ;
	do_snp_block_read_test( 6, 6, "SNP 01", "RS 01", genfile::Chromosome1, 1000001, 'A', 'C' ) ;
	do_snp_block_read_test( 1000, 50, "01234567890123456789012345678901234567890123456789", "01234567890123456789012345678901234567890123456789", genfile::Chromosome22, 4294967295u, 'G', 'T' ) ;
}

AUTO_TEST_CASE( test_snp_block_output ) {
	std::cout << "test_snp_block_output\n" ;
	do_snp_block_write_test( 6, 6, "SNP 01", "RS 01", genfile::Chromosome1, 1000001, 'A', 'C' ) ;
	do_snp_block_write_test( 1000, 50, "01234567890123456789012345678901234567890123456789", "01234567890123456789012345678901234567890123456789", genfile::Chromosome22, 4294967295u, 'G', 'T' ) ;
}
