
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cassert>
#include <fstream>
#include "test_case.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/GenFileSNPDataSource.hpp"
#include "genfile/BGenFileSNPDataSource.hpp"
#include "genfile/SNPDataSink.hpp"
#include "genfile/GenFileSNPDataSink.hpp"
#include "genfile/BGenFileSNPDataSink.hpp"
#include "genfile/bgen/bgen.hpp"
#include "stdint.h"

// #define DEBUG 1

AUTO_TEST_SUITE( test_genfile_io )

// The following section contains a simple snp block writer.
namespace data {
	unsigned int const number_of_samples = 7 ;
	unsigned int const number_of_snps = 16 ;
	std::string data =
		"--- rs11089130 14431347 C G 0.11 0.44 0.44 0.11 0.44 0.44 0.11 0.44 0.44 0.11 0.44 0.44 0.11 0.44 0.44 0.11 0.44 0.44 0.11 0.44 0.44\n"
		"--- rs738829 14432618 A G 0.04 0.31 0.65 0.04 0.31 0.65 0.04 0.31 0.65 0.04 0.31 0.65 0.04 0.31 0.65 0.04 0.31 0.65 0.04 0.31 0.65\n"
		"--- rs915674 14433624 A G 0.02 0.23 0.75 0.02 0.23 0.75 0.02 0.23 0.75 0.02 0.23 0.75 0.02 0.23 0.75 0.02 0.23 0.75 0.02 0.23 0.75\n"
		"--- rs915675 14433659 A C 0.04 0.33 0.63 0.04 0.33 0.63 0.04 0.33 0.63 0.04 0.33 0.63 0.04 0.33 0.63 0.04 0.33 0.63 0.04 0.33 0.63\n"
		"--- rs915677 14433758 A G 0 0.13 0.87 0 0.13 0.87 0 0.13 0.87 0 0.13 0.87 0 0.13 0.87 0 0.13 0.87 0 0.13 0.87\n"
		"--- rs9604721 14434713 C T 0.93 0.07 0 0.93 0.07 0 0.93 0.07 0 0.93 0.07 0 0.93 0.07 0 0.93 0.07 0 0.93 0.07 0\n"
		"--- rs4389403 14435070 A G 0.01 0.19 0.79 0.01 0.19 0.79 0.01 0.19 0.79 0.01 0.19 0.79 0.01 0.19 0.79 0.01 0.19 0.79 0.01 0.19 0.79\n"
		"--- rs5746356 14439734 C T 0.55 0.38 0.07 0.55 0.38 0.07 0.55 0.38 0.07 0.55 0.38 0.07 0.55 0.38 0.07 0.55 0.38 0.07 0.55 0.38 0.07\n"
		"--- rs9617528 14441016 C T 0.07 0.38 0.55 0.07 0.38 0.55 0.07 0.38 0.55 0.07 0.38 0.55 0.07 0.38 0.55 0.07 0.38 0.55 0.07 0.38 0.55\n"
		"--- rs2154787 14449374 C T 0.47 0.43 0.10 0.47 0.43 0.10 0.47 0.43 0.10 0.47 0.43 0.10 0.47 0.43 0.10 0.47 0.43 0.10 0.47 0.43 0.10\n"
		"--- rs12484041 14452292 A G 0 0.03 0.97 0 0.03 0.97 0 0.03 0.97 0 0.03 0.97 0 0.03 0.97 0 0.03 0.97 0 0.03 0.97\n"
		"--- rs10154731 14479437 A G 0 0.15 0.84 0 0.15 0.84 0 0.15 0.84 0 0.15 0.84 0 0.15 0.84 0 0.15 0.84 0 0.15 0.84\n"
		"--- rs11913813 14480059 C G 0.92 0.08 0 0.92 0.08 0 0.92 0.08 0 0.92 0.08 0 0.92 0.08 0 0.92 0.08 0 0.92 0.08 0\n"
		"--- rs2260460 14482325 C T 0.90 0.10 0 0.90 0.10 0 0.90 0.10 0 0.90 0.10 0 0.90 0.10 0 0.90 0.10 0 0.90 0.10 0\n"
		"--- rs1964966 14483902 A G 0.74 0.24 0.02 0.74 0.24 0.02 0.74 0.24 0.02 0.74 0.24 0.02 0.74 0.24 0.02 0.74 0.24 0.02 0.74 0.24 0.02\n"
		"SNP_A-1928576 rs11705026 14490036 G T 0 1 0 0.01 0.99 0 0 1 0 0 0 1 0.06 0.94 0 0 1 0 0 0 1\n" ;
}

// The following section defines the needed objects for use with the bgen.hpp implementation.
template< typename T >
struct Setter
{
	Setter( T& field ): m_field( field ) {} ;
	void operator()( T const& value ) { m_field = value ; }
private:
	T& m_field ;
} ;

template< typename T >
Setter< T > make_setter( T& field ) { return Setter<T>( field ) ; }


struct probabilities_t {
	double AA, AB, BB ;
	bool operator==( probabilities_t const& other ) const {
		return AA == other.AA && AB == other.AB && BB == other.BB ;
	}
} ;

std::ostream& operator<<( std::ostream& o, probabilities_t p ) {
	return o << p.AA << ", " << p.AB << ", " << p.BB ;
}

struct ProbabilitySetter {
	ProbabilitySetter( std::vector< probabilities_t >& probabilities ): m_probabilities( probabilities ) {}
	void operator() ( std::size_t i, double aa, double ab, double bb ) {
		m_probabilities[i].AA = aa ;
		m_probabilities[i].AB = ab ;
		m_probabilities[i].BB = bb ;
	}

private:
	std::vector<probabilities_t>& m_probabilities ;
} ;

struct ProbabilityGetter {
	ProbabilityGetter( std::vector< probabilities_t > const& probabilities, int index ): m_probabilities( probabilities ), m_index( index ) {}
	double operator() ( std::size_t i ) const {
		switch( m_index ) {
			case 0: return m_probabilities[i].AA ; break ;
			case 1: return m_probabilities[i].AB ; break ;
			case 2: return m_probabilities[i].BB ; break ;
			default:
			assert(0);
			break ;
		}
		return 0.0 ;
	}
private:
	std::vector<probabilities_t> const& m_probabilities ;
	int m_index ;
} ;

struct SnpData {
	
	SnpData(): probabilities( data::number_of_samples ) {} ;
	
	genfile::VariantIdentifyingData snp ;
	std::vector< probabilities_t > probabilities ;
	
	bool operator==( SnpData const& other ) const {
		return snp == other.snp && probabilities == other.probabilities ;
	}
} ;

namespace {
	std::string get_test_sample_name( std::size_t i ) {
		return "test_sample_" + genfile::string_utils::to_string( i ) ;
	}
}

void copy_gen_file( genfile::SNPDataSource& snp_data_source, genfile::SNPDataSink& snp_data_sink ) {
	SnpData snp_data ;
	snp_data_sink.set_sample_names( snp_data_source.number_of_samples(), &get_test_sample_name ) ;
	while( snp_data_source.get_snp_identifying_data( &snp_data.snp ) ) {
		snp_data_source.read_snp_probability_data( ProbabilitySetter( snp_data.probabilities ) ) ;
		snp_data_sink.write_snp(
			snp_data_source.number_of_samples(),
			snp_data.snp.get_alternate_identifiers_as_string(),
			snp_data.snp.get_rsid(),
			snp_data.snp.get_position().chromosome(),
			snp_data.snp.get_position().position(),
			snp_data.snp.get_allele(0),
			snp_data.snp.get_allele(1),
			ProbabilityGetter( snp_data.probabilities, 0 ),
			ProbabilityGetter( snp_data.probabilities, 1 ),
			ProbabilityGetter( snp_data.probabilities, 2 )
		) ;
	}
}

void copy_gen_file( std::string original, genfile::SNPDataSink& target ) {
	genfile::GenFileSNPDataSource gen_file_snp_data_source( original, genfile::UnidentifiedChromosome ) ;
	copy_gen_file( gen_file_snp_data_source, target ) ;
}

void copy_gen_file( std::string const& original, std::string const& target ) {
	genfile::SNPDataSink::UniquePtr snp_data_sink( genfile::SNPDataSink::create( target )) ;
	copy_gen_file( original, *snp_data_sink ) ;
}

// The following section contains the main tests.
void create_files( std::string original, std::string gen, std::string bgen_v11, std::string bgen_v12, std::string zgen ) {
	// set up our original file.
	{
#if DEBUG
		std::cerr << "Creating original...\n" ;
#endif
		std::ofstream original_file( original.c_str() ) ;
		original_file << data::data ;
	}

	// In the first version of this function,
	// we set up the data sinks by hand.  In create_files2 below we'll use the
	// SNPDataSink::create() factory instead.
	
	// Construct plain gen file.
	{
#if DEBUG
		std::cerr << "Creating gen...\n" ;
#endif
		genfile::GenFileSNPDataSink gen_file_snp_data_sink( gen, genfile::UnidentifiedChromosome ) ;
		copy_gen_file( original, gen_file_snp_data_sink ) ;
	}

	// Construct bgen v1.1
	genfile::SNPDataSink::Metadata metadata ;
	{
#if DEBUG
		std::cerr << "Creating v11 bgen...\n" ;
#endif
		genfile::BGenFileSNPDataSink bgen_file_snp_data_sink( bgen_v11, metadata, "v11" ) ;
		copy_gen_file( original, bgen_file_snp_data_sink ) ;
	}

	// Construct bgen v1.2
	{
#if DEBUG
		std::cerr << "Creating bgen v1.1...\n" ;
#endif
		genfile::BGenFileSNPDataSink bgen_file_snp_data_sink( bgen_v12, metadata, "v11" ) ;
		copy_gen_file( original, bgen_file_snp_data_sink ) ;
	}

	// Construct zipped gen file
	{
#if DEBUG
		std::cerr << "Creating gen.gz...\n" ;
#endif
		genfile::GenFileSNPDataSink zipped_gen_file_snp_data_sink( zgen, genfile::UnidentifiedChromosome, "gzip_compression" ) ;
		copy_gen_file( original, zipped_gen_file_snp_data_sink ) ;
	}
}

// The following section contains the main tests.
void create_files2( std::string original, std::string gen, std::string bgen_v11, std::string bgen_v12, std::string zgen ) {
	std::cerr << "Creating gen...\n" ;
	copy_gen_file( original, gen ) ;
	std::cerr << "Creating bgen v11...\n" ;
	copy_gen_file( original, bgen_v11 ) ;
	std::cerr << "Creating bgen v12...\n" ;
	copy_gen_file( original, bgen_v12 ) ;
	std::cerr << "Creating gen.gz...\n" ;
	copy_gen_file( original, zgen ) ;
}

std::vector< SnpData > read_gen_file( genfile::SNPDataSource& snp_data_source ) {
	std::vector< SnpData > result ;
	SnpData snp_data ;
	while( snp_data_source.get_snp_identifying_data( &snp_data.snp ) ) {
		snp_data_source.read_snp_probability_data( ProbabilitySetter( snp_data.probabilities ) ) ;
		result.push_back( snp_data ) ;
	}
	return result ;
}

std::vector< SnpData > read_gen_file( std::string filename ) {
	std::auto_ptr< genfile::SNPDataSource > snp_data_source_ptr( genfile::SNPDataSource::create( filename, genfile::Chromosome( "NA" ) )) ;
	return read_gen_file( *snp_data_source_ptr ) ;
}

AUTO_TEST_CASE( test_formats ) {
	std::string original = tmpnam(0) + std::string( ".gen" );
	std::string gen = tmpnam(0) + std::string( ".gen" );
	std::string gen2 = tmpnam(0) + std::string( ".gen" );
	std::string bgen_v11 = tmpnam(0) + std::string( ".v11.bgen" );
	std::string bgen_v112 = tmpnam(0) + std::string( ".v11.bgen" );
	std::string bgen_v12 = tmpnam(0) + std::string( ".v12.bgen" );
	std::string bgen_v122 = tmpnam(0) + std::string( ".v12.bgen" );
	std::string zgen = tmpnam(0) + std::string( ".gen.gz" );
	std::string zgen2 = tmpnam(0) + std::string( ".gen.gz" );

#if DEBUG
	std::cerr << "Creating files...\n" ;
	std::cerr << "Master = \"" << original << "\".\n" ;
#endif
	
	create_files( original, gen, bgen_v11, bgen_v12, zgen ) ;
	create_files2( original, gen2, bgen_v112, bgen_v122, zgen2 ) ;
	
	genfile::GenFileSNPDataSource original_file_snp_data_source( original, genfile::UnidentifiedChromosome ) ;
	genfile::GenFileSNPDataSource gen_file_snp_data_source( gen, genfile::UnidentifiedChromosome ) ;
	genfile::BGenFileSNPDataSource bgen_v11_file_snp_data_source( bgen_v11 ) ;
	genfile::BGenFileSNPDataSource bgen_v12_file_snp_data_source( bgen_v12 ) ;
	genfile::BGenFileSNPDataSource bgen_v11_file_snp_data_source2( bgen_v112 ) ;
	genfile::BGenFileSNPDataSource bgen_v12_file_snp_data_source2( bgen_v122 ) ;
	genfile::GenFileSNPDataSource zgen_file_snp_data_source( zgen, genfile::UnidentifiedChromosome ) ;
	genfile::GenFileSNPDataSource gen_file_snp_data_source2( gen2, genfile::UnidentifiedChromosome ) ;
	genfile::GenFileSNPDataSource zgen_file_snp_data_source2( zgen, genfile::UnidentifiedChromosome ) ;

	std::vector< std::vector< SnpData > > results ;

	results.push_back(read_gen_file( original_file_snp_data_source )) ;
	
	results.push_back( read_gen_file( gen_file_snp_data_source )) ;
	results.push_back( read_gen_file( bgen_v11_file_snp_data_source )) ;
	results.push_back( read_gen_file( bgen_v12_file_snp_data_source )) ;
	results.push_back( read_gen_file( zgen_file_snp_data_source )) ;

	results.push_back( read_gen_file( gen )) ;
	results.push_back( read_gen_file( bgen_v11 )) ;
	results.push_back( read_gen_file( bgen_v12 )) ;
	results.push_back( read_gen_file( zgen )) ;

	results.push_back( read_gen_file( gen_file_snp_data_source2 )) ;
	results.push_back( read_gen_file( bgen_v11_file_snp_data_source2 )) ;
	results.push_back( read_gen_file( bgen_v12_file_snp_data_source2 )) ;
	results.push_back( read_gen_file( zgen_file_snp_data_source2 )) ;
	
	results.push_back( read_gen_file( gen2 )) ;
	results.push_back( read_gen_file( bgen_v112 )) ;
	results.push_back( read_gen_file( bgen_v122 )) ;
	results.push_back( read_gen_file( zgen2 )) ;

	TEST_ASSERT( gen_file_snp_data_source.number_of_samples() == data::number_of_samples ) ;
	TEST_ASSERT( bgen_v11_file_snp_data_source.number_of_samples() == data::number_of_samples ) ;
	TEST_ASSERT( bgen_v12_file_snp_data_source.number_of_samples() == data::number_of_samples ) ;
	TEST_ASSERT( zgen_file_snp_data_source.number_of_samples() == data::number_of_samples ) ;

	TEST_ASSERT( gen_file_snp_data_source2.number_of_samples() == data::number_of_samples ) ;
	TEST_ASSERT( bgen_v11_file_snp_data_source2.number_of_samples() == data::number_of_samples ) ;
	TEST_ASSERT( bgen_v12_file_snp_data_source2.number_of_samples() == data::number_of_samples ) ;
	TEST_ASSERT( zgen_file_snp_data_source2.number_of_samples() == data::number_of_samples ) ;

	for( std::size_t i = 0; i < results.size(); ++i ) {
#if DEBUG
		std::cerr << "Testing :" << i << "\n" ;
#endif
		TEST_ASSERT( results[i].size() == data::number_of_snps ) ;
			for( std::size_t j = 0; j < results[i].size(); ++j ) {
				TEST_ASSERT( results[i][j].snp.get_alternate_identifiers_as_string() == results[0][j].snp.get_alternate_identifiers_as_string() ) ;
				TEST_ASSERT( results[i][j].snp.get_rsid() == results[0][j].snp.get_rsid() ) ;
				TEST_ASSERT( results[i][j].snp.get_position() == results[0][j].snp.get_position() ) ;
				TEST_ASSERT( results[i][j].snp.get_allele(0) == results[0][j].snp.get_allele(0) ) ;
				TEST_ASSERT( results[i][j].snp.get_allele(1) == results[0][j].snp.get_allele(1) ) ;
				if( results[i][j].probabilities != results[0][j].probabilities ) {
					for( std::size_t k = 0; k < results[i][j].probabilities.size(); ++k ) {
#if DEBUG
						std::cerr << results[0][j].probabilities[k] << " : " << results[i][j].probabilities[k] << ".\n" ;
#endif
						// BGEN with 16 bits gets accuracy to 1/32768 which is ~4dps of accuracy
						BOOST_CHECK_SMALL( results[0][j].probabilities[k].AA - results[i][j].probabilities[k].AA, 0.00005 ) ;
						BOOST_CHECK_SMALL( results[0][j].probabilities[k].AB - results[i][j].probabilities[k].AB, 0.00005 ) ;
						BOOST_CHECK_SMALL( results[0][j].probabilities[k].BB - results[i][j].probabilities[k].BB, 0.00005 ) ;
					}
				}
			}
		//TEST_ASSERT( results[i] == results[0] ) ;
	}
}

AUTO_TEST_SUITE_END()
