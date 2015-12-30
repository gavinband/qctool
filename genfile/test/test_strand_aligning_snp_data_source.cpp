
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
#include "genfile/GenFileSNPDataSource.hpp"
#include "genfile/GenFileSNPDataSink.hpp"
#include "genfile/StrandAligningSNPDataSource.hpp"
#include "genfile/SNPDataSink.hpp"
#include "stdint.h"

// #define DEBUG 1

AUTO_TEST_SUITE( StrandAligningSNPDataSourceTest )

// The following section contains a simple snp block writer.
namespace data {
	namespace {
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
}

namespace {
	std::string complement( std::string allele ) {
		assert( allele.size() == 1 ) ;
		switch( allele[0] ) {
			case 'A': return "T" ; break ;
			case 'T': return "A" ; break ;
			case 'C': return "G" ; break ;
			case 'G': return "C" ; break ;
			case '?': return "?" ; break ;
			default: assert(0) ; break ;
		}
	}
}

AUTO_TEST_CASE( test_strand_aligning_snp_data_source ) {
	std::cerr << "test_strand_aligning_snp_data_source()..." ;
	
	std::size_t const N = 100 ; // do 100 iterations.
	std::string const question_mark = "?" ;
	
	std::vector< genfile::VariantIdentifyingData > snps ;
	{
		std::auto_ptr< std::istream > istr( new std::istringstream( data::data ) ) ;
		std::auto_ptr< genfile::SNPDataSource > plain_source(
			new genfile::GenFileSNPDataSource( istr, genfile::Chromosome( "01" ) )
		) ;
		genfile::VariantIdentifyingData snp ;
		while( plain_source->get_snp_identifying_data( &snp )) {
			snps.push_back( snp ) ;
			plain_source->ignore_snp_probability_data() ;
		}
		assert( snps.size() == data::number_of_snps ) ;
	}
	
	for( std::size_t i = 0; i < N; ++i ) {
		genfile::StrandAligningSNPDataSource::StrandAlignments strand_alignments ;
		
		std::vector< std::size_t > aligned_unaligned_map ;
		{
#if DEBUG
			std::cerr << i << " " ;
#endif
			for( std::size_t j = 0; j < snps.size(); ++j ) {
				switch( (i*j) % 3 ) {
					case 0:
						strand_alignments[ snps[j] ].strand = '+' ;
						strand_alignments[ snps[j] ].flip = '+' ;
						aligned_unaligned_map.push_back( j ) ;
					break ;
					case 1:
						strand_alignments[ snps[j] ].strand = '-' ;
						strand_alignments[ snps[j] ].flip = '+' ;
						aligned_unaligned_map.push_back( j ) ;
					break ;
					case 2:
						strand_alignments[ snps[j] ].strand = '?' ;
						strand_alignments[ snps[j] ].flip = '+' ;
					break ;
				}
			}
			assert( strand_alignments.size() == data::number_of_snps ) ;
		}
		
		genfile::VariantIdentifyingData snp ;

		// Set up and read from a plain, unaligned source.
		std::auto_ptr< genfile::SNPDataSource > plain_source ;
		{
			std::auto_ptr< std::istream > istr( new std::istringstream( data::data ) ) ;
			plain_source.reset( new genfile::GenFileSNPDataSource( istr, genfile::Chromosome( "01" ) ) );
		}
		std::vector< genfile::VariantIdentifyingData > snps1 ;
		std::vector< genfile::SingleSNPGenotypeProbabilities > probabilities1 ;
		{
			while( plain_source->get_snp_identifying_data( &snp )) {
				snps1.push_back( snp ) ;
				genfile::SingleSNPGenotypeProbabilities probs( data::number_of_samples ) ;
				plain_source->read_snp_probability_data(
					genfile::set_genotypes( probs )
				) ;
				probabilities1.push_back( probs ) ;
			}
			TEST_ASSERT( plain_source->number_of_snps_read() == data::number_of_snps ) ;
			TEST_ASSERT( snps1.size() == probabilities1.size() ) ;
			TEST_ASSERT( snps1.size() == data::number_of_snps ) ;
		}

		// Set up and read from an aligned source.
		std::auto_ptr< genfile::SNPDataSource > aligned_source ;
		std::vector< genfile::VariantIdentifyingData > snps2 ;
		std::vector< genfile::SingleSNPGenotypeProbabilities > probabilities2 ;
		{
			std::auto_ptr< std::istream > istr( new std::istringstream( data::data ) ) ;
			aligned_source.reset( new genfile::GenFileSNPDataSource( istr, genfile::Chromosome( "01" ) ) );
			aligned_source.reset( new genfile::StrandAligningSNPDataSource( aligned_source, strand_alignments )) ;
		}
		
		{
			for( std::size_t count = 0; aligned_source->get_snp_identifying_data( &snp ); ++count ) {
				snps2.push_back( snp ) ;
				genfile::SingleSNPGenotypeProbabilities probs( data::number_of_samples ) ;
				aligned_source->read_snp_probability_data(
					genfile::set_genotypes( probs )
				) ;
				probabilities2.push_back( probs ) ;
			}
#if DEBUG
			std::cerr << "\n" << aligned_source->number_of_snps_read() << " : " << aligned_unaligned_map.size() << "\n" ;
#endif
			TEST_ASSERT( aligned_source->number_of_snps_read() == aligned_unaligned_map.size() ) ;
			TEST_ASSERT( snps2.size() == probabilities2.size() ) ;
			TEST_ASSERT( snps2.size() == aligned_unaligned_map.size() ) ;
		}

		for( std::size_t aligned_i = 0; aligned_i < aligned_unaligned_map.size(); ++aligned_i ) {
			TEST_ASSERT( probabilities2[aligned_i].size() == data::number_of_samples ) ;
			std::size_t i = aligned_unaligned_map[ aligned_i ] ;

			for( std::size_t j = 0; j < data::number_of_samples; ++j ) {
				for( std::size_t g = 0; g < 3; ++g ) {
#if DEBUG
					std::cerr << i << " " << j << " " << g << ": " << probabilities1[i]( j, g ) << " : " << probabilities2[aligned_i]( j, g ) << "\n" ;
#endif
					TEST_ASSERT( probabilities1[i]( j, g ) == probabilities2[aligned_i]( j, g ) ) ;
				}
			}

			genfile::VariantIdentifyingData expected_aligned_snp = snps1[i] ;

			switch( strand_alignments[ snps1[i] ].strand ) {
				case '+':
					break ;
				case '-':
					expected_aligned_snp.set_allele( 0, complement( expected_aligned_snp.get_allele( 0 ) ) ) ;
					expected_aligned_snp.set_allele( 1, complement( expected_aligned_snp.get_allele( 1 ) ) ) ;
					break ;
				case '?':
					expected_aligned_snp.set_allele( 0, question_mark ) ;
					expected_aligned_snp.set_allele( 1, question_mark ) ;
					break ;
				default:
					assert(0) ;
					break ;
			}
			if( expected_aligned_snp != snps2[aligned_i] ) {
				std::cerr << "Oh dear, got " << snps2[aligned_i] << ", expected " << expected_aligned_snp << ".\n" ;
			}
			TEST_ASSERT( expected_aligned_snp == snps2[aligned_i] ) ;
		}
	}
	
	std::cerr << "ok.\n" ;
}

AUTO_TEST_SUITE_END()
