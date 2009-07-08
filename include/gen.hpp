#ifndef GEN_REFERENCE_IMPLEMENTATION_HPP
#define GEN_REFERENCE_IMPLEMENTATION_HPP

/* This file contains a reference implementation of the GEN file format
* described here:
* http://www.stats.ox.ac.uk/%7Emarchini/software/gwas/file_format.html
*/

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <stdint.h>
#include "snp_data_utils.hpp"

namespace genfile {
	namespace gen {
		namespace impl {
			typedef ::uint32_t uint32_t ;
		}

		typedef impl::uint32_t uint32_t ;
		
		// Get "header" information from the given stream, which must represent
		// a newly-opened GEN file.  This reports the number of snps and the number
		// of samples represented in the file.  This function consumes the
		// stream, which must be re-opened afterwards if further processing is required.
		template<
			typename NumberOfSNPBlocksSetter,
			typename NumberOfSamplesSetter,
			typename FlagsSetter
		>
		void read_header_information(
			std::istream& aStream,
			NumberOfSNPBlocksSetter set_number_of_snp_blocks,
			NumberOfSamplesSetter set_number_of_samples,
			FlagsSetter set_flags
		) ;

		// Read a SNP block from the (plain) gen file
		template<
			typename IntegerSetter,
			typename StringSetter,
			typename AlleleSetter,
			typename SNPPositionSetter,
			typename GenotypeProbabilitySetter
		>
		void read_snp_block(
			std::istream& in,
			IntegerSetter set_number_of_samples,
			StringSetter set_SNPID,
			StringSetter set_RSID,
			SNPPositionSetter set_SNP_position,
			AlleleSetter set_allele1,
			AlleleSetter set_allele2,
			GenotypeProbabilitySetter set_genotype_probabilities
		) ;

		/*
		* Function: write_snp_block()
		* Write a snp block with the given information to the given ostream object.
		* Genotype probabilities must be supplied by the given GenotypeProbabilityGetter
		* objects, which must be callable as
		* - get_AA_probability( index )
		* - get_AB_probability( index )
		* - get_BB_probability( index )
		* where index is the index of the individual in the SNP block.
		*/
		template< typename GenotypeProbabilityGetter >
		void write_snp_block(
			std::ostream& aStream,
			uint32_t number_of_samples,
			std::string SNPID,
			std::string RSID,
			uint32_t SNP_position,
			char first_allele,
			char second_allele,
			GenotypeProbabilityGetter get_AA_probability,
			GenotypeProbabilityGetter get_AB_probability,
			GenotypeProbabilityGetter get_BB_probability
		) ;



	/* IMPLEMENTATION */

		template<
			typename NumberOfSNPBlocksSetter,
			typename NumberOfSamplesSetter,
			typename FlagsSetter
		>
		void read_header_information(
			std::istream& aStream,
			NumberOfSNPBlocksSetter set_number_of_snp_blocks,
			NumberOfSamplesSetter set_number_of_samples,
			FlagsSetter set_flags
		) {
			uint32_t
				number_of_samples,
				number_of_snp_blocks ;
			std::string line ;

			read_snp_block( aStream, set_value( number_of_samples ), ignore(), ignore(), ignore(), ignore(), ignore(), ignore() ) ;
			for( number_of_snp_blocks = 1; std::getline( aStream, line ); ++number_of_snp_blocks ) ;

			// We should have reached eof.
			// If so, report the data now.
			if( aStream.eof() ) {
				set_number_of_snp_blocks( number_of_snp_blocks ) ;
				set_number_of_samples( number_of_samples ) ;
				set_flags( 0 ) ;
			}
		}
		

		template<
			typename IntegerSetter,
			typename StringSetter,
			typename AlleleSetter,
			typename SNPPositionSetter,
			typename GenotypeProbabilitySetter
		>
		void read_snp_block(
			std::istream& in,
			IntegerSetter set_number_of_samples,
			StringSetter set_SNPID,
			StringSetter set_RSID,
			SNPPositionSetter set_SNP_position,
			AlleleSetter set_allele1,
			AlleleSetter set_allele2,
			GenotypeProbabilitySetter set_genotype_probabilities
		) {
			// pull a line out of the stream and use a stringstream to read the line.
			std::string line ;
			std::getline( in, line ) ;
			std::istringstream inStream( line ) ;

			std::string SNPID, RSID ;
			double position ;
			char allele1, allele2 ;
			inStream >> SNPID >> RSID >> position >> allele1 >> allele2;
	
			if( inStream ) {
				// All good so far; take the plunge and write the data using the supplied setters.
				set_SNPID( SNPID ) ;
				set_RSID( RSID ) ;
				set_SNP_position( position ) ;
				set_allele1( allele1 ) ;
				set_allele2( allele2 ) ;

				// Now read genotype proportions.
				std::size_t number_of_samples = 0;
				std::size_t count = 0 ;
				std::vector< double > d(3) ;
				while( inStream >> d[count] ) {
					count = (count+1) % 3 ;
					if( count == 0 ) {
						set_genotype_probabilities( number_of_samples++, d[0], d[1], d[2] ) ;
					}
				}
				if( count % 3 != 0 ) {
					in.setstate( std::ios::failbit ) ;
				}
				else {
					set_number_of_samples( number_of_samples ) ;
				}
			}
		}
	}
}

#endif