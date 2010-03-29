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
#include <cmath>
#include <algorithm>
#include <stdint.h>
#include "snp_data_utils.hpp"

namespace genfile {
	namespace gen {
		namespace impl {
			typedef ::uint32_t uint32_t ;
			
			template< typename FloatType >
			std::istream& read_float( std::istream& aStream, FloatType* number ) {
#ifndef GENFILE_USE_FAST_PARSE_METHODS
				return aStream >> (*number) ;
#else
				std::string float_str ;
				aStream >> float_str ;
				std::string::const_iterator i = float_str.begin() ,
					end_i = float_str.end() ;

				unsigned long integer_part = 0.0, fractional_part = 0ul;
				unsigned long* this_part = &integer_part ;
				int fractional_part_length = 0 ;
				bool good = true ;

				for( ; good && i != end_i; ++i ) {
					if( *i == '.' ) {
						this_part = &fractional_part ;
						fractional_part_length = std::distance( i, end_i ) - 1 ;
					}
					else if(( '0' <= *i ) && ( '9' >= *i ) ) {
						(*this_part) *= 10 ;
						(*this_part) += (*i - '0') ;
					}
					else {
						good = false ;
					}
				}

				if( good ) {
					*number = static_cast< FloatType >( integer_part ) + 
						+ ( static_cast< FloatType >( fractional_part ) / std::pow( 10.0, fractional_part_length )) ;
				}
				else {
					// Failed to parse.  Try to parse using usual method.
					std::istringstream inStream( float_str ) ;
					inStream >> (*number) ;
					if( inStream.fail() ) {
						aStream.setstate( std::ios::failbit ) ;
					}
				}
#endif
				return aStream ;
			}
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

		namespace impl {
	        void read_snp_identifying_data(
	            std::istream& aStream,
	            std::string* SNPID,
	            std::string* RSID,
	            uint32_t* SNP_position,
	            char* first_allele,
	            char* second_allele
			) ;
			
			template<
				typename IntegerSetter,
				typename GenotypeProbabilitySetter
			>
			void read_snp_probability_data(
				std::istream& inStream,
				IntegerSetter set_number_of_samples,
				GenotypeProbabilitySetter set_genotype_probabilities
			) {
				std::string line ;
				std::getline( inStream, line ) ;
				std::istringstream lineStream ;
				lineStream.str( line ) ;

				// Now read genotype proportions.
				std::size_t number_of_samples = 0;
				std::size_t count = 0 ;
				std::vector< double > d(3) ;
				while( impl::read_float(lineStream, &(d[count])) ) {
					count = (count+1) % 3 ;
					if( count == 0 ) {
						set_genotype_probabilities( number_of_samples++, d[0], d[1], d[2] ) ;
					}
				}
				if( count % 3 != 0 ) {
					inStream.setstate( std::ios::failbit ) ;
				}
				else {
					set_number_of_samples( number_of_samples ) ;
				}
			}
		}

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
			uint32_t number_of_samples ;
			read_snp_block( aStream, set_value( number_of_samples ), ignore(), ignore(), ignore(), ignore(), ignore(), ignore() ) ;
			if( !aStream ) {
				throw FileStructureInvalidError() ;
			}

			uint32_t number_of_snp_blocks = 1 ;
			std::vector< char > buffer( 10000000 ) ;
			do {
				aStream.read( &(buffer[0]), 10000000 ) ;
				number_of_snp_blocks += std::count( buffer.begin(), buffer.begin() + aStream.gcount(), '\n' ) ;
				// A gen file can't contain a blank line.
				// Because popular editors (vim, nano, ..., but not emacs) typically add a trailing newline,
				// we might get in the situation where the GEN file has two trailing newlines thus messing
				// with our count.
				// Therefore we check here for the special case where what we've read ends in two newlines.
				if( (aStream.gcount() > 1) && (buffer[ aStream.gcount() - 1] == '\n') && (buffer[ aStream.gcount() - 2] == '\n') ) {
					throw FileHasTwoConsecutiveNewlinesError() ;
				}
			}
			while( aStream ) ;

			// Most editors (vim, nano, but not emacs) automatically add a newline to the end of the file.
			// If the file has a trailing newline, we already have the correct count.
			// But if not, we've undercounted by one.
			if( aStream.gcount() > 0 ) {
				std::size_t pos = aStream.gcount() - 1 ;
				if( buffer[pos] != '\n' ) {
					++number_of_snp_blocks ;
				}
			}

			// We should have reached eof.
			// If so, report the data now.
			if( aStream.eof() ) {
				set_number_of_snp_blocks( number_of_snp_blocks ) ;
				set_number_of_samples( number_of_samples ) ;
				set_flags( 0 ) ;
			} else {
				throw FileStructureInvalidError() ;
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
			std::istream& inStream,
			IntegerSetter set_number_of_samples,
			StringSetter set_SNPID,
			StringSetter set_RSID,
			SNPPositionSetter set_SNP_position,
			AlleleSetter set_allele1,
			AlleleSetter set_allele2,
			GenotypeProbabilitySetter set_genotype_probabilities
		) {
			std::string SNPID, RSID ;
			uint32_t position ;
			char allele1, allele2 ;
			impl::read_snp_identifying_data( inStream, &SNPID, &RSID, &position, &allele1, &allele2 ) ;

			if( inStream ) {
				// All good so far; take the plunge and write the data using the supplied setters.
				set_SNPID( SNPID ) ;
				set_RSID( RSID ) ;
				set_SNP_position( position ) ;
				set_allele1( allele1 ) ;
				set_allele2( allele2 ) ;
				
				impl::read_snp_probability_data( inStream, set_number_of_samples, set_genotype_probabilities ) ;
			}
		}
		
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
		) {
			aStream
				<< SNPID << " "
				<< RSID << " "
				<< SNP_position << " "
				<< first_allele << " "
				<< second_allele << " " ;

			for( std::size_t i = 0 ; i < number_of_samples ; ++i ) {
				if( i > 0 ) {
					aStream << " " ;
				}
				aStream
					<< get_AA_probability(i) << " "
					<< get_AB_probability(i) << " "
					<< get_BB_probability(i) ;
			}

			aStream << "\n" ;
		}		
	}
}

#endif