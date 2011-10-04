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
#include "genfile/snp_data_utils.hpp"
#include "genfile/Error.hpp"
#include "genfile/get_set.hpp"

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

			void read_snp_identifying_data(
				std::istream& aStream,
				Chromosome* chromosome,
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
			int flags = 0 ;
			{
				std::string line ;
				std::getline( aStream, line ) ;
				// count spaces.
				std::string elt ;
				std::istringstream instr( line ) ;
				std::size_t count = 0 ;
				for( ; instr >> elt; ++count ) ;
				if(( count - 5 ) % 3 == 0 ) {
					// no chromosome column.
					number_of_samples = ( count - 5 ) / 3 ;
				}
				else if(( count - 6 ) % 3 == 0 ) {
					// chromosome column present.
					number_of_samples = ( count - 6 ) / 3 ;
					flags |= 0x1 ;
				}
				else {
					throw MalformedInputError( "(unknown)", 0 ) ;
				}
			}
			if( !aStream ) {
				throw MalformedInputError( "(unknown)", 0 ) ;
			}

			set_number_of_samples( number_of_samples ) ;
			set_flags( flags ) ;
		}

		uint32_t count_snp_blocks(
			std::istream& aStream
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
			genfile::Chromosome const& chromosome,
			std::string SNPID,
			std::string RSID,
			uint32_t SNP_position,
			char first_allele,
			char second_allele,
			GenotypeProbabilityGetter get_AA_probability,
			GenotypeProbabilityGetter get_AB_probability,
			GenotypeProbabilityGetter get_BB_probability
		) {
			aStream << chromosome << " " ;
			write_snp_block(
				aStream,
				number_of_samples,
				SNPID,
				RSID,
				SNP_position,
				first_allele,
				second_allele,
				get_AA_probability,
				get_AB_probability,
				get_BB_probability
			) ;
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