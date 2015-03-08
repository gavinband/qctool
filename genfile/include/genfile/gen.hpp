
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

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
#include <boost/format.hpp>
#include "genfile/snp_data_utils.hpp"
#include "genfile/string_utils/slice.hpp"
#include "genfile/string_utils/strtod.hpp"
#include "genfile/Error.hpp"
#include "genfile/get_set.hpp"

#ifndef GENFILE_USE_FAST_PARSE_METHODS
#include "genfile/string_utils/string_utils.hpp"
#endif
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
				try {
					*number = string_utils::to_repr< FloatType >( float_str ) ;
				}
				catch( string_utils::StringConversionError const& e ) {
					aStream.setstate(  std::ios::failbit ) ;
				}
				return aStream ;
#endif	
			}
		}

		typedef impl::uint32_t uint32_t ;
		
		// Get "header" information from the given stream, which must represent
		// a newly-opened GEN file.  This reports the number of snps and the number
		// of samples represented in the file.  This function consumes the
		// stream, which must be re-opened afterwards if further processing is required.
		template<
			typename NumberOfSamplesSetter,
			typename FlagsSetter
		>
		void read_header_information(
			std::istream& aStream,
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
			std::string first_allele,
			std::string second_allele,
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
				std::string* first_allele,
				std::string* second_allele
			) ;

			void read_snp_identifying_data(
				std::istream& aStream,
				Chromosome* chromosome,
				std::string* SNPID,
				std::string* RSID,
				uint32_t* SNP_position,
				std::string* first_allele,
				std::string* second_allele
			) ;
			
			template<
				typename IntegerSetter,
				typename GenotypeProbabilitySetter
			>
			void read_snp_probability_data(
				std::istream& inStream,
				IntegerSetter set_number_of_samples,
				GenotypeProbabilitySetter set_genotype_probabilities,
				std::string& line
			) {
				using namespace string_utils ;
				//std::string line ;
				if( inStream.get() != ' ' ) {
					throw InputError(
						"inStream",
						"Expected a space before probabilities."
					) ;
				}
				std::getline( inStream, line ) ;
				std::vector< string_utils::slice > elts ;
				string_utils::slice( line ).split( " ", &elts ) ;
				if( elts.size() % 3 != 0 ) {
					throw InputError(
						"inStream",
						( boost::format( "Wrong number of elements (%d, not a multiple of 3)" ) % elts.size() ).str()
					) ;
				}
				std::size_t number_of_samples = 0 ;
				for( std::size_t i = 0; i < elts.size(); i += 3 ) {
					set_genotype_probabilities(
						number_of_samples++,
						strtod( elts[i] ),
						strtod( elts[i+1] ),
						strtod( elts[i+2] )
					) ;
				}
				set_number_of_samples( number_of_samples ) ;
			}
		}

		template<
			typename NumberOfSamplesSetter,
			typename FlagsSetter
		>
		void read_header_information(
			std::istream& aStream,
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
				double float_elt ;
				std::istringstream instr( line ) ;
				std::size_t count = 0 ;
				for( ; count < 6 && instr >> elt; ++count ) ;
				for( ; instr >> float_elt; ++count ) ;

				if( !instr ) {
					instr.peek() ;
					if( !instr.eof() ) {
						throw MalformedInputError( "(unknown)", 0 ) ;
					}
				}
				if( count >= 5 && ( count - 5 ) % 3 == 0 ) {
					// no chromosome column.
					number_of_samples = ( count - 5 ) / 3 ;
				}
				else if( count >= 6 && ( count - 6 ) % 3 == 0 ) {
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
			std::string first_allele,
			std::string second_allele,
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
			std::string first_allele,
			std::string second_allele,
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
