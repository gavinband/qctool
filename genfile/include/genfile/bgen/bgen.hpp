
//			Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//	  (See accompanying file LICENSE_1_0.txt or copy at
//			http://www.boost.org/LICENSE_1_0.txt)

#ifndef BGEN_REFERENCE_IMPLEMENTATION_HPP
#define BGEN_REFERENCE_IMPLEMENTATION_HPP

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <stdint.h>
//#include "genfile/snp_data_utils.hpp"
//#include "genfile/get_set.hpp"
#include "genfile/zlib.hpp"
#include "genfile/types.hpp"
#include "genfile/MissingValue.hpp"

/*
* This file contains a reference implementation of the BGEN file format
* specification described at:
* http://www.well.ox.ac.uk/~gav/bgen_format/bgen_format.html
*
*/

// #define DEBUG_BGEN_FORMAT 1

#if DEBUG_BGEN_FORMAT
#include <iostream>
#include <iomanip>
#endif

///////////////////////////////////////////////////////////////////////////////////////////
// INTERFACE
///////////////////////////////////////////////////////////////////////////////////////////

namespace genfile {
	namespace bgen {
#if DEBUG_BGEN_FORMAT
		namespace impl {
			std::string to_hex( std::string const& str ) ;

			template< typename I, typename I2 >
			std::string to_hex( I i, I2 end_i ) {
				std::ostringstream o ;
				for( std::size_t count = 0; i < end_i; ++i, ++count ) {
					if( count % 4 == 0 )
						o << "|" ;
					o << std::hex << std::setw(2) << std::setfill('0') << static_cast<int> ( static_cast<unsigned char>( *i ) ) ;
				}
				return o.str() ;
			}
		}
#endif

		// class thrown when errors are detected
		struct BGenError: public virtual std::exception {
			~BGenError() throw() {}
			char const* what() const throw() { return "BGenError" ; }
		} ;

		// integer types
		typedef ::uint32_t uint32_t ;
		typedef ::uint16_t uint16_t ;

		// Header flag definitions
		enum FlagMask { e_NoFlags = 0, e_CompressedSNPBlocks = 0x1, e_Layout = 0x3C } ;
		enum Layout { e_v10Layout = 0x0, e_v11Layout = 0x4, e_v12Layout = 0x8 } ;
		enum Structure { e_SampleIdentifiers = 0x80000000 } ;
		
		// Structure containing information from the header block.
		struct Context {
			Context() ;
			Context( Context const& other ) ;
			Context& operator=( Context const& other ) ;
			uint32_t header_size() const ;
		public:	
			uint32_t number_of_samples ;
			uint32_t number_of_variants ;
			std::string magic ;
			std::string free_data ;
			uint32_t flags ;
		} ;

		// Read the offset from the start of the stream.
		void read_offset( std::istream& iStream, uint32_t* offset ) ;
		// Write an offset value to the stream.
		void write_offset( std::ostream& oStream, uint32_t const offset ) ;

		// Read a header block from the supplied stream,
		// filling the fields of the supplied context object.
		std::size_t read_header_block(
			std::istream& aStream,
			Context* context
		) ;

		// Write a bgen header block to the supplied stream,
		// taking data from the fields of the supplied context object.
		void write_header_block(
			std::ostream& aStream,
			Context const& context
		) ;

		// Read a sample identifier block from the given stream.
		// The setter object passed in must be a unary function or
		// function object that takes a string.  It must be callable as
		// setter( value ) ;
		// where value is of type std::string.  It will be called once for
		// each sample identifier in the block, in the order they occur.
		template< typename SampleSetter >
		std::size_t read_sample_identifier_block(
			std::istream& aStream,
			Context const& context,
			SampleSetter setter
		) ;

		// Write the sample identifiers contained in the
		// given vector to the stream.
		std::size_t write_sample_identifier_block(
			std::ostream& aStream,
			Context const& context,
			std::vector< std::string > const& sample_ids
		) ;

		// Attempt to read identifying information for the next variant in the file.
		// This function will return true if SNP data was successfully read or false if the initial
		// reads met an EOF.  It will throw an BGenError if only a subset of fields can be read before EOF.
		// The two setter objects must be callable as
		// set_number_of_alleles( n )
		// set_allele( i, a )
		// where n is an unsigned integer representing the number of alleles,
		// i is an unsigned integer in the range 0..n-1, and a is of type std::string
		// representing the ith allele.
		template<
			typename NumberOfAllelesSetter,
			typename AlleleSetter
		>
		bool read_snp_identifying_data(
			std::istream& aStream,
			Context const& context,
			std::string* SNPID,
			std::string* RSID,
			std::string* chromosome,
			uint32_t* SNP_position,
			NumberOfAllelesSetter set_number_of_alleles,
			AlleleSetter set_allele
		) ;

		// Read identifying data fields for the next variant in the file, assuming 2 alleles.
		// This function forwards to the generic multi-allele version, above, and will throw
		// a BGenError if the number of alleles is different than 2.
		bool read_snp_identifying_data(
			std::istream& aStream,
			Context const& context,
			std::string* SNPID,
			std::string* RSID,
			std::string* chromosome,
			uint32_t* SNP_position,
			std::string* first_allele,
			std::string* second_allele
		) ;
			
		// Write identifying data fields for the given variant.
		void write_snp_identifying_data(
			std::ostream& aStream,
			Context const& context,
			unsigned char max_id_size,
			std::string SNPID,
			std::string RSID,
			std::string chromosome,
			uint32_t SNP_position,
			std::string first_allele,
			std::string second_allele
		) ;

		// Ignore (and seek forward past) the genotype data block contained in the given stream.
		// (The main purpose of this function is to encapsulate the slightly complicated rules
		// for reading or computing the compressed data size needed to skip the right number of bytes).
		void ignore_genotype_data_block(
			std::istream& aStream,
			Context const& context
		) ;

		// Low-level function which reads raw probability data from a genotype data block
		// contained in the input stream into a supplied buffer. The buffer will be resized
		// to fit the data (incurring an allocation if the buffer is not large enough.)
		// (The main purpose of this function is to encapsulate the slightly complicated rules
		// for reading or computing the compressed data size.  Where applicable this function first
		// reads the four bytes indicating the compressed data size; these are discarded
		// and do not appear in the buffer).
		void read_genotype_data_block(
			std::istream& aStream,
			Context const& context,
			std::vector< char >* buffer1
		) ;

		// Low-level function which uncompresses probability data stored in the genotype data block
		// contained in the first buffer into a second buffer (or just copies it over if the probability
		// data is not compressed.) The second buffer will be resized to fit the result (incurring an
		// allocation if the buffer is not large enough).
		// (The main purpose of this function is to encapsulate the rules for reading or computing
		// the uncompressed data size.  Where applicable this function first consumes four bytes of
		// the input buffer and interprets them as the uncompressed data size, before uncompressing
		// the rest.)
		// Usually bgen files are stored compressed.  If the data is not compressed, this function
		// simply copies the source buffer to the target buffer.
		void uncompress_probability_data(
			Context const& context,
			std::vector< char > const& buffer1,
			std::vector< char >* buffer2
		) ;

		// template< typename Setter >
		// parse uncompressed genotype probability data stored in the given buffer.
		// Values are returned as doubles or as missing values using the
		// setter object provided.  The setter must support the following expressions:
		//
		// - setter.initialise( N, K ) ;
		// where N and K are convertible from std::size_t and represent the number of samples and number of alleles
		// present for the variant.
		//
		// - setter.set_sample( i )
		// where i is convertible from std::size_t and represents the index of a sample between 0 and N-1.
		// This function should return a value convertible to bool, representing a hint as to whether the setter
		// wants to use the values for this sample or not.
		// If the return value is false, the implementation may choose not to report any data for this sample,
		// whence the next call, if any, will be set_sample(i+1).
		//
		// - setter.set_number_of_entries( Z, order_type, value_type ) ;
		// where Z is an unsigned integer reflecting the number of probability values comprising the data
		// for this sample; order_type is of type genfile::OrderType and value_type is of type genfile::ValueType.
		// This is called once per sample (for which set_sample is true)
		// and informs the setter the number of probability values present for this sample (depending on the
		// ploidy, the number of alleles, and whether the data is phased.)  E.g. In the common case of a diploid sample
		// at a biallelic variant, this will be called with Z=3.
		// The order_type is ePerUnorderedGenotype or ePerPhasedHaplotypePerAllele depending on whether or
		// not the data is phased.  For bgen <= 1.2, all data is interpreted as probabilities so the value_type
		// is always eProbability.
		//
		// - setter( value ) ;
		// where value is of type double or genfile::MissingValue.
		// This is called Z times and reflects the probability data stored for this sample.
		// Samples with missing data have Z missing values; those without missing data have
		// Z double values.
		//
		template< typename Setter >
		void parse_probability_data(
			char const* buffer,
			char const* const end,
			Context const& context,
			Setter& setter
		) ;
		
		// Utility function which wraps the above steps for reading probability data into a single function.
		// Concretely this function:
		// 1: calls read_genotype_data_block(), reading probability data from the given stream.
		// 2: calls uncompress_probability_data() to uncompress the data where necessary.
		// 3: calls parse_probability_data() to parse it, returning values using the setter object provided.
		// The buffers are used as intermediate storage and will be resized to fit data as needed.
		template< typename Setter >
		void read_and_parse_genotype_data_block(
			std::istream& aStream,
			Context const& context,
			Setter& setter,
			std::vector< char >* buffer1,
			std::vector< char >* buffer2
		) ;
	}
}	

///////////////////////////////////////////////////////////////////////////////////////////
// IMPLEMENTATION
///////////////////////////////////////////////////////////////////////////////////////////

namespace genfile {
	namespace bgen {
		// Read an integer stored in little-endian format into an integer stored in memory.
		template< typename IntegerType >
		char const* read_little_endian_integer( char const* buffer, char const* const end, IntegerType* integer_ptr ) {
			assert( end >= buffer + sizeof( IntegerType )) ;
			*integer_ptr = 0 ;
			for( std::size_t byte_i = 0; byte_i < sizeof( IntegerType ); ++byte_i ) {
				(*integer_ptr) |= IntegerType( *reinterpret_cast< unsigned char const* >( buffer++ )) << ( 8 * byte_i ) ;
			}
			return buffer ;
		}

		// Read an integer stored in little-endian format into an integer stored in memory.
		// The stream is assumed to have sizeof( Integertype ) readable bytes.
		template< typename IntegerType >
		void read_little_endian_integer( std::istream& in_stream, IntegerType* integer_ptr ) {
			char buffer[ sizeof( IntegerType ) ] ;
			in_stream.read( buffer, sizeof( IntegerType )) ;
			if( in_stream ) {
				read_little_endian_integer( buffer, buffer + sizeof( IntegerType ), integer_ptr ) ;
			}
		}

		template< typename IntegerType >
		void read_length_followed_by_data( std::istream& in_stream, IntegerType* length_ptr, std::string* string_ptr ) {
			IntegerType& length = *length_ptr ;
			read_little_endian_integer( in_stream, length_ptr ) ;
			std::vector< char >buffer ( length ) ;
			in_stream.read( &buffer[0], length ) ;
			string_ptr->assign( buffer.begin(), buffer.end() ) ;
		}

		// Write an integer to the buffer in little-endian format.
		template< typename IntegerType >
		char* write_little_endian_integer( char* buffer, char* const end, IntegerType const integer ) {
			assert( end >= buffer + sizeof( IntegerType )) ;
			for( std::size_t byte_i = 0; byte_i < sizeof( IntegerType ); ++byte_i ) {
				*buffer++ = ( integer >> ( 8 * byte_i ) ) & 0xff ;
			}
			return buffer ;
		}
		
		// Write an integer to the stream in little-endian format.
		// The stream is assumed to have sizeof( Integertype ) bytes writeable.
		template< typename IntegerType >
		void write_little_endian_integer( std::ostream& out_stream, IntegerType const integer ) {
			char buffer[ sizeof( IntegerType ) ] ;
			write_little_endian_integer( buffer, buffer + sizeof( IntegerType ), integer ) ;
			out_stream.write( buffer, sizeof( IntegerType )) ;
		}
		
		template< typename IntegerType >
		void write_length_followed_by_data( std::ostream& out_stream, IntegerType length, std::string const data_string ) {
			assert( length <= data_string.size() ) ;
			write_little_endian_integer( out_stream, length ) ;
			out_stream.write( data_string.data(), length ) ;
		}

		template< typename SampleSetter >
		std::size_t read_sample_identifier_block(
			std::istream& aStream,
			Context const& context,
			SampleSetter setter
		) {
			uint32_t block_size = 0 ;
			uint32_t number_of_samples = 0 ;
			uint16_t identifier_size ;
			std::string identifier ;
			std::size_t bytes_read = 0 ;

			read_little_endian_integer( aStream, &block_size ) ;
			read_little_endian_integer( aStream, &number_of_samples ) ;
			bytes_read += 8 ;
			assert( number_of_samples == context.number_of_samples ) ;

			for( uint32_t i = 0; i < number_of_samples; ++i ) {
				read_length_followed_by_data( aStream, &identifier_size, &identifier ) ;
				if( aStream ) {
					bytes_read += sizeof( identifier_size ) + identifier_size ;
					setter( identifier ) ;
				} else {
					throw BGenError() ;
				}
			}
			assert( bytes_read == block_size ) ;
			return bytes_read ;
		}

		template<
			typename NumberOfAllelesSetter,
			typename AlleleSetter
		>
		bool read_snp_identifying_data(
			std::istream& aStream,
			Context const& context,
			std::string* SNPID,
			std::string* RSID,
			std::string* chromosome,
			uint32_t* SNP_position,
			NumberOfAllelesSetter set_number_of_alleles,
			AlleleSetter set_allele
		) {
			uint16_t SNPID_size = 0;
			uint16_t RSID_size = 0;
			uint16_t numberOfAlleles = 0 ;
			uint16_t chromosome_size = 0 ;
			uint32_t allele_size = 0;
			std::string allele ;
			uint32_t const layout = context.flags & e_Layout ;
			
			if( layout == e_v11Layout || layout == e_v10Layout ) {
				uint32_t number_of_samples ;
				read_little_endian_integer( aStream, &number_of_samples ) ;
				if( !aStream ) {
					return false ;
				}
				if( number_of_samples != context.number_of_samples ) {
					throw BGenError() ;
				}
			}
			std::string chromosome_string ;
			read_length_followed_by_data( aStream, &SNPID_size, SNPID ) ;
			if( layout == e_v12Layout && !aStream ) {
				return false ;
			}
			read_length_followed_by_data( aStream, &RSID_size, RSID ) ;
			read_length_followed_by_data( aStream, &chromosome_size, chromosome ) ;
			read_little_endian_integer( aStream, SNP_position ) ;
			if( layout == e_v12Layout ) {
				read_little_endian_integer( aStream, &numberOfAlleles ) ;
			} else {
				numberOfAlleles = 2 ;
			}
			set_number_of_alleles( numberOfAlleles ) ;
			for( uint16_t i = 0; i < numberOfAlleles; ++i ) {
				read_length_followed_by_data( aStream, &allele_size, &allele ) ;
				set_allele( i, allele ) ;
			}
			if( !aStream ) {
				throw BGenError() ;
			}
			return true ;
		}

		namespace v11 {
			namespace impl {
				template< typename FloatType >
				FloatType convert_from_integer_representation( uint16_t number, FloatType factor ) {
					FloatType result = number ;
					result /= factor ;
					return result ;
				}

				template< typename FloatType >
				uint16_t convert_to_integer_representation( FloatType number, FloatType factor ) {
					number *= factor ;
					number = std::min( std::max( number, 0.0 ), 65535.0 ) ;
					return static_cast< uint16_t > ( std::floor( number + 0.5 ) ) ;
				}

				double get_probability_conversion_factor( uint32_t flags ) ;
			}

			template< typename GenotypeProbabilityGetter >
			char* write_uncompressed_snp_probability_data(
				char* buffer,
				char* const end,
				Context const& context,
				GenotypeProbabilityGetter get_AA_probability,
				GenotypeProbabilityGetter get_AB_probability,
				GenotypeProbabilityGetter get_BB_probability
			) {
				double const factor = impl::get_probability_conversion_factor( context.flags ) ;
				for ( uint32_t i = 0 ; i < context.number_of_samples ; ++i ) {
					uint16_t
						AA = impl::convert_to_integer_representation( get_AA_probability( i ), factor ),
						AB = impl::convert_to_integer_representation( get_AB_probability( i ), factor ),
						BB = impl::convert_to_integer_representation( get_BB_probability( i ), factor ) ;
					assert( ( buffer + 6 ) <= end ) ;
					buffer = write_little_endian_integer( buffer, end, AA ) ;
					buffer = write_little_endian_integer( buffer, end, AB ) ;
					buffer = write_little_endian_integer( buffer, end, BB ) ;
				}
				return buffer ;
			}

			template< typename Setter >
			void parse_probability_data(
				char const* buffer,
				char const* const end,
				Context const& context,
				Setter& setter
			) {
				setter.initialise( context.number_of_samples, 2 ) ;
				
				double const probability_conversion_factor = impl::get_probability_conversion_factor( context.flags ) ;
				for ( uint32_t i = 0 ; i < context.number_of_samples ; ++i ) {
					setter.set_sample( i ) ;
					setter.set_number_of_entries( 3, ePerUnorderedGenotype, eProbability ) ;
					assert( end >= buffer + 6 ) ;
					for( std::size_t g = 0; g < 3; ++g ) {
						uint16_t prob ;
						buffer = read_little_endian_integer( buffer, end, &prob ) ;
						setter( impl::convert_from_integer_representation( prob, probability_conversion_factor ) ) ;
					}
				}
			}
		}

		namespace v12{
			namespace impl {
				// utility function to fill
				char const* read_bits_from_buffer(
					char const* buffer,
					char const* const end,
					uint64_t* data,
					int* size,
					uint8_t const bits
				) ;

				double parse_bit_representation(
					uint64_t* data,
					int* size,
					int const bits
				) ;
					
				// n choose k implementation
				// a faster implementation is of course possible, (e.g. table lookup)
				// but this is not a bottleneck.
				uint32_t n_choose_k( uint32_t n, uint32_t k ) ;

				void round_probs_to_scaled_simplex( double* p, std::size_t* index, std::size_t const n, int const number_of_bits ) ;

				// Write data encoding n probabilities, given in probs, that sum to 1,
				// starting at the given offset in data, to the given buffer.
				char* write_scaled_probs(
					uint64_t* data,
					std::size_t* offset,
					double const* probs,
					std::size_t const n,
					int const number_of_bits,
					char* buffer,
					char* const end
				) ;
			}

			template< typename Setter >
			void parse_probability_data(
				char const* buffer,
				char const* const end,
				Context const& context,
				Setter& setter
			) {
				uint32_t numberOfSamples ;
				uint16_t numberOfAlleles ;
				unsigned char ploidyExtent[2] ;
				enum { ePhased = 1, eUnphased = 0 } ;
				
				buffer = read_little_endian_integer( buffer, end, &numberOfSamples ) ;
				buffer = read_little_endian_integer( buffer, end, &numberOfAlleles ) ;
				buffer = read_little_endian_integer( buffer, end, &ploidyExtent[0] ) ;
				buffer = read_little_endian_integer( buffer, end, &ploidyExtent[1] ) ;

				if( numberOfSamples != context.number_of_samples ) {
					throw BGenError() ;
				}
				setter.initialise( numberOfSamples, uint32_t( numberOfAlleles ) ) ;

				// Keep a pointer to the ploidy and move buffer past the ploidy information
				char const* ploidy_p = buffer ;
				buffer += numberOfSamples ;
				// Get the phased flag and number of bits
				bool const phased = ((*buffer++) & 0x1 ) ;
				int const bits = int( *reinterpret_cast< unsigned char const *>( buffer++ ) ) ;
			
	#if DEBUG_BGEN_FORMAT
				std::cerr << "parse_probability_data_v12(): numberOfSamples = " << numberOfSamples
					<< ", phased = " << phased << ".\n" ;
				std::cerr << "parse_probability_data_v12(): *buffer: "
					<< bgen::impl::to_hex( buffer, end ) << ".\n" ;
	#endif

				{
					uint64_t data = 0 ;
					int size = 0 ;
					for( uint32_t i = 0; i < numberOfSamples; ++i, ++ploidy_p ) {
						uint32_t const ploidy = uint32_t(*reinterpret_cast< unsigned char const* >( ploidy_p ) & 0x3F) ;
						bool const missing = (*reinterpret_cast< unsigned char const* >( ploidy_p ) & 0x80) ;
						uint32_t const valueCount
							= phased
							? (ploidy * numberOfAlleles)
							: ((numberOfAlleles == 2) ? (ploidy+1) : impl::n_choose_k( ploidy + numberOfAlleles - 1, numberOfAlleles - 1 )) ;

						uint32_t const storedValueCount = valueCount - ( phased ? ploidy : 1 ) ;
					
	#if DEBUG_BGEN_FORMAT > 1
						std::cerr << "parse_probability_data_v12(): sample " << i
							<< ", ploidy = " << ploidy
							<< ", missing = " << missing
							<< ", valueCount = " << valueCount
							<< ", storedValueCount = " << storedValueCount
							<< ", data = " << bgen::impl::to_hex( buffer, end )
							<< ".\n" ;
	#endif

						if( setter.set_sample( i ) ) {
							setter.set_number_of_entries(
								valueCount, 
								phased ? ePerPhasedHaplotypePerAllele : ePerUnorderedGenotype,
								eProbability
							) ;
							if( missing ) {
								for( std::size_t h = 0; h < valueCount; ++h ) {
									setter( genfile::MissingValue() ) ;
								}
							} else {
								double sum = 0.0 ;
								for( uint32_t h = 0; h < storedValueCount; ++h ) {
									buffer = impl::read_bits_from_buffer( buffer, end, &data, &size, bits ) ;
									double const value = impl::parse_bit_representation( &data, &size, bits ) ;
									setter( value ) ;
									sum += value ;
	#if DEBUG_BGEN_FORMAT
									std::cerr << "parse_probability_data_v12(): i = " << i << ", h = " << h << ", size = " << size << ", bits = " << bits << ", parsed value = " << value
										<< ", sum = " << std::setprecision( 10 ) << sum << ".\n" ;
	#endif
								
									if(
										( phased && ((h+1) % (numberOfAlleles-1) ) == 0 )
										|| ((!phased) && (h+1) == storedValueCount )
									) {
										assert( sum <= 1.00000001 ) ;
										setter( 1.0 - sum ) ;
										sum = 0.0 ;
									}
								}
							}
						} else {
							// just consume data, don't set anything.
							for( uint32_t h = 0; h < storedValueCount; ++h ) {
								buffer = impl::read_bits_from_buffer( buffer, end, &data, &size, bits ) ;
								impl::parse_bit_representation( &data, &size, bits ) ;
							}
						}
					}
				}
			}

			template< typename GenotypeProbabilityGetter >
			char* write_uncompressed_snp_probability_data(
				char* buffer,
				char* const end,
				Context const& context,
				GenotypeProbabilityGetter get_AA_probability,
				GenotypeProbabilityGetter get_AB_probability,
				GenotypeProbabilityGetter get_BB_probability,
				int const number_of_bits
			) {
				assert( number_of_bits > 0 ) ;
				assert( number_of_bits <= 32 ) ;
#if DEBUG_BGEN_FORMAT
				std::cerr << "genfile::bgen::impl::v12::write_uncompressed_snp_probability_data(): number_of_bits = " << number_of_bits << ", buffer = " << reinterpret_cast< void* >( buffer ) << ", (end-buffer) = " << (end-buffer) << ".\n" ;
#endif
				buffer = write_little_endian_integer( buffer, end, context.number_of_samples ) ;
				// Write ploidy
				uint16_t const numberOfAlleles = 2 ;
				uint8_t const ploidy = 2 ;
				buffer = write_little_endian_integer( buffer, end, numberOfAlleles ) ;
				buffer = write_little_endian_integer( buffer, end, ploidy ) ;
				buffer = write_little_endian_integer( buffer, end, ploidy ) ;
				char* ploidy_p = buffer ;
				buffer += context.number_of_samples ;
				buffer = write_little_endian_integer( buffer, end, uint8_t( 0 ) ) ;
				buffer = write_little_endian_integer( buffer, end, uint8_t( number_of_bits ) ) ;
				// We use an array of three doubles to compute the rounded probabilities.
				// We use a single 64-bit integer to marshall the data to be written.
				double v[3] ;
				std::size_t index[3] ;
				uint64_t data = 0 ;
				std::size_t offset = 0 ;
				for( std::size_t i = 0; i < context.number_of_samples; ++i ) {
					v[0] = get_AA_probability(i) ;
					v[1] = get_AB_probability(i) ;
					v[2] = get_BB_probability(i) ;
					double sum = v[0] + v[1] + v[2] ;
					bool missing = ( sum == 0 ) ;
					v[0] = v[0] / sum ;
					v[1] = v[1] / sum ;
					v[2] = v[2] / sum ;
					uint8_t ploidy = 2 | ( missing ? 0x80 : 0 ) ;
					*(ploidy_p++) = ploidy ;

					if( !missing ) {
						impl::round_probs_to_scaled_simplex( &v[0], &index[0], 3, number_of_bits ) ;
						
						buffer = impl::write_scaled_probs( &data, &offset, &v[0], 3, number_of_bits, buffer, end ) ;
#if DEBUG_BGEN_FORMAT
						double sum = v[0] + v[1] + v[2] ;
						std::cerr << "genfile::bgen::impl::v12::write_uncompressed_snp_probability_data(): scaled probs are:"
							<< v[0] << ", " << v[1] << ", " << v[2] << ".\n" ;
						std::cerr << "genfile::bgen::impl::v12::write_uncompressed_snp_probability_data(): sum is:"
							<< sum << ".\n" ;
						std::cerr << "genfile::bgen::impl::v12::write_uncompressed_snp_probability_data(): after write, data = "
							<< bgen::impl::to_hex(
								reinterpret_cast< unsigned char const* >( &data ),
								reinterpret_cast< unsigned char const* >( &data ) + 8
							) << ".\n" ;
#endif
					}
				}
				// Get any leftover bytes.
				if( offset > 0 ) {
					int const nBytes = (offset+7)/8 ;
#if DEBUG_BGEN_FORMAT
					std::cerr << "genfile::bgen::impl::v12::write_uncompressed_snp_probability_data(): final offset = "
						<< offset << ", number of bits = " << number_of_bits << ", writing " << nBytes << " final bytes (space = " << (end-buffer) << ".\n" ;
#endif
					assert( (buffer+nBytes) <= end ) ;
					buffer = std::copy(
						reinterpret_cast< char const* >( &data ),
						reinterpret_cast< char const* >( &data ) + nBytes,
						buffer
					) ;
				}
				return buffer ;
			}
		}
		
		template< typename Setter >
		void parse_probability_data(
			char const* buffer,
			char const* const end,
			Context const& context,
			Setter& setter
		) {
			if( (context.flags & e_Layout) == e_v10Layout || (context.flags & e_Layout) == e_v11Layout ) {
				v11::parse_probability_data( buffer, end, context, setter ) ;
			} else {
				v12::parse_probability_data( buffer, end, context, setter ) ;
			}
		}

		template< typename Setter >
		void read_and_parse_genotype_data_block(
			std::istream& aStream,
			Context const& context,
			Setter& setter,
			std::vector< char >* buffer1,
			std::vector< char >* buffer2
		) {
			read_genotype_data_block( aStream, context, buffer1 ) ;
			uncompress_probability_data( context, *buffer1, buffer2 ) ;
			parse_probability_data(
				&(*buffer2)[0],
				&(*buffer2)[0] + buffer2->size(),
				context,
				setter
			) ;
		}

		template< typename GenotypeProbabilityGetter >
		char* write_uncompressed_snp_probability_data(
			char* buffer,
			char* const bufferEnd,
			Context const& context,
			GenotypeProbabilityGetter get_AA_probability,
			GenotypeProbabilityGetter get_AB_probability,
			GenotypeProbabilityGetter get_BB_probability,
			int const number_of_bits = 16
		) {
			uint32_t const layout = context.flags & e_Layout ;
			if( layout == e_v11Layout ) {
				buffer = v11::write_uncompressed_snp_probability_data(
					buffer,
					bufferEnd,
					context,
					get_AA_probability, get_AB_probability, get_BB_probability
				) ;
				assert( buffer == bufferEnd ) ;
			} else if( layout == e_v12Layout ) {
				buffer = v12::write_uncompressed_snp_probability_data(
					buffer,
					bufferEnd,
					context,
					get_AA_probability, get_AB_probability, get_BB_probability,
					number_of_bits
				) ;
			} else {
				assert(0) ;
			}
			return buffer ;
		}

		template< typename GenotypeProbabilityGetter >
		void write_snp_probability_data(
			std::ostream& aStream,
			Context const& context,
			GenotypeProbabilityGetter get_AA_probability,
			GenotypeProbabilityGetter get_AB_probability,
			GenotypeProbabilityGetter get_BB_probability,
			int const number_of_bits,
			std::vector< char >* buffer,
			std::vector< char >* compression_buffer
		) {
			uint32_t const layout = context.flags & e_Layout ;
			// Write the data to a buffer, which we then compress and write to the stream
			uLongf uncompressed_data_size =
				( layout == e_v11Layout )
					? (6 * context.number_of_samples)
					: ( 10 + context.number_of_samples + ((( context.number_of_samples * number_of_bits * 2 )+7) / 8) ) ;

			buffer->resize( uncompressed_data_size ) ;
			char* p = write_uncompressed_snp_probability_data(
				&(*buffer)[0],
				&(*buffer)[0] + uncompressed_data_size,
				context,
				get_AA_probability, get_AB_probability, get_BB_probability,
				number_of_bits
			) ;
			assert( p = &(*buffer)[0] + uncompressed_data_size ) ;

			if( context.flags & e_CompressedSNPBlocks ) {
	#if HAVE_ZLIB
				uLongf compression_buffer_size = 12 + (1.1 * uncompressed_data_size) ;		// calculated according to zlib manual.
				compression_buffer->resize( compression_buffer_size ) ;
				int result = compress(
					reinterpret_cast< Bytef* >( &(*compression_buffer)[0] ), &compression_buffer_size,
					reinterpret_cast< Bytef* >( &(*buffer)[0] ), uncompressed_data_size
				) ;
				assert( result == Z_OK ) ;
				// write total payload size (compression_buffer_size is now the compressed length of the data).
				// Account for a 4-byte uncompressed data size if we are in layout 1.2
				if( layout == e_v12Layout ) {
					write_little_endian_integer( aStream, uint32_t( compression_buffer_size ) + 4 ) ;
					write_little_endian_integer( aStream, uint32_t( uncompressed_data_size )) ;
				} else {
					write_little_endian_integer( aStream, uint32_t( compression_buffer_size ) ) ;
				}

				// and write the data
				aStream.write( &(*compression_buffer)[0], compression_buffer_size ) ;
	#else
				assert(0) ;
	#endif
			} else {
				if( layout == e_v12Layout ) {
					write_little_endian_integer( aStream, uint32_t( uncompressed_data_size )) ;
				}
				aStream.write( &(*buffer)[0], uncompressed_data_size ) ;
			}
		}
	}
}

#endif
