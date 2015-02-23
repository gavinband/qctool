
//			Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//	  (See accompanying file LICENSE_1_0.txt or copy at
//			http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <string>
#include "genfile/endianness_utils.hpp"
#include "genfile/bgen/bgen.hpp"
#include "genfile/bgen/impl.hpp"
#include "genfile/string_utils/hex.hpp"

// #define DEBUG_BGEN_FORMAT 2

namespace genfile {
	namespace bgen {
		namespace impl {
			void check_for_two_alleles( uint16_t numberOfAlleles ) {
				if( numberOfAlleles != 2 ) {
					std::cerr << "genfile::bgen::impl::check_for_two_alleles: only biallelic variants are currently supported.\n" ;
					assert(0) ;
				}
			}
			
			struct TwoAlleleSetter {
				TwoAlleleSetter( std::string* allele1, std::string* allele2 ):
					m_allele1( allele1 ),
					m_allele2( allele2 )
				{
					assert( allele1 != 0 ) ;
					assert( allele2 != 0 ) ;
				}

				void operator()( uint16_t i, std::string const& value ) {
					if( i == 0 ) {
						*m_allele1 = value ;
					} else if( i == 1 ) {
						*m_allele2 = value ;
					} else {
						assert(0) ;
					}
				}
				std::string* m_allele1 ;
				std::string* m_allele2 ;
			} ;
		}


		void read_offset( std::istream& iStream, uint32_t* offset ) {
			read_little_endian_integer( iStream, offset ) ;
		}

		void write_offset( std::ostream& oStream, uint32_t const offset ) {
			write_little_endian_integer( oStream, offset ) ;
		}
		
		std::size_t read_header_block(
			std::istream& aStream,
			Context* context
		) {
			assert( context != 0 ) ;
			uint32_t
				header_size = 0,
				number_of_snp_blocks = 0,
				number_of_samples = 0,
				flags = 0 ;

			char magic[4] ;
			std::size_t fixed_data_size = 20 ;
			std::vector<char> free_data ;

			read_little_endian_integer( aStream, &header_size ) ;
			assert( header_size >= fixed_data_size ) ;
			read_little_endian_integer( aStream, &number_of_snp_blocks ) ;
			read_little_endian_integer( aStream, &number_of_samples ) ;
			aStream.read( &magic[0], 4 ) ;
			free_data.resize( header_size - fixed_data_size ) ;
			aStream.read( &free_data[0], free_data.size() ) ;
			read_little_endian_integer( aStream, &flags ) ;

			if(
				( magic[0] != 'b' || magic[1] != 'g' || magic[2] != 'e' || magic[3] != 'n' )
				&& ( magic[0] != 0 || magic[1] != 0 || magic[2] != 0 || magic[3] != 0 )
			) {
				throw BGenError() ;
			}

			if( aStream ) {
				context->number_of_samples = number_of_samples ;
				context->number_of_variants = number_of_snp_blocks ;
				context->magic.assign( &magic[0], &magic[0] + 4 ) ;
				context->free_data.assign( free_data.begin(), free_data.end() ) ;
				context->flags = flags ;

				return( header_size ) ;
			} else {
				throw BGenError() ;
			}
		}
		

		bool read_snp_identifying_data(
			std::istream& aStream,
			Context const& context,
			std::string* SNPID,
			std::string* RSID,
			std::string* chromosome,
			uint32_t* SNP_position,
			std::string* first_allele,
			std::string* second_allele
		) {
#if DEBUG_BGEN_FORMAT
			std::cerr << "genfile::bgen::impl::read_snp_identifying_data(): flags = 0x" << std::hex << context.flags << ".\n" ;
#endif
			uint32_t const layout = context.flags & e_Layout ;
			if( layout == e_v11Layout || layout == e_v12Layout ) {
				// forward to v12 version which handles multiple alleles.
				impl::TwoAlleleSetter allele_setter( first_allele, second_allele ) ;
				return bgen::read_snp_identifying_data(
					aStream, context,
					SNPID, RSID, chromosome, SNP_position,
					&impl::check_for_two_alleles,
					allele_setter
				) ;
			} else if( layout == e_v10Layout ) {
				// v1.0-style layout, deprecated
				uint32_t number_of_samples ;
				unsigned char max_id_size = 0;
				unsigned char SNPID_size = 0;
				unsigned char RSID_size = 0;
				if( aStream ) {
					genfile::read_little_endian_integer( aStream, &number_of_samples ) ;
					if( !aStream ) {
						return false ;
					}
					if( number_of_samples != context.number_of_samples ) {
						throw BGenError() ;
					}
				}
				if( aStream ) {
					genfile::read_little_endian_integer( aStream, &max_id_size ) ;
				}
				if( aStream ) {
					read_length_followed_by_data( aStream, &SNPID_size, SNPID ) ;
					assert( SNPID_size <= max_id_size ) ;
					aStream.ignore( max_id_size - SNPID_size ) ;
				}
				if( aStream ) {
					read_length_followed_by_data( aStream, &RSID_size, RSID ) ;
					assert( RSID_size <= max_id_size ) ;
					aStream.ignore( max_id_size - RSID_size ) ;
				}
				if( aStream ) {
					unsigned char chromosome_char ;
					genfile::read_little_endian_integer( aStream, &chromosome_char ) ;
					genfile::read_little_endian_integer( aStream, SNP_position ) ;

					switch( chromosome_char ) {
						case 1: *chromosome = "01" ; break ;
						case 2: *chromosome = "02" ; break ;
						case 3: *chromosome = "03" ; break ;
						case 4: *chromosome = "04" ; break ;
						case 5: *chromosome = "05" ; break ;
						case 6: *chromosome = "06" ; break ;
						case 7: *chromosome = "07" ; break ;
						case 8: *chromosome = "08" ; break ;
						case 9: *chromosome = "09" ; break ;
						case 10: *chromosome = "10" ; break ;
						case 11: *chromosome = "11" ; break ;
						case 12: *chromosome = "12" ; break ;
						case 13: *chromosome = "13" ; break ;
						case 14: *chromosome = "14" ; break ;
						case 15: *chromosome = "15" ; break ;
						case 16: *chromosome = "16" ; break ;
						case 17: *chromosome = "17" ; break ;
						case 18: *chromosome = "18" ; break ;
						case 19: *chromosome = "19" ; break ;
						case 20: *chromosome = "20" ; break ;
						case 21: *chromosome = "21" ; break ;
						case 22: *chromosome = "22" ; break ;
						case 23: *chromosome = "0X" ; break ;
						case 24: *chromosome = "0Y" ; break ;
						case 253: *chromosome = "XY" ; break ;
						case 254: *chromosome = "MT" ; break ;
						default: *chromosome = "NA" ;
					}
					unsigned char allele1_size = 0, allele2_size = 0 ;
					genfile::read_length_followed_by_data( aStream, &allele1_size, first_allele ) ;
					genfile::read_length_followed_by_data( aStream, &allele2_size, second_allele ) ;
				}
			} else {
				assert(0) ;
			}
			if( !aStream ) {
				throw BGenError() ;
			}
			return true ;
		}
		
		void write_snp_identifying_data(
			std::ostream& aStream,
			Context const& context,
			unsigned char max_id_size,
			std::string SNPID,
			std::string RSID,
			unsigned char chromosome,
			uint32_t SNP_position,
			std::string first_allele,
			std::string second_allele
		) {
			uint32_t const layout = context.flags & e_Layout ;
			assert( layout == e_v11Layout || layout == e_v12Layout ) ;

			if( layout == e_v11Layout ) {
				write_little_endian_integer( aStream, context.number_of_samples ) ;
				// otherwise this is moved to the probability data block.
			}

			std::size_t const max_allele_length = std::numeric_limits< uint32_t >::max() ;
			std::size_t const max_id_length = std::numeric_limits< uint16_t >::max() ;
			assert( SNPID.size() <= static_cast< std::size_t >( max_id_length )) ;
			assert( RSID.size() <= static_cast< std::size_t >( max_id_length )) ;
			if( first_allele.size() > static_cast< std::size_t >( max_allele_length ) ) {
				std::cerr << "Warning: at SNP " << SNPID << " " << RSID << " pos=" << SNP_position << ", truncating first allele of size " << first_allele.size() << ".\n" ;
				first_allele.resize( max_allele_length - 3 ) ;
				first_allele += "..." ;
			}
			if( second_allele.size() > static_cast< std::size_t >( max_allele_length ) ) {
				std::cerr << "Warning: at SNP " << SNPID << " " << RSID << " pos=" << SNP_position << ", truncating second allele of size " << second_allele.size() << ".\n" ;
				second_allele.resize( max_allele_length - 3 ) ;
				second_allele += "..." ;
			}
			write_length_followed_by_data( aStream, uint16_t( SNPID.size() ), SNPID.data() ) ;
			write_length_followed_by_data( aStream, uint16_t( RSID.size() ), RSID.data() ) ;
			std::string const chromosome_string = Chromosome( chromosome ) ;
			write_length_followed_by_data( aStream, uint16_t( chromosome_string.size() ), chromosome_string ) ;
			write_little_endian_integer( aStream, SNP_position ) ;
			
			if( layout == e_v12Layout ) {
				// v12 has an explicit allele count, here equal to 2.
				write_little_endian_integer( aStream, uint16_t(2) ) ;
			}
			write_length_followed_by_data( aStream, uint32_t( first_allele.size() ), first_allele.data() ) ;
			write_length_followed_by_data( aStream, uint32_t( second_allele.size() ), second_allele.data() ) ;
		}

		namespace impl {
			double get_probability_conversion_factor( uint32_t flags ) {
				uint32_t layout = flags & e_Layout ;
				if( layout == e_v10Layout ) {
					// v1.0-style blocks, deprecated
					return 10000.0 ;
				} else if( layout == e_v11Layout ) {
					// v1.1-style blocks
					return 32768.0 ;
				} else {
					// v1.2 style (or other) blocks, these are treated differently and this function does not apply.
					assert(0) ;
				}
			}
		}

		void write_header_block(
			std::ostream& aStream,
			Context const& context
		) {
			uint32_t header_size = context.header_size() ;
			genfile::write_little_endian_integer( aStream, header_size ) ;
			genfile::write_little_endian_integer( aStream, context.number_of_variants ) ;
			genfile::write_little_endian_integer( aStream, context.number_of_samples ) ;
			aStream.write( context.magic.data(), 4 ) ;
			aStream.write( context.free_data.data(), context.free_data.size() ) ;
			genfile::write_little_endian_integer( aStream, context.flags ) ;
		}
		
		std::size_t write_sample_identifier_block(
			std::ostream& aStream,
			Context const& context,
			std::vector< std::string > const& sample_ids
		) {
			assert( sample_ids.size() == context.number_of_samples ) ;
			uint32_t block_size = 8 ;
			for( uint32_t i = 0; i < sample_ids.size(); ++i ) {
				block_size += 2 + sample_ids[i].size() ;
			}
			write_little_endian_integer( aStream, block_size ) ;
			write_little_endian_integer( aStream, context.number_of_samples ) ;
#if DEBUG_BGEN_FORMAT
			std::cerr << "genfile::bgen::write_sample_identifier_block(): writing " << sample_ids.size() << " samples...\n" ;
#endif
			for( uint32_t i = 0; i < sample_ids.size(); ++i ) {
#if DEBUG_BGEN_FORMAT
				std::cerr << "genfile::bgen::write_sample_identifier_block(): writing sample " << sample_ids[i] << ".\n" ;
#endif
				std::string const& identifier = sample_ids[i] ;
				assert( identifier.size() <= std::size_t( std::numeric_limits< uint16_t >::max() ) ) ;
				uint16_t const id_size = identifier.size() ;
				write_length_followed_by_data( aStream, id_size, identifier ) ;
			}
			return block_size ;
		}
		
		void ignore_snp_probability_data(
			std::istream& aStream,
			Context const& context
		) {
			if( context.flags & bgen::e_CompressedSNPBlocks ) {
				uint32_t compressed_data_size ;
				genfile::read_little_endian_integer( aStream, &compressed_data_size ) ;
				aStream.ignore( compressed_data_size ) ;
			}
			else {
				aStream.ignore( 6 * context.number_of_samples ) ;
			}
		}

		void read_uncompressed_snp_probability_data_v12(
			char const* buffer,
			char const* const end,
			Context const& context,
			VariantDataReader::PerSampleSetter& setter
		) ;

		void read_uncompressed_snp_probability_data(
			char const* buffer,
			char const* const end,
			Context const& context,
			VariantDataReader::PerSampleSetter& setter
		) {
			if( (context.flags & e_Layout) == e_v10Layout || (context.flags & e_Layout) == e_v11Layout ) {
				setter.set_number_of_samples( context.number_of_samples ) ;
				setter.set_number_of_alleles( 2 ) ;
				
				double const probability_conversion_factor = impl::get_probability_conversion_factor( context.flags ) ;
				for ( uint32_t i = 0 ; i < context.number_of_samples ; ++i ) {
					setter.set_sample( i ) ;
					setter.set_number_of_entries( 3 ) ;
					setter.set_order_type(
						VariantDataReader::PerSampleSetter::ePerUnorderedGenotype,
						VariantDataReader::PerSampleSetter::eProbability
					) ;
					assert( end >= buffer + 6 ) ;
					for( std::size_t g = 0; g < 3; ++g ) {
						uint16_t prob ;
						buffer = read_little_endian_integer( buffer, end, &prob ) ;
						setter( impl::v11::convert_from_integer_representation( prob, probability_conversion_factor ) ) ;
					}
				}
			} else {
				read_uncompressed_snp_probability_data_v12(
					buffer, end, context, setter
				) ;
			}
		}
		
		namespace impl {
			uint32_t n_choose_k( uint32_t n, uint32_t k ) {
				if( k == 0 )  {
					return 1 ;
				}
				return ( n * n_choose_k(n - 1, k - 1) ) / k ;
			}
			
			namespace v12 {
				char const* fill_data(
					char const* buffer,
					char const* const end,
					uint64_t* data,
					int* size,
					uint8_t const bits
				) {
					// fill data with up to 8 bytes.
					while( (*size) < bits && buffer < end ) {
						(*data) |= uint64_t( *(reinterpret_cast< unsigned char const* >( buffer++ ))) << (*size) ;
						(*size) += sizeof( unsigned char ) * 8 ;
	#if DEBUG_BGEN_FORMAT > 1
						std::cerr << "genfile::impl::v12::fill_data(): size = " << (*size)
							<< ", buffer = " << reinterpret_cast< void const* >( buffer )
							<< ", data = "
							<< string_utils::to_hex(
								reinterpret_cast< unsigned char const* >( data ),
								reinterpret_cast< unsigned char const* >( data ) + 8
							) << ".\n" ;
	#endif
					}
					assert( (*size) >= bits ) ;
					return buffer ;
				}
			
				double consume_value(
					uint64_t* data,
					int* size,
					int const bits
				) {
#if DEBUG_BGEN_FORMAT > 1
					std::cerr << "genfile::impl::v12::consume_value(): size = " << (*size)
						<< ", data = "
						<< string_utils::to_hex(
							reinterpret_cast< unsigned char const* >( data ),
							reinterpret_cast< unsigned char const* >( data ) + 8
						) << ".\n" ;
#endif
					uint64_t bitMask = (0xFFFFFFFFFFFFFFFF >> ( 64 - bits )) ;
					double const result = ( *data & bitMask ) / double( bitMask ) ;
					(*size) -= bits ;
					(*data) >>= bits ;
					return result ; 
				}
			}
		}
		
		void read_uncompressed_snp_probability_data_v12(
			char const* buffer,
			char const* const end,
			Context const& context,
			VariantDataReader::PerSampleSetter& setter
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
			setter.set_number_of_samples( numberOfSamples ) ;
			setter.set_number_of_alleles( uint32_t( numberOfAlleles ) ) ;

			// Keep a pointer to the ploidy and move buffer past the ploidy information
			char const* ploidy_p = buffer ;
			buffer += numberOfSamples ;
			// Get the phased flag and number of bits
			bool const phased = ((*buffer++) & 0x1 ) ;
			int const bits = int( *reinterpret_cast< unsigned char const *>( buffer++ ) ) ;
			
#if DEBUG_BGEN_FORMAT
			std::cerr << "read_uncompressed_snp_probability_data_v12(): numberOfSamples = " << numberOfSamples
				<< ", phased = " << phased << ".\n" ;
			std::cerr << "read_uncompressed_snp_probability_data_v12(): *buffer: "
				<< string_utils::to_hex( buffer, end ) << ".\n" ;
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
						: impl::n_choose_k( ploidy + numberOfAlleles - 1, numberOfAlleles - 1 ) ;
					uint32_t const storedValueCount = valueCount - ( phased ? ploidy : 1 ) ;
					
#if DEBUG_BGEN_FORMAT > 1
					std::cerr << "read_uncompressed_snp_probability_data_v12(): sample " << i
						<< ", ploidy = " << ploidy
						<< ", missing = " << missing
						<< ", valueCount = " << valueCount
						<< ", storedValueCount = " << storedValueCount
						<< ", data = " << string_utils::to_hex( buffer, end )
						<< ".\n" ;
#endif

					setter.set_sample( i ) ;
					setter.set_number_of_entries( valueCount ) ;
					setter.set_order_type(
						(
							phased
								? VariantDataReader::PerSampleSetter::ePerPhasedHaplotypePerAllele
								: VariantDataReader::PerSampleSetter::ePerUnorderedGenotype
						),
						VariantDataReader::PerSampleSetter::eProbability
					) ;

					if( missing ) {
						for( std::size_t h = 0; h < valueCount; ++h ) {
							setter( genfile::MissingValue() ) ;
						}
					} else {
						double sum = 0.0 ;
						for( uint32_t h = 0; h < storedValueCount; ++h ) {
							buffer = impl::v12::fill_data( buffer, end, &data, &size, bits ) ;
							double const value = impl::v12::consume_value( &data, &size, bits ) ;
							setter( value ) ;
							sum += value ;
							
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
				}
			}
		}
		
		void read_snp_probability_data(
			std::istream& aStream,
			Context const& context,
			VariantDataReader::PerSampleSetter& setter,
			std::vector< char >* buffer1,
			std::vector< char >* buffer2
		) {
			uint32_t uncompressed_data_size = 0 ;
			if( (context.flags & e_Layout) == e_v12Layout ) {
				genfile::read_little_endian_integer( aStream, &uncompressed_data_size ) ;
			} else {
				uncompressed_data_size = 6 * context.number_of_samples ;
			}
			buffer1->resize( uncompressed_data_size ) ;
			if( context.flags & bgen::e_CompressedSNPBlocks ) {
				uint32_t compressed_data_size ;
				genfile::read_little_endian_integer( aStream, &compressed_data_size ) ;
				buffer2->resize( compressed_data_size ) ;
				aStream.read( &(*buffer2)[0], compressed_data_size ) ;
				zlib_uncompress( *buffer2, buffer1 ) ;
				assert( buffer1->size() == uncompressed_data_size ) ;
			}
			else {
				aStream.read( &(*buffer1)[0], uncompressed_data_size ) ;
			}
			read_uncompressed_snp_probability_data(
				&(*buffer1)[0],
				&(*buffer1)[0] + uncompressed_data_size,
				context,
				setter
			) ;
		}
	}
}
