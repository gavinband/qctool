
//			Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//	  (See accompanying file LICENSE_1_0.txt or copy at
//			http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <string>
#include "genfile/Chromosome.hpp"
#include "genfile/endianness_utils.hpp"
#include "genfile/bgen.hpp"
#include "genfile/bgen/impl.hpp"

#define DEBUG_BGEN_FORMAT 1

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
		
		std::size_t get_header_block_size(
			std::string const& free_data
		) {
			std::size_t fixed_data_size = 5 * sizeof( uint32_t ) ;
			return fixed_data_size + free_data.size() ;
		}
		
		void write_header_block(
			std::ostream& aStream,
			uint32_t number_of_snp_blocks,
			uint32_t number_of_samples,
			std::string const& free_data,
			uint32_t flags
		) {
			uint32_t reserved = 0u ;
			uint32_t header_size = get_header_block_size( free_data ) ;

			genfile::write_little_endian_integer( aStream, header_size ) ;
			genfile::write_little_endian_integer( aStream, number_of_snp_blocks ) ;
			genfile::write_little_endian_integer( aStream, number_of_samples ) ;
			genfile::write_little_endian_integer( aStream, reserved ) ;
			aStream.write( free_data.data(), free_data.size() ) ;
			genfile::write_little_endian_integer( aStream, flags ) ;
		}
		
		void read_snp_identifying_data(
			std::istream& aStream,
			uint32_t const flags,
			std::string* SNPID,
			std::string* RSID,
			unsigned char* chromosome,
			uint32_t* SNP_position,
			std::string* first_allele,
			std::string* second_allele
		) {
            uint32_t number_of_samples ;
			if( aStream ) {
				genfile::read_little_endian_integer( aStream, &number_of_samples ) ;
			}

#if DEBUG_BGEN_FORMAT
			std::cerr << "genfile::bgen::impl::read_snp_identifying_data(): flags = 0x" << std::hex << flags << ".\n" ;
#endif
			uint32_t const layout = flags & e_Layout ;
			if( layout == e_v11Layout ) {
				// v1.1-style layout
				uint16_t SNPID_size = 0;
				uint16_t RSID_size = 0;
				uint32_t allele1_size = 0;
				uint32_t allele2_size = 0;
				uint16_t chromosome_size = 0 ;
				std::string chromosome_string ;
				read_length_followed_by_data( aStream, &SNPID_size, SNPID ) ;
				read_length_followed_by_data( aStream, &RSID_size, RSID ) ;
				read_length_followed_by_data( aStream, &chromosome_size, &chromosome_string ) ;
				genfile::read_little_endian_integer( aStream, SNP_position ) ;
				read_length_followed_by_data( aStream, &allele1_size, first_allele ) ;
				read_length_followed_by_data( aStream, &allele2_size, second_allele ) ;
				*chromosome = ChromosomeEnum( Chromosome( chromosome_string ) ) ;
			}
			else if( layout == e_v10Layout ) {
				// v1.0-style layout, deprecated
				unsigned char max_id_size = 0;
				unsigned char SNPID_size = 0;
				unsigned char RSID_size = 0;
				
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
					genfile::read_little_endian_integer( aStream, chromosome ) ;
					genfile::read_little_endian_integer( aStream, SNP_position ) ;

					first_allele->resize( 1 ) ;
					(*first_allele)[0] = aStream.get() ;
					second_allele->resize( 1 ) ;
					(*second_allele)[0] = aStream.get() ;
				}
			} else if( layout == e_v12Layout ) {
                // forward to v12 version which handles multiple alleles.
                impl::TwoAlleleSetter allele_setter( first_allele, second_allele ) ;
                
                bgen::read_snp_identifying_data_v12(
                    aStream, flags,
                    SNPID, RSID, chromosome, SNP_position,
                    &impl::check_for_two_alleles,
                    allele_setter
                ) ;
			} else {
			    assert(0) ;
			}
		}
        
		void write_snp_identifying_data(
			std::ostream& aStream,
			uint32_t const flags,
			uint32_t number_of_samples,
			unsigned char max_id_size,
			std::string SNPID,
			std::string RSID,
			unsigned char chromosome,
			uint32_t SNP_position,
			std::string first_allele,
			std::string second_allele
		) {
			write_little_endian_integer( aStream, number_of_samples ) ;
			uint32_t const layout = flags & e_Layout ;

			if( layout == e_v11Layout ) {
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
				write_length_followed_by_data( aStream, uint32_t( first_allele.size() ), first_allele.data() ) ;
				write_length_followed_by_data( aStream, uint32_t( second_allele.size() ), second_allele.data() ) ;
			}
			else if( layout == e_v10Layout ) {
				// v1.0-style layout is deprecated.
				assert( SNPID.size() <= static_cast< std::size_t >( max_id_size )) ;
				assert( RSID.size() <= static_cast< std::size_t >( max_id_size )) ;
				unsigned char SNPID_size = SNPID.size() ;
				unsigned char RSID_size = RSID.size() ;
				SNPID.resize( max_id_size, ' ' ) ;
				RSID.resize( max_id_size, ' ' ) ;

				write_little_endian_integer( aStream, max_id_size ) ;
				write_length_followed_by_data( aStream, SNPID_size, SNPID.data() ) ;
				aStream.write( SNPID.data() + SNPID_size, max_id_size - SNPID_size ) ;
				write_length_followed_by_data( aStream, RSID_size, RSID.data() ) ;
				aStream.write( RSID.data() + RSID_size, max_id_size - RSID_size ) ;
				write_little_endian_integer( aStream, chromosome ) ;
				write_little_endian_integer( aStream, SNP_position ) ;

				assert( first_allele.size() == 1 && second_allele.size() == 1 ) ;
				aStream.put( first_allele[0] ) ;
				aStream.put( second_allele[0] ) ;
			} else {
				assert(0) ;
			}
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

		void ignore_snp_probability_data(
			std::istream& aStream,
			uint32_t const flags,
			uint32_t number_of_samples
		) {
			if( flags & bgen::e_CompressedSNPBlocks ) {
				uint32_t compressed_data_size ;
				genfile::read_little_endian_integer( aStream, &compressed_data_size ) ;
				aStream.ignore( compressed_data_size ) ;
			}
			else {
				aStream.ignore( 6 * number_of_samples ) ;
			}
		}
	}
}
