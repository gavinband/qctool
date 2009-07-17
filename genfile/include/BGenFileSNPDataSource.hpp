#ifndef BGENFILESNPDATAPROVIDER_HPP
#define BGENFILESNPDATAPROVIDER_HPP

#include <iostream>
#include <string>
#include "snp_data_utils.hpp"
#include "SNPDataSource.hpp"
#include "bgen.hpp"

namespace genfile {
	// This class represents a SNPDataSource which reads its data
	// from a BGEN file.
	class BGenFileSNPDataSource: public IdentifyingDataCachingSNPDataSource
	{
	public:
		BGenFileSNPDataSource( std::string const& filename )
			: m_filename( filename )
		{
			setup( filename, get_compression_type_indicated_by_filename( filename )) ;
		}

		BGenFileSNPDataSource( std::string const& filename, CompressionType compression_type )
			: m_filename( filename )
		{
			setup( filename, compression_type ) ;
		}

		unsigned int number_of_samples() const { return m_number_of_samples ; }
		unsigned int total_number_of_snps() const { return m_total_number_of_snps ; }
		operator bool() const { return *m_stream_ptr ; }

		std::istream& stream() { return *m_stream_ptr ; }
		std::istream const& stream() const { return *m_stream_ptr ; }

	private:

		void read_snp_identifying_data_impl( 
			uint32_t* number_of_samples,
			std::string* SNPID,
			std::string* RSID,
			uint32_t* SNP_position,
			char* allele1,
			char* allele2
		) {
			bgen::impl::read_snp_identifying_data( stream(), number_of_samples, SNPID, RSID, SNP_position, allele1, allele2 ) ;
		}

		void read_snp_probability_data_impl(
			uint32_t* number_of_samples,
			GenotypeProbabilitySetter const& set_genotype_probabilities
		) {
			if( m_flags & bgen::e_CompressedSNPBlocks ) {
				bgen::impl::read_compressed_snp_probability_data( stream(), *number_of_samples, set_genotype_probabilities ) ;
			}
			else {
				bgen::impl::read_snp_probability_data( stream(), *number_of_samples, set_genotype_probabilities ) ;
			}
		}

		void ignore_snp_probability_data_impl(
			uint32_t number_of_samples
		) {
			if( m_flags & bgen::e_CompressedSNPBlocks ) {
				bgen::impl::read_compressed_snp_probability_data( stream(), number_of_samples, ignore() ) ;
			}
			else {
				bgen::impl::read_snp_probability_data( stream(), number_of_samples, ignore() ) ;
			}
		}

	private:

		std::string m_filename ;
		unsigned int m_number_of_samples, m_total_number_of_snps ;
		std::auto_ptr< std::istream > m_stream_ptr ;
		bgen::uint32_t m_flags ;

		void setup( std::string const& filename, CompressionType compression_type ) {
			m_stream_ptr = open_binary_file_for_input( filename, compression_type ) ;
			bgen::uint32_t offset ;
			bgen::read_offset( (*m_stream_ptr), &offset ) ;
			bgen::uint32_t header_size = read_header_data() ;

			if( offset < header_size ) {
				throw FileStructureInvalidError() ;
			}
			// skip any remaining bytes before the first snp block
			m_stream_ptr->ignore( offset - header_size ) ;
		}
	
		bgen::uint32_t read_header_data() {
			bgen::uint32_t header_size ;

			bgen::read_header_block(
				(*m_stream_ptr),
				set_value( header_size ),
				set_value( m_total_number_of_snps ),
				set_value( m_number_of_samples ),
				ignore(),
				set_value( m_flags )
			) ;

			return header_size ;
		}
	} ;
}	

#endif