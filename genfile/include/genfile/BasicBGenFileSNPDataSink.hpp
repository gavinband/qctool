#ifndef GENFILE_BASICBGENFILESNPDATASINK_HPP
#define GENFILE_BASICBGENFILESNPDATASINK_HPP

#include <iostream>
#include <string>
#include "snp_data_utils.hpp"
#include "SNPDataSink.hpp"
#include "bgen.hpp"

namespace genfile {

	// This class encapsulates the basic method of writing a BGen file.
	class BasicBGenFileSNPDataSink: public SNPDataSink
	{
	protected:
		// This class is intended to be used via a derived class.
		BasicBGenFileSNPDataSink(
			std::string const& filename,
			std::string const& free_data,
			CompressionType compression_type,
			bgen::uint32_t flags
		)
		: 	m_filename( filename ),
			m_free_data( free_data ),
			m_flags( flags )
		{
			setup( filename, compression_type ) ;
		}

	public:
		// Methods required by SNPDataSink
		operator bool() const { return *m_stream_ptr ; }
		
		void write_snp_impl(
			uint32_t number_of_samples,
			std::string SNPID,
			std::string RSID,
			Chromosome chromosome,
			uint32_t SNP_position,
			char first_allele,
			char second_allele,
			GenotypeProbabilityGetter const& get_AA_probability,
			GenotypeProbabilityGetter const& get_AB_probability,
			GenotypeProbabilityGetter const& get_BB_probability
		) {
			std::size_t id_field_size = std::min( std::max( SNPID.size(), RSID.size() ), static_cast< std::size_t >( 255 )) ;
			if( m_flags & bgen::e_CompressedSNPBlocks ) {
				bgen::write_compressed_snp_block( *stream_ptr(), number_of_samples, id_field_size, SNPID, RSID, chromosome, SNP_position, first_allele, second_allele, get_AA_probability, get_AB_probability, get_BB_probability ) ;				
			}
			else {
				bgen::write_snp_block( *stream_ptr(), number_of_samples, id_field_size, SNPID, RSID, chromosome, SNP_position, first_allele, second_allele, get_AA_probability, get_AB_probability, get_BB_probability ) ;
			}
		}

	protected:
		// Other methods.
		std::auto_ptr< std::ostream >& stream_ptr() { return m_stream_ptr ; }
		std::string const& filename() const { return m_filename ; }
		
		void write_header_data( std::ostream& stream ) {
			bgen::write_header_block(
				stream,
				number_of_snps_written(),
				number_of_samples(),
				m_free_data,
				m_flags
			) ;
		}
		
	private:

		void setup( std::string const& filename, CompressionType compression_type ) {
			m_stream_ptr = open_binary_file_for_output( filename, compression_type ) ;
			bgen::uint32_t offset = bgen::get_header_block_size( m_free_data ) ;
			bgen::write_offset( (*m_stream_ptr), offset ) ;
			write_header_data( *m_stream_ptr ) ;
		}

		std::string m_filename ;
		std::string m_free_data ;
		std::auto_ptr< std::ostream > m_stream_ptr ;
		bgen::uint32_t m_flags ;
	} ;
}

#endif