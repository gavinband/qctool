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
	class BGenFileSNPDataSource: public SNPDataSource
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
		FormatType format() const {
			if (m_flags & bgen::e_CompressedSNPBlocks) {
				return e_BGenCompressedFormat ;
			}
			else {
				return e_BGenFormat ;
			}
		}
		std::istream& stream() { return *m_stream_ptr ; }
		std::istream const& stream() const { return *m_stream_ptr ; }

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
				ignore(),
				set_value( m_flags )
			) ;

			return header_size ;
		}
	} ;
}	

#endif