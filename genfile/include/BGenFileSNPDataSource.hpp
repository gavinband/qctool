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
			setup( filename, genfile::filename_indicates_file_is_gzipped( filename )) ;
		}

		BGenFileSNPDataSource( std::string const& filename, bool file_is_gzipped )
			: m_filename( filename )
		{
			setup( filename, file_is_gzipped ) ;
		}

		unsigned int number_of_samples() const { return m_number_of_samples ; }
		unsigned int total_number_of_snps() const { return m_total_number_of_snps ; }
		FormatType format() const { return e_BGenFormat ; }
		std::istream& stream() { return *m_stream_ptr ; }
		std::istream const& stream() const { return *m_stream_ptr ; }

	private:

		std::string m_filename ;
		unsigned int m_number_of_samples, m_total_number_of_snps ;
		std::auto_ptr< std::istream > m_stream_ptr ;

		void setup( std::string const& filename, bool file_is_gzipped ) {
			m_stream_ptr = genfile::open_binary_file_for_input( filename, file_is_gzipped ) ;
			genfile::bgen::uint32_t offset ;
			genfile::bgen::read_offset( (*m_stream_ptr), &offset ) ;
			genfile::bgen::uint32_t header_size = read_header_data() ;

			if( offset < header_size ) {
				throw FileStructureInvalidError() ;
			}
			// skip any remaining bytes before the first snp block
			m_stream_ptr->ignore( offset - header_size ) ;
		}
	
		genfile::bgen::uint32_t read_header_data() {
			genfile::bgen::uint32_t header_size ;
		
			genfile::bgen::read_header_block(
				(*m_stream_ptr),
				genfile::set_value( header_size ),
				genfile::set_value( m_total_number_of_snps ),
				genfile::set_value( m_number_of_samples ),
				genfile::ignore(),
				genfile::ignore(),
				genfile::ignore()
			) ;

			return header_size ;
		}
	} ;
}	

#endif