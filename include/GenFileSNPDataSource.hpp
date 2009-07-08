#ifndef GENFILESNPDATAPROVIDER_HPP
#define GENFILESNPDATAPROVIDER_HPP

#include <iostream>
#include <string>
#include "snp_data_utils.hpp"
#include "gen.hpp"
#include "SNPDataSource.hpp"

namespace genfile {
	// This class represents a SNPDataSource which reads its data
	// from a plain GEN file.
	class GenFileSNPDataSource: public SNPDataSource
	{
	public:
		GenFileSNPDataSource( std::string const& filename )
			: m_filename( filename )
		{
			setup( filename, filename_indicates_file_is_gzipped( filename )) ;
		}

		GenFileSNPDataSource( std::string const& filename, bool file_is_gzipped )
			: m_filename( filename )
		{
			setup( filename, file_is_gzipped ) ;
		}

		unsigned int number_of_samples() const { return m_number_of_samples ; }
		unsigned int total_number_of_snps() const { return m_total_number_of_snps ; }
		FormatType format() const { return e_GenFormat ; }
		std::istream& stream() { return *m_stream_ptr ; }
		std::istream const& stream() const { return *m_stream_ptr ; }

	private:
		std::string m_filename ;
		unsigned int m_number_of_samples, m_total_number_of_snps ;
		std::auto_ptr< std::istream > m_stream_ptr ;

		void setup( std::string const& filename, bool file_is_gzipped ) {
			m_stream_ptr = open_text_file_for_input( filename, file_is_gzipped ) ;
			if( !(*m_stream_ptr)) {
				throw FileNotOpenedError() ;
			}
			read_header_data() ;
			// That will have consumed the file, so re-open it.
			m_stream_ptr = open_text_file_for_input( filename, file_is_gzipped ) ;
			assert( *m_stream_ptr ) ;
		}

		void read_header_data() {
			gen::read_header_information(
				*m_stream_ptr,
				set_value( m_total_number_of_snps ),
				set_value( m_number_of_samples ),
				ignore()
			) ;
		}
	} ;
}

#endif