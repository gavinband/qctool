#ifndef GENFILESNPDATASINK_HPP
#define GENFILESNPDATASINK_HPP

#include <iostream>
#include <string>
#include "snp_data_utils.hpp"
#include "gen.hpp"
#include "SNPDataSink.hpp"
#include "SNPDataSource.hpp"

namespace genfile {
	// This class represents a SNPDataSink which writes its data
	// to a plain GEN file.
	class GenFileSNPDataSink: public SNPDataSink
	{
	public:
		GenFileSNPDataSink( std::string const& filename )
			: m_filename( filename )
		{
			setup( filename, get_compression_type_indicated_by_filename( filename )) ;
		}

		GenFileSNPDataSink( std::string const& filename, CompressionType compression_type )
			: m_filename( filename )
		{
			setup( filename, compression_type ) ;
		}

		FormatType format() const { return e_GenFormat ; }
		std::ostream& stream() { return *m_stream_ptr ; }
		std::ostream const& stream() const { return *m_stream_ptr ; }

	private:
		std::string m_filename ;
		std::auto_ptr< std::ostream > m_stream_ptr ;

		void setup( std::string const& filename, CompressionType compression_type ) {
			m_stream_ptr = open_text_file_for_output( filename, compression_type ) ;
			if( !(*m_stream_ptr)) {
				throw genfile::FileNotOpenedError() ;
			}
			assert( *m_stream_ptr ) ;
		}
	} ;
}

#endif