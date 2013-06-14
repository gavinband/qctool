
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <string>
#include <Eigen/Core>
#include "genfile/snp_data_utils.hpp"
#include "genfile/gen.hpp"
#include "genfile/SNPDataSink.hpp"
#include "genfile/GenLikeSNPDataSink.hpp"
#include "genfile/vcf/get_set.hpp"
#include "genfile/vcf/get_set_eigen.hpp"

namespace genfile {
	
	GenLikeSNPDataSink::GenLikeSNPDataSink( std::string const& filename, Chromosome chromosome, CompressionType compression_type ):
		m_filename( filename ),
		m_write_chromosome_column( true ),
		m_compression_type( compression_type )
	{
		setup( filename, m_compression_type ) ;
	}

	GenLikeSNPDataSink::SinkPos GenLikeSNPDataSink::get_stream_pos() const {
		if( m_compression_type == CompressionType( "no_compression" )) {
			return SinkPos( this, m_stream_ptr->tellp() ) ;
		}
		else {
			throw OperationUnsupportedError( "genfile::GenLikeSNPDataSink::get_stream_pos()", "get stream position", "GenLikeSnpDataSink( \"" + m_filename + "\" )" ) ;
		}
	}

	void GenLikeSNPDataSink::omit_chromosome() { m_write_chromosome_column = false ; }

	std::string GenLikeSNPDataSink::get_spec() const { return m_filename ; }

	GenLikeSNPDataSink::operator bool() const { return *m_stream_ptr ; }
	
	void GenLikeSNPDataSink::write_variant( std::ostream& out, genfile::SNPIdentifyingData const& variant ) {
		if( m_write_chromosome_column ) {
			out << variant.get_position().chromosome() << " " ;
		}
		else {
			if( number_of_snps_written() == 0 ) {
				m_chromosome = variant.get_position().chromosome() ;
			}
			else if( m_chromosome != variant.get_position().chromosome() ) {
				throw BadArgumentError( "GenLikeSNPDataSink::write_variant()", "chromosome=" + string_utils::to_string( variant.get_position().chromosome() ) ) ;
			}
		}
		out << variant.get_SNPID() << " "
			<< variant.get_rsid() << " "
			<< variant.get_position().position() << " "
			<< variant.get_first_allele() << " "
			<< variant.get_second_allele() ;
	}

	void GenLikeSNPDataSink::setup( std::string const& filename, CompressionType compression_type ) {
		m_stream_ptr = open_text_file_for_output( filename, compression_type ) ;
		if( !(*m_stream_ptr)) {
			throw ResourceNotOpenedError( m_filename ) ;
		}
		assert( *m_stream_ptr ) ;
	}
}

