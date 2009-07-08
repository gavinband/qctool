#ifndef BGENFILESNPDATASINK_HPP
#define BGENFILESNPDATASINK_HPP

#include <iostream>
#include <string>
#include "snp_data_utils.hpp"
#include "SNPDataSink.hpp"
#include "bgen.hpp"

namespace genfile {

	// This class encapsulates the basic method of writing a BGen file.
	// However, to get the header block correct, it's necessary to rewrite the header
	// block after all of the snp blocks have been written.  The mechanism to do this
	// depends on whether the file is compressed or not (because we can't seek in a compressed file).
	// Therefore we provide two specialisations below which handle the uncompressed (easy) or compressed (hard) case.
	class BasicBGenFileSNPDataSink: public SNPDataSink
	{
	protected:
		BasicBGenFileSNPDataSink(
			std::string const& filename,
			std::string const& free_data,
			bool file_is_gzipped
		)
		: 	m_filename( filename ),
			m_free_data( free_data )
		{
			setup( filename, file_is_gzipped ) ;
		}

	public:
		FormatType format() const { return e_BGenFormat ; }
		std::ostream& stream() { return *m_stream_ptr ; }
		std::ostream const& stream() const { return *m_stream_ptr ; }

		std::string const& filename() const { return m_filename ; }

	private:

		std::string m_filename ;
		genfile::bgen::uint32_t m_number_of_samples, m_number_of_snps_written, m_snp_block_size ;
		std::string m_free_data ;
		std::auto_ptr< std::ostream > m_stream_ptr ;

		void setup( std::string const& filename, bool file_is_gzipped ) {
			m_stream_ptr = genfile::open_binary_file_for_output( filename, file_is_gzipped ) ;
			genfile::bgen::uint32_t offset = genfile::bgen::get_header_block_size( m_free_data ) ;
			bgen::write_offset( (*m_stream_ptr), offset ) ;
			write_header_data( *m_stream_ptr ) ;

			// we will count the number of snps written and so 
			m_number_of_snps_written = 0 ;
		}

	protected:
		void write_header_data( std::ostream& stream ) {
			bgen::write_header_block(
				stream,
				number_of_snps_written(),
				number_of_samples(),
				m_snp_block_size,
				m_free_data
			) ;
		}
	} ;


	// This class represents a SNPDataSink which writes its data
	// to an unzipped BGEN file.
	class UnzippedBGenFileSNPDataSink: public BasicBGenFileSNPDataSink
	{
	public:
		UnzippedBGenFileSNPDataSink(
			std::string const& filename,
			std::string const& free_data
		)
		: 	BasicBGenFileSNPDataSink( filename, free_data, false )
		{
		}

		~UnzippedBGenFileSNPDataSink() {
			// We are about to close the file.
			// To write the correct header info, we seek back to the start and rewrite the header block
			stream().seekp(0, std::ios_base::beg ) ;
			if( stream().bad() ) {
				throw genfile::FormatUnsupportedError() ;
			}

			write_header_data( stream() ) ;
		}
	} ;



	// This class represents a SNPDataSink which writes its data
	// to an unzipped BGEN file.
	class ZippedBGenFileSNPDataSink: public BasicBGenFileSNPDataSink
	{
	public:
		ZippedBGenFileSNPDataSink(
			std::string const& filename,
			std::string const& free_data,
			std::size_t const buffer_size = 2000
		)
		: 	BasicBGenFileSNPDataSink(
				genfile::create_temporary_filename(),
				free_data,
				true
			),
			m_filename( filename ),
			m_buffer_size( buffer_size )
		{
			std::cout << "ZippedBGenFileSNPDataSink(): filename = \"" << m_filename << "\", temp filename = \"" << temp_filename() << "\".\n" ;
		}
	
		std::string const& filename() const { return m_filename ; }
		std::string const& temp_filename() const { return BasicBGenFileSNPDataSink::filename() ; }

		~ZippedBGenFileSNPDataSink() {
			// Copy the temp file created by our base class, zipping it as we go.
			std::auto_ptr< std::istream > input_file_ptr = genfile::open_binary_file_for_input( temp_filename(), false ) ;
			std::auto_ptr< std::ostream > output_file_ptr = genfile::open_binary_file_for_output( m_filename, true ) ;

			bgen::uint32_t offset ;
			bgen::read_offset( *input_file_ptr, &offset ) ;
			bgen::read_header_block( *input_file_ptr, genfile::ignore(), genfile::ignore(), genfile::ignore(), genfile::ignore(), genfile::ignore(), genfile::ignore()) ;

			// Write the offset and header into the output file.
			bgen::write_offset( *output_file_ptr, offset ) ;
			BasicBGenFileSNPDataSink::write_header_data( *output_file_ptr ) ;

			// Copy the rest of the file contents.
			std::vector< char > buffer( m_buffer_size ) ;
			do {
				input_file_ptr->read( &buffer[0], m_buffer_size ) ;
				output_file_ptr->write( &buffer[0], input_file_ptr->gcount() ) ;
			}
			while( *input_file_ptr ) ;
		}
	
	private:
	
		std::string m_filename ;
		std::size_t const m_buffer_size ;
	} ;
}

#endif