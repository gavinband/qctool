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
		FormatType format() const {
			if( m_flags & bgen::e_CompressedSNPBlocks ) {
				return e_BGenCompressedFormat ;
			}
			else {
				return e_BGenFormat ;
			}
		}
		std::ostream& stream() { return *m_stream_ptr ; }
		std::ostream const& stream() const { return *m_stream_ptr ; }

		std::string const& filename() const { return m_filename ; }

	private:

		std::string m_filename ;
		bgen::uint32_t m_number_of_samples, m_number_of_snps_written ;
		std::string m_free_data ;
		std::auto_ptr< std::ostream > m_stream_ptr ;
		bgen::uint32_t m_flags ;
		
		void setup( std::string const& filename, CompressionType compression_type ) {
			m_stream_ptr = open_binary_file_for_output( filename, compression_type ) ;
			bgen::uint32_t offset = bgen::get_header_block_size( m_free_data ) ;
			bgen::write_offset( (*m_stream_ptr), offset ) ;
			write_header_data( *m_stream_ptr ) ;

			// we will count the number of snps written and so 
			m_number_of_snps_written = 0 ;
		}

	protected:
		
		std::auto_ptr< std::ostream >& stream_ptr() { return m_stream_ptr ; }
		
		void write_header_data( std::ostream& stream ) {
			bgen::write_header_block(
				stream,
				number_of_snps_written(),
				number_of_samples(),
				m_free_data,
				m_flags
			) ;
		}
	} ;


	// This class represents a SNPDataSink which writes its data
	// to an unzipped BGEN file.
	class BGenFileSNPDataSink: public BasicBGenFileSNPDataSink
	{
	public:
		BGenFileSNPDataSink(
			std::string const& filename,
			std::string const& free_data,
			bgen::uint32_t flags
		)
		: 	BasicBGenFileSNPDataSink( filename, free_data, e_NoCompression, flags )
		{
		}

		~BGenFileSNPDataSink() {
			// We are about to close the file.
			// To write the correct header info, we seek back to the start and rewrite the header block
			// The header comes after the offset which is 4 bytes.
			stream().seekp(4, std::ios_base::beg ) ;
			if( stream().bad() ) {
				throw FormatUnsupportedError() ;
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
			std::size_t const buffer_size = 500000
		)
		: 	BasicBGenFileSNPDataSink(
				create_temporary_filename(),
				free_data,
				e_NoCompression,
				bgen::e_NoFlags
			),
			m_filename( filename ),
			m_buffer_size( buffer_size )
		{
			std::cout << "ZippedBGenFileSNPDataSink(): filename = \"" << m_filename << "\", temp filename = \"" << temp_filename() << "\".\n" ;
		}
	
		std::string const& filename() const { return m_filename ; }
		std::string const& temp_filename() const { return BasicBGenFileSNPDataSink::filename() ; }

		~ZippedBGenFileSNPDataSink() {
			// Close the file we were writing.
			stream_ptr().reset() ;
			// Copy the temp file created by our base class, zipping it as we go.
			std::auto_ptr< std::istream > input_file_ptr = open_binary_file_for_input( temp_filename(), e_NoCompression ) ;
			std::auto_ptr< std::ostream > output_file_ptr = open_binary_file_for_output( m_filename, e_GzipCompression ) ;

			bgen::uint32_t offset ;
			bgen::read_offset( *input_file_ptr, &offset ) ;
			bgen::read_header_block( *input_file_ptr, ignore(), ignore(), ignore(), ignore(), ignore()) ;

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