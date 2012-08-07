
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <string>
#include "genfile/snp_data_utils.hpp"
#include "genfile/SNPDataSink.hpp"
#include "genfile/bgen.hpp"
#include "genfile/Error.hpp"
#include "genfile/BGenFileSNPDataSink.hpp"

namespace genfile {
	// This class is intended to be used via a derived class.
	BasicBGenFileSNPDataSink::BasicBGenFileSNPDataSink(
		std::string const& filename,
		Metadata const& metadata,
		CompressionType compression_type,
		bgen::uint32_t flags
	)
	: 	m_filename( filename ),
		m_metadata( metadata ),
		m_free_data( serialise( m_metadata )),
		m_stream_ptr( open_binary_file_for_output( m_filename, compression_type ) ),
		m_have_written_header( false ),
		m_flags( flags )
		
	{
		setup() ;
	}

	BasicBGenFileSNPDataSink::BasicBGenFileSNPDataSink(
		std::auto_ptr< std::ostream > stream_ptr,
		std::string const& filename,
		Metadata const& metadata,
		CompressionType compression_type,
		bgen::uint32_t flags
	)
	: 	m_filename( filename ),
		m_metadata( metadata ),
		m_free_data( serialise( m_metadata )),
		m_stream_ptr( stream_ptr ),
		m_have_written_header( false ),
		m_flags( flags )
	{
		setup() ;
	}

	std::ostream::streampos BasicBGenFileSNPDataSink::get_stream_pos() const {
		return m_stream_ptr->tellp() ;
	}
	
	std::string BasicBGenFileSNPDataSink::get_spec() const { return m_filename ; }

	BasicBGenFileSNPDataSink::operator bool() const { return *m_stream_ptr ; }
	
	void BasicBGenFileSNPDataSink::write_snp_impl(
		uint32_t number_of_samples,
		std::string SNPID,
		std::string RSID,
		Chromosome chromosome,
		uint32_t SNP_position,
		std::string first_allele,
		std::string second_allele,
		GenotypeProbabilityGetter const& get_AA_probability,
		GenotypeProbabilityGetter const& get_AB_probability,
		GenotypeProbabilityGetter const& get_BB_probability,
		Info const&
	) {
		assert( m_have_written_header ) ;
		std::size_t id_field_size = std::min( std::max( SNPID.size(), RSID.size() ), static_cast< std::size_t >( 255 )) ;
		bgen::write_snp_block( *stream_ptr(), m_flags, number_of_samples, id_field_size, SNPID, RSID, chromosome, SNP_position, first_allele, second_allele, get_AA_probability, get_AB_probability, get_BB_probability ) ;
	}

	std::auto_ptr< std::ostream >& BasicBGenFileSNPDataSink::stream_ptr() { return m_stream_ptr ; }
	std::string const& BasicBGenFileSNPDataSink::filename() const { return m_filename ; }
	
	void BasicBGenFileSNPDataSink::write_header_data( std::ostream& stream, std::size_t const number_of_samples ) {
		bgen::write_header_block(
			stream,
			number_of_snps_written(),
			number_of_samples,
			m_free_data,
			m_flags
		) ;
	}
	
	void BasicBGenFileSNPDataSink::setup() {
		bgen::uint32_t offset = bgen::get_header_block_size( m_free_data ) ;
		bgen::write_offset( (*m_stream_ptr), offset ) ;
	}

	void BasicBGenFileSNPDataSink::set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter sample_name_getter ) {
		assert( sample_name_getter ) ;
		assert( !m_have_written_header ) ;
		write_header_data( *m_stream_ptr, number_of_samples ) ;
		m_have_written_header = true ;
	}

	std::string BasicBGenFileSNPDataSink::serialise( Metadata const& metadata ) const {
		return "" ;
	}

	BGenFileSNPDataSink::BGenFileSNPDataSink(
		std::string const& filename,
		Metadata const& metadata
	)
	: 	BasicBGenFileSNPDataSink( filename, metadata, "no_compression", bgen::e_CompressedSNPBlocks | bgen::e_LongIds )
	{
	}

	BGenFileSNPDataSink::BGenFileSNPDataSink(
		std::auto_ptr< std::ostream > stream_ptr,
		std::string const& filename,
		Metadata const& metadata,
		uint32_t const flags
	)
	: 	BasicBGenFileSNPDataSink( stream_ptr, filename, metadata, "no_compression", flags )
	{
	}

	BGenFileSNPDataSink::BGenFileSNPDataSink(
		std::string const& filename,
		Metadata const& metadata,
		uint32_t const flags
	)
	: 	BasicBGenFileSNPDataSink( filename, metadata, "no_compression", flags )
	{
	}

	BGenFileSNPDataSink::~BGenFileSNPDataSink() {
		// We are about to close the file.
		// To write the correct header info, we seek back to the start and rewrite the header block
		// The header comes after the offset which is 4 bytes.
		stream_ptr()->seekp( 4, std::ios_base::beg ) ;
		if( stream_ptr()->bad() ) {
			throw OutputError( filename() ) ;
		}

		write_header_data( *stream_ptr(), number_of_samples() ) ;
	}
}

