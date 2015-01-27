
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <string>
#include <boost/bind.hpp>
#include "genfile/snp_data_utils.hpp"
#include "genfile/SNPDataSink.hpp"
#include "genfile/bgen/bgen.hpp"
#include "genfile/bgen/impl.hpp"
#include "genfile/Error.hpp"
#include "genfile/BGenFileSNPDataSink.hpp"

namespace genfile {
	// This class is intended to be used via a derived class.
	BasicBGenFileSNPDataSink::BasicBGenFileSNPDataSink(
		std::string const& filename,
		Metadata const& metadata,
		CompressionType compression_type,
		bgen::uint32_t flags,
		int const number_of_bits
	)
	: 	m_filename( filename ),
		m_metadata( metadata ),
		m_offset(0),
		m_stream_ptr( open_binary_file_for_output( m_filename, compression_type ) ),
		m_have_written_header( false ),
		m_number_of_bits( number_of_bits )
	{
		m_bgen_context.flags = flags ;
		m_bgen_context.free_data = serialise( metadata ) ;
		setup() ;
	}

	BasicBGenFileSNPDataSink::BasicBGenFileSNPDataSink(
		std::auto_ptr< std::ostream > stream_ptr,
		std::string const& filename,
		Metadata const& metadata,
		CompressionType compression_type,
		bgen::uint32_t flags,
		int const number_of_bits
	)
	: 	m_filename( filename ),
		m_metadata( metadata ),
		m_offset(0),
		m_stream_ptr( stream_ptr ),
		m_have_written_header( false ),
		m_number_of_bits( number_of_bits )
	{
		m_bgen_context.flags = flags ;
		m_bgen_context.free_data = serialise( metadata ) ;
		setup() ;
	}

	void BasicBGenFileSNPDataSink::set_number_of_bits( int const bits ) {
		assert( bits > 0 ) ;
		assert( bits <= 32 ) ;
		m_number_of_bits = bits ;
	}

	SNPDataSink::SinkPos BasicBGenFileSNPDataSink::get_stream_pos() const {
		return SinkPos( this, m_stream_ptr->tellp() ) ;
	}
	
	std::string BasicBGenFileSNPDataSink::get_spec() const { return m_filename ; }

	BasicBGenFileSNPDataSink::operator bool() const { return m_stream_ptr->good() ; }
	
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
		bgen::write_snp_identifying_data( *stream_ptr(), m_bgen_context, id_field_size, SNPID, RSID, chromosome, SNP_position, first_allele, second_allele ) ;
		bgen::write_snp_probability_data(
			*m_stream_ptr,
			m_bgen_context, get_AA_probability, get_AB_probability, get_BB_probability, m_number_of_bits,
			&m_buffer,
			&m_compression_buffer
		) ;
	}

	std::auto_ptr< std::ostream >& BasicBGenFileSNPDataSink::stream_ptr() { return m_stream_ptr ; }
	std::string const& BasicBGenFileSNPDataSink::filename() const { return m_filename ; }
	
	void BasicBGenFileSNPDataSink::update_offset_and_header_block() {
		m_bgen_context.number_of_variants = number_of_snps_written() ;
		m_stream_ptr->seekp( 0, std::ios_base::beg ) ;
		if( !stream_ptr()->bad() ) {
			bgen::write_offset( *m_stream_ptr, m_offset ) ;
			bgen::write_header_block(
				*m_stream_ptr,
				m_bgen_context
			) ;
		}
	}
	
	void BasicBGenFileSNPDataSink::setup() {
	}

	void BasicBGenFileSNPDataSink::set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter sample_name_getter ) {
		assert( sample_name_getter ) ;
		assert( !m_have_written_header ) ;
		// Start by writing dummy offset and header block, and the sample identifier block.
		// Offset and header blocks get overwritten with correct values on destruction.
		std::vector< std::string > sample_ids ;
		for( std::size_t i = 0; i < number_of_samples; ++i ) {
			sample_ids.push_back( sample_name_getter(i).as< std::string >() ) ;
		}
		m_bgen_context.number_of_samples = number_of_samples ;
		update_offset_and_header_block() ;
		m_offset = m_bgen_context.header_size() ;
//		m_offset += bgen::write_sample_identifier_block(
//			*m_stream_ptr,
//			m_bgen_context,
//			sample_ids
//		) ;
		m_have_written_header = true ;
	}

	std::string BasicBGenFileSNPDataSink::serialise( Metadata const& metadata ) const {
		return "" ;
	}

	BGenFileSNPDataSink::BGenFileSNPDataSink(
		std::string const& filename,
		Metadata const& metadata,
		std::string const& version,
		int const number_of_bits
	)
		: BasicBGenFileSNPDataSink(
			filename,
			metadata,
			"no_compression",
			bgen::impl::get_flags( version ),
			number_of_bits
		)
	{
	}

	BGenFileSNPDataSink::BGenFileSNPDataSink(
		std::auto_ptr< std::ostream > stream_ptr,
		std::string const& filename,
		Metadata const& metadata,
		uint32_t const flags,
		int const number_of_bits
	)
	: 	BasicBGenFileSNPDataSink( stream_ptr, filename, metadata, "no_compression", flags, number_of_bits )
	{
	}

	BGenFileSNPDataSink::BGenFileSNPDataSink(
		std::string const& filename,
		Metadata const& metadata,
		uint32_t const flags,
		int const number_of_bits
	)
	: 	BasicBGenFileSNPDataSink( filename, metadata, "no_compression", flags, number_of_bits )
	{
	}

	BGenFileSNPDataSink::~BGenFileSNPDataSink() {
		// We are about to close the file.
		// To write the correct header info, we seek back to the start and rewrite the header block
		// The header comes after the offset which is 4 bytes.
		update_offset_and_header_block() ;
	}
}

