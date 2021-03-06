
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
#include "genfile/Error.hpp"
#include "genfile/BGenFileSNPDataSink.hpp"
#include "genfile/ToGP.hpp"

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
		m_number_of_bits( number_of_bits ),
		// assume stored probabilities are accurate to 3dps by default
		m_permitted_rounding_error( 0.0005 ) 
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
		m_number_of_bits( number_of_bits ),
		// assume stored probabilities are accurate to 3dps by default
		m_permitted_rounding_error( 0.0005 ) 
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

	void BasicBGenFileSNPDataSink::set_compression_type( std::string const& compression_type ) {
		if( compression_type == "none" ) {
			m_bgen_context.flags = ( m_bgen_context.flags & ~bgen::e_CompressedSNPBlocks ) | bgen::e_NoCompression ;
		} else if( compression_type == "zlib" ) {
			m_bgen_context.flags = ( m_bgen_context.flags & ~bgen::e_CompressedSNPBlocks ) | bgen::e_ZlibCompression ;
		} else if( compression_type == "zstd" ) {
			m_bgen_context.flags = ( m_bgen_context.flags & ~bgen::e_CompressedSNPBlocks ) | bgen::e_ZstdCompression ;
		} else {
			throw BadArgumentError(
				"genfile::BasicBGenFileSNPDataSink::set_compression_type()",
				"compression_type=\"" + compression_type + "\"",
				"Unrecognised compression type"
			) ;
		}
	}
	
	void BasicBGenFileSNPDataSink::set_permitted_input_rounding_error( double const accuracy ) {
		assert( accuracy >= 0 ) ;
		m_permitted_rounding_error = accuracy ;
	}

	SNPDataSink::SinkPos BasicBGenFileSNPDataSink::get_stream_pos() const {
		return SinkPos( this, m_stream_ptr->tellp() ) ;
	}
	
	std::string BasicBGenFileSNPDataSink::get_spec() const { return m_filename ; }

	BasicBGenFileSNPDataSink::operator bool() const { return m_stream_ptr->good() ; }

	namespace {
		std::string get_allele( VariantIdentifyingData const* id_data, std::size_t i ) {
			return id_data->get_allele( i ) ;
		}
	}

	void BasicBGenFileSNPDataSink::write_variant_data_impl(
		VariantIdentifyingData const& id_data,
		VariantDataReader& data_reader,
		Info const& info
	) {
		// std::cerr << id_data << ".\n" ;
		assert( m_have_written_header ) ;
		{
			std::string const& SNPID = ( id_data.number_of_identifiers() > 1 ? id_data.get_identifiers_as_string(",", 1) : id_data.get_identifiers_as_string(",", 0,1) ) ;
			std::string chromosome ;
			if( !id_data.get_position().chromosome().is_missing() ) {
				chromosome = id_data.get_position().chromosome() ;
			}
			byte_t const* const end = bgen::write_snp_identifying_data(
				&m_buffer1,
				m_bgen_context,
				//id_data.get_identifiers_as_string(",", 1),
				SNPID,
				id_data.get_primary_id(),
				chromosome,
				id_data.get_position().position(),
				id_data.number_of_alleles(),
				boost::bind( &get_allele, &id_data, _1 )
			) ;
			stream_ptr()->write( reinterpret_cast< char* >( &(m_buffer1[0]) ), end - &(m_buffer1[0]) ) ;
		}

		{
			bgen::GenotypeDataBlockWriter writer(
				&m_buffer1, &m_buffer2,
				m_bgen_context,
				m_number_of_bits,
				m_permitted_rounding_error
			) ;
			
			data_reader.get( ":genotypes:", to_GP( writer ) ) ;
			
			stream_ptr()->write( reinterpret_cast< char const* >( writer.repr().first ), writer.repr().second  - writer.repr().first ) ;
		}
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

	void BasicBGenFileSNPDataSink::set_free_data( std::string const& free_data ) {
		m_bgen_context.free_data = free_data ;
	}

	void BasicBGenFileSNPDataSink::set_write_sample_identifier_block( bool write ) {
		uint32_t const layout = m_bgen_context.flags & bgen::e_Layout ;
		if( write && (layout == bgen::e_Layout2) ) {
			m_bgen_context.flags |= bgen::e_SampleIdentifiers ;
		} else {
			m_bgen_context.flags &= ~bgen::e_SampleIdentifiers ;
		}
	}
	
	void BasicBGenFileSNPDataSink::setup() {
		m_offset = m_bgen_context.header_size() ;
		update_offset_and_header_block() ;
		// Probably don't need this here.  This will be overwritten in set_sample_names_impl().
	}

	void BasicBGenFileSNPDataSink::set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter sample_name_getter ) {
		assert( sample_name_getter ) ;
		assert( !m_have_written_header ) ;
		// Seems dumb but right now we have to write the header block, then the sample identifier block,
		// then the header block again (to get the offset right).
		// The header is rewritten when this class is deconstructed, but doing this here means at
		// least that the header is well-formed even if an exception occurs.
		std::vector< std::string > sample_ids ;
		for( std::size_t i = 0; i < number_of_samples; ++i ) {
			sample_ids.push_back( sample_name_getter(i).as< std::string >() ) ;
		}
		m_bgen_context.number_of_samples = number_of_samples ;
		m_offset = m_bgen_context.header_size() ;
		update_offset_and_header_block() ;
		if( m_bgen_context.flags & bgen::e_SampleIdentifiers ) {
			m_offset += bgen::write_sample_identifier_block(
				*m_stream_ptr,
				m_bgen_context,
				sample_ids
			) ;
		}
		update_offset_and_header_block() ;
		m_stream_ptr->seekp( m_offset+4, std::ios_base::beg ) ;

		m_have_written_header = true ;
	}

	std::string BasicBGenFileSNPDataSink::serialise( Metadata const& metadata ) const {
		return "" ;
	}

	namespace {
		uint32_t get_flags( std::string const& version ) {
			if( version == "v11" ) {
				return bgen::e_ZlibCompression | bgen::e_Layout1 ;
			} else if( version == "v12" ) {
				return bgen::e_ZlibCompression | bgen::e_Layout2 ;
			} else {
				assert(0) ;
			}
		}
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
			get_flags( version ),
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

