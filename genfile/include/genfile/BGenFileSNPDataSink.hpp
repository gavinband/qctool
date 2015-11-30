
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef BGENFILESNPDATASINK_HPP
#define BGENFILESNPDATASINK_HPP

#include <iostream>
#include <string>
#include "genfile/snp_data_utils.hpp"
#include "genfile/SNPDataSink.hpp"
#include "genfile/bgen/bgen.hpp"
#include "genfile/Error.hpp"

namespace genfile {

	// This class encapsulates the basic method of writing a BGen file.
	class BasicBGenFileSNPDataSink: public SNPDataSink
	{
	protected:
		// This class is intended to be used via a derived class.
		BasicBGenFileSNPDataSink(
			std::string const& filename,
			Metadata const& metadata,
			CompressionType compression_type,
			bgen::uint32_t flags,
			int const number_of_bits = 16
		) ;

		BasicBGenFileSNPDataSink(
			std::auto_ptr< std::ostream > stream_ptr,
			std::string const& filename,
			Metadata const& metadata,
			CompressionType compression_type,
			bgen::uint32_t flags,
			int const number_of_bits = 16
		) ;


		SinkPos get_stream_pos() const ;
		std::string get_spec() const ;

	public:
		// Methods required by SNPDataSink
		operator bool() const ;
		
		void write_snp_impl(
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
		) ;

#if 0		
		void write_variant_data_impl(
			SNPIdentifyingData const& id_data,
			VariantDataReader& data_reader,
			Info const& info
		) ;
#endif
		void set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter sample_name_getter ) ;

		void set_number_of_bits( int const bits ) ;
		void set_free_data( std::string const& free_data ) ;
		void set_write_sample_identifier_block( bool write ) ;

		bgen::Context const& bgen_context() const { return m_bgen_context ; }

	protected:
		// Other methods.
		std::auto_ptr< std::ostream >& stream_ptr() ;
		std::string const& filename() const ;
		
		void update_offset_and_header_block() ;
		
	private:

		void setup() ;
		std::string serialise( Metadata const& metadata ) const ;

		std::string m_filename ;
		Metadata m_metadata ;
		bgen::Context m_bgen_context ;
		std::size_t m_offset ;
		std::auto_ptr< std::ostream > m_stream_ptr ;
		bool m_have_written_header ;
		int m_number_of_bits ;
		
		std::vector< byte_t > m_buffer ;
		std::vector< byte_t > m_compression_buffer ;
	} ;


	// This class represents a SNPDataSink which writes its data
	// to an unzipped BGEN file.
	class BGenFileSNPDataSink: public BasicBGenFileSNPDataSink
	{
	public:
		BGenFileSNPDataSink(
			std::string const& filename,
			Metadata const& metadata = Metadata(),
			std::string const& version = "v11",
			int const number_of_bits = 16
		) ;

		BGenFileSNPDataSink(
			std::auto_ptr< std::ostream > stream_ptr,
			std::string const& filename,
			Metadata const& metadata,
			uint32_t const flags,
			int const number_of_bits = 16
		) ;

		BGenFileSNPDataSink(
			std::string const& filename,
			Metadata const& metadata,
			uint32_t const flags,
			int const number_of_bits = 16
		) ;

		~BGenFileSNPDataSink() ;

		friend class SortingBGenFileSNPDataSink ;
	} ;
}

#endif

