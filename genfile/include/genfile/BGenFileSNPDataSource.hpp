
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef BGENFILESNPDATAPROVIDER_HPP
#define BGENFILESNPDATAPROVIDER_HPP

#include <iostream>
#include <string>
#include "snp_data_utils.hpp"
#include "SNPDataSource.hpp"
#include "bgen/bgen.hpp"
#include "Chromosome.hpp"

namespace genfile {
	namespace impl {
		struct BGenFileSNPDataReader ;
	}

	// This class represents a SNPDataSource which reads its data
	// from a BGEN file.
	class BGenFileSNPDataSource: public IdentifyingDataCachingSNPDataSource
	{
		friend struct impl::BGenFileSNPDataReader ;
	public:
		BGenFileSNPDataSource( std::string const& filename, Chromosome missing_chromosome = Chromosome() ) ;

		Metadata get_metadata() const ;

		unsigned int number_of_samples() const { return m_bgen_context.number_of_samples ; }
		void get_sample_ids( GetSampleIds ) const ;
		OptionalSnpCount total_number_of_snps() const { return m_bgen_context.number_of_variants ; }
		operator bool() const { return m_stream_ptr->good() ; }

		std::istream& stream() { return *m_stream_ptr ; }
		std::istream const& stream() const { return *m_stream_ptr ; }
		std::string get_source_spec() const { return m_filename ; }
		bgen::Context const& bgen_context() const { return m_bgen_context ; }

	private:

		void reset_to_start_impl() ;

		void read_snp_identifying_data_impl( VariantIdentifyingData* result ) ;

		VariantDataReader::UniquePtr read_variant_data_impl() ;

		void ignore_snp_probability_data_impl() ;

	private:

		std::string m_filename ;
		Chromosome m_missing_chromosome ;
		bgen::Context m_bgen_context ;
		boost::optional< std::vector< std::string > > m_sample_ids ;
		std::auto_ptr< std::istream > m_stream_ptr ;

		void setup( std::string const& filename, CompressionType compression_type ) ;
		uint32_t read_header_data() ;
		std::vector< byte_t > m_compressed_data_buffer ;
		std::vector< byte_t > m_uncompressed_data_buffer ;
		
		typedef std::istream_iterator<char> StreamIterator ;
	} ;
}	

#endif
