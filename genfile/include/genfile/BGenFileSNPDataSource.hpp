
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
		BGenFileSNPDataSource( std::string const& filename ) ;
		BGenFileSNPDataSource( std::string const& filename, CompressionType compression_type ) ;

		Metadata get_metadata() const ;

		unsigned int number_of_samples() const { return m_number_of_samples ; }
		OptionalSnpCount total_number_of_snps() const { return m_total_number_of_snps ; }
		operator bool() const { return *m_stream_ptr ; }

		std::istream& stream() { return *m_stream_ptr ; }
		std::istream const& stream() const { return *m_stream_ptr ; }
		std::string get_source_spec() const { return m_filename ; }
		uint32_t flags() const { return m_flags ; }

	private:

		void reset_to_start_impl() ;

		void read_snp_identifying_data_impl( 
			uint32_t* number_of_samples,
			std::string* SNPID,
			std::string* RSID,
			Chromosome* chromosome,
			uint32_t* SNP_position,
			std::string* allele1,
			std::string* allele2
		) ;

		VariantDataReader::UniquePtr read_variant_data_impl() ;

		void ignore_snp_probability_data_impl() ;

	private:

		std::string m_filename ;
		unsigned int m_number_of_samples, m_total_number_of_snps ;
		std::auto_ptr< std::istream > m_stream_ptr ;
		uint32_t m_flags ;
		std::string m_version ;

		void setup( std::string const& filename, CompressionType compression_type ) ;
		uint32_t read_header_data() ;
		std::vector< char > m_compressed_data_buffer ;
		std::vector< char > m_uncompressed_data_buffer ;
	} ;
}	

#endif
