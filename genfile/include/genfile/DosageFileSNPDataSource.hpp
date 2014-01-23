
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_DOSAGE_FILE_SNP_DATA_SOURCE_HPP
#define GENFILE_DOSAGE_FILE_SNP_DATA_SOURCE_HPP

#include <iostream>
#include <string>
#include "snp_data_utils.hpp"
#include "gen.hpp"
#include "SNPDataSource.hpp"
#include "string_utils/slice.hpp"
#include "vcf/MetadataParser.hpp"

namespace genfile {
	namespace impl {
		struct DosageFileSNPDataReader ;
	}
	// This class represents a SNPDataSource which reads its data
	// from a dosage file like the ones qctool outputs.
    // These have one entry per sample counting the number of 'B' alleles.
	class DosageFileSNPDataSource: public IdentifyingDataCachingSNPDataSource
	{
		friend struct impl::DosageFileSNPDataReader ;
	public:
		DosageFileSNPDataSource( std::string const& filename, Chromosome chromosome ) ;
		DosageFileSNPDataSource(
			std::string const& filename,
			Chromosome chromosome,
			CompressionType compression_type
		) ;

		unsigned int number_of_samples() const { return m_number_of_samples ; }
		OptionalSnpCount total_number_of_snps() const { return m_total_number_of_snps ; }
		
		operator bool() const { return *m_stream_ptr ; }
		std::istream& stream() { return *m_stream_ptr ; }
		std::istream const& stream() const { return *m_stream_ptr ; }

		Chromosome chromosome() const { return m_chromosome ; }
		std::string get_source_spec() const { return m_filename ; }

	private:
		void reset_to_start_impl() ;
		
		void read_snp_identifying_data_impl( 
			uint32_t*, // number_of_samples is unused.
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
		CompressionType m_compression_type ;
		unsigned int m_number_of_samples ;
		OptionalSnpCount m_total_number_of_snps ;
		std::auto_ptr< std::istream > m_stream_ptr ;
		Chromosome m_chromosome ;
        
        std::string m_line ;
        typedef string_utils::slice slice ;
        std::vector< slice > m_elts ;

    private:
		void setup( std::string const& filename, CompressionType compression_type ) ;
		void setup( std::auto_ptr< std::istream > stream_ptr ) ;
		void read_header_data() ;
	} ;
}

#endif
