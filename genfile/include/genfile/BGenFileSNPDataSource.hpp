#ifndef BGENFILESNPDATAPROVIDER_HPP
#define BGENFILESNPDATAPROVIDER_HPP

#include <iostream>
#include <string>
#include "snp_data_utils.hpp"
#include "SNPDataSource.hpp"
#include "bgen.hpp"

namespace genfile {
	namespace impl {
		struct BGenFileSNPDataReader ;
	}

	// This class represents a SNPDataSource which reads its data
	// from a BGEN file.
	class BGenFileSNPDataSource: public IdentifyingDataCachingSNPDataSource
	{
		friend class impl::BGenFileSNPDataReader ;
	public:
		BGenFileSNPDataSource( std::string const& filename ) ;
		BGenFileSNPDataSource( std::string const& filename, CompressionType compression_type ) ;

		unsigned int number_of_samples() const { return m_number_of_samples ; }
		unsigned int total_number_of_snps() const { return m_total_number_of_snps ; }
		operator bool() const { return *m_stream_ptr ; }

		std::istream& stream() { return *m_stream_ptr ; }
		std::istream const& stream() const { return *m_stream_ptr ; }
		std::string get_source_spec() const { return m_filename ; }

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

		void read_snp_probability_data_impl(
			GenotypeProbabilitySetter const& set_genotype_probabilities
		) ;

		void ignore_snp_probability_data_impl() ;

	private:

		std::string m_filename ;
		unsigned int m_number_of_samples, m_total_number_of_snps ;
		std::auto_ptr< std::istream > m_stream_ptr ;
		bgen::uint32_t m_flags ;

		void setup( std::string const& filename, CompressionType compression_type ) ;
		bgen::uint32_t read_header_data() ;
	} ;
}	

#endif