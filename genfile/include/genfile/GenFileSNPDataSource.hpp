#ifndef GENFILESNPDATAPROVIDER_HPP
#define GENFILESNPDATAPROVIDER_HPP

#include <iostream>
#include <string>
#include "snp_data_utils.hpp"
#include "gen.hpp"
#include "SNPDataSource.hpp"

namespace genfile {
	// This class represents a SNPDataSource which reads its data
	// from a plain GEN file.
	class GenFileSNPDataSource: public IdentifyingDataCachingSNPDataSource
	{
	public:
		GenFileSNPDataSource( std::auto_ptr< std::istream > stream, Chromosome chromosome ) ;
		GenFileSNPDataSource( std::string const& filename, Chromosome chromosome ) ;
		GenFileSNPDataSource( std::string const& filename, Chromosome chromosome, CompressionType compression_type ) ;

		unsigned int number_of_samples() const { return m_number_of_samples ; }
		unsigned int total_number_of_snps() const { return m_total_number_of_snps ; }
		
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
			char* allele1,
			char* allele2
		) ;

		void read_snp_probability_data_impl(
			GenotypeProbabilitySetter const& set_genotype_probabilities
		) ;
		
		void ignore_snp_probability_data_impl() ;

	private:
		std::string m_filename ;
		CompressionType m_compression_type ;
		unsigned int m_number_of_samples, m_total_number_of_snps ;
		std::auto_ptr< std::istream > m_stream_ptr ;
		Chromosome m_chromosome ;

		void setup( std::string const& filename, CompressionType compression_type ) ;
		void setup( std::auto_ptr< std::istream > stream_ptr ) ;
		void read_header_data() ;
	} ;
}

#endif