#ifndef GENFILE_IMPUTE_HAPLOTYPES_SNP_DATA_SOURCE_HPP
#define GENFILE_IMPUTE_HAPLOTYPES_SNP_DATA_SOURCE_HPP

#include <string>
#include <vector>
#include "genfile/SNPDataSource.hpp"
#include "genfile/string_utils/slice.hpp"

namespace genfile {
	// Read genotypes from a .legend file and a .haps file, in the format that can be downloaded
	// From the IMPUTE2 website http://mathgen.stats.ox.ac.uk/impute/impute_v2.html
	// It is assumed that haplotypes come in pairs, corresponding to haplotypes from the same individual.
	class ImputeHaplotypesSNPDataSource: public SNPDataSource
	{
	public:
		ImputeHaplotypesSNPDataSource( std::string const& haplotypes_filename, Chromosome chromosome ) ;
		ImputeHaplotypesSNPDataSource( std::string const& haplotypes_filename, Chromosome chromosome, CompressionType compression_type ) ;

		unsigned int number_of_samples() const { return m_number_of_samples ; }
		unsigned int total_number_of_snps() const { return m_snps.size() ; }
		
		operator bool() const { return m_good ; }
		std::istream& stream() { return *m_stream_ptr ; }
		std::istream const& stream() const { return *m_stream_ptr ; }

		Chromosome chromosome() const { return m_chromosome ; }
		std::string get_source_spec() const { return m_haplotypes_filename ; }

	private:
		void reset_to_start_impl() ;
		
		void get_snp_identifying_data_impl( 
			IntegerSetter const& set_number_of_samples,
			StringSetter const& set_SNPID,
			StringSetter const& set_RSID,
			ChromosomeSetter const& set_chromosome,
			SNPPositionSetter const& set_SNP_position,
			AlleleSetter const& set_allele1,
			AlleleSetter const& set_allele2
		) ;

		VariantDataReader::UniquePtr read_variant_data_impl() ;

		void read_snp_probability_data_impl(
			GenotypeProbabilitySetter const& set_genotype_probabilities
		) ;
		
		void ignore_snp_probability_data_impl() ;

	private:
		std::string m_legend_filename, m_haplotypes_filename ;
		CompressionType m_compression_type ;
		unsigned int m_number_of_samples ;
		std::vector< SNPIdentifyingData > m_snps ;
		std::auto_ptr< std::istream > m_stream_ptr ;
		Chromosome m_chromosome ;
		bool m_good ;

		void setup( std::string const& filename, CompressionType compression_type ) ;
		void setup( std::auto_ptr< std::istream > stream_ptr ) ;
		void read_header_data() ;
	} ;
}

#endif