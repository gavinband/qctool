#ifndef GENFILE_VCF_FORMAT_SNP_DATA_SOURCE_HPP
#define GENFILE_VCF_FORMAT_SNP_DATA_SOURCE_HPP

#include <boost/ptr_container/ptr_map.hpp>
#include "genfile/VariantEntry.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/vcf/MetadataParser.hpp"
#include "genfile/vcf/Types.hpp"
#include "genfile/VariantDataReader.hpp"

namespace genfile {
	class VCFFormatSNPDataSource: public SNPDataSource
	// SNPDataSource which obtains data from a file in VCF format.
	// I used the 4.1 spec, available here:
	// http://www.1000genomes.org/wiki/Analysis/Variant Call Format/vcf-variant-call-format-version-41
	// to implement this.
	{
	public:
		VCFFormatSNPDataSource(
			std::auto_ptr< std::istream > stream_ptr,
			std::string const& genotype_probability_field
		) ;
		VCFFormatSNPDataSource(
			std::string const& filename,
			std::string const& genotype_probability_field
		) ;
		VCFFormatSNPDataSource(
			std::string const& filename,
			std::string const& index_filename,
			std::string const& genotype_probability_field
		) ;
	public:
		typedef vcf::MetadataParser::Metadata Metadata ;
		void update_metadata( Metadata const& metadata ) ;
		
		operator bool() const ;
		unsigned int number_of_samples() const ;
		unsigned int total_number_of_snps() const ;
		std::string get_source_spec() const ;
		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const ;
		void set_genotype_probability_field( std::string const& ) ;

		std::size_t get_index_of_first_data_line() const { return m_metadata_parser.get_number_of_lines() + 1 ; }
		std::size_t get_index_of_first_data_column() const { return 9 ; }
	protected:

		void get_snp_identifying_data_impl( 
			IntegerSetter const& set_number_of_samples,
			StringSetter const& set_SNPID,
			StringSetter const& set_RSID,
			ChromosomeSetter const& set_chromosome,
			SNPPositionSetter const& set_SNP_position,
			AlleleSetter const& set_allele1,
			AlleleSetter const& set_allele2
		) ;	

		void read_snp_probability_data_impl(
			GenotypeProbabilitySetter const& set_genotype_probabilities
		) ;
		
		VariantDataReader::UniquePtr read_variant_data_impl() ;

		void ignore_snp_probability_data_impl() ;
		void reset_to_start_impl() ;

	private:
		std::string const m_spec ;
		CompressionType m_compression_type ;
		std::auto_ptr< std::istream > m_stream_ptr ;
		vcf::MetadataParser m_metadata_parser ;
		vcf::MetadataParser::Metadata m_metadata  ;
		typedef boost::ptr_map< std::string, vcf::VCFEntryType > EntryTypeMap ;
		EntryTypeMap m_info_types ;
		EntryTypeMap m_format_types ;
		std::string m_genotype_probability_field ;

		std::vector< std::string > const m_column_names ;
		std::size_t const m_number_of_samples ;
		std::size_t const m_number_of_lines ;
		
		// We record the alleles per SNP, so that they can be used on subsequent SNPs.
		std::vector< std::string > m_variant_alleles ; 
	private:
		void setup() ;
		void check_genotype_probability_field( std::string const& field ) const ;
		std::vector< std::string > read_column_names( std::istream& stream ) const ;
		char read_format_and_get_trailing_char( std::string& format, std::size_t column ) const ;
		std::size_t determine_number_of_lines(
			std::istream& vcf_file_stream,
			vcf::MetadataParser::Metadata const& metadata,
			std::auto_ptr< std::istream > index_file = std::auto_ptr< std::istream >()
		) const ;
		std::size_t count_lines( std::istream& ) const ;
		void reset_stream() ;
		void read_element( std::string& elt, char delim, std::size_t column ) const ;
	private:
		VCFFormatSNPDataSource( VCFFormatSNPDataSource const& other ) ;
	} ;
	
	
}

#endif
