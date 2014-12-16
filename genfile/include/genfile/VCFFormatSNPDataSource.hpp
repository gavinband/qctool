
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_VCF_FORMAT_SNP_DATA_SOURCE_HPP
#define GENFILE_VCF_FORMAT_SNP_DATA_SOURCE_HPP

#include <boost/optional.hpp>
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/bimap.hpp>
#include <boost/optional.hpp>
#include "genfile/VariantEntry.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/vcf/MetadataParser.hpp"
#include "genfile/vcf/Types.hpp"
#include "genfile/VariantDataReader.hpp"

namespace genfile {
	
	namespace impl {
		struct VCFFormatDataReader ;
	}
	
	class VCFFormatSNPDataSource: public SNPDataSource
	// SNPDataSource which obtains data from a file in VCF format.
	// I used the 4.1 spec, available here:
	// http://www.1000genomes.org/wiki/Analysis/Variant Call Format/vcf-variant-call-format-version-41
	// to implement this.
	{
		friend struct impl::VCFFormatDataReader ;
	public:
		typedef std::auto_ptr< VCFFormatSNPDataSource > UniquePtr ;
		typedef vcf::MetadataParser::Metadata Metadata ;
		
		VCFFormatSNPDataSource(
			std::auto_ptr< std::istream > stream_ptr,
			boost::optional< Metadata > metadata = boost::optional< Metadata >()
		) ;
		VCFFormatSNPDataSource(
			std::string const& filename,
			boost::optional< Metadata > metadata = boost::optional< Metadata >()
		) ;
	public:
		operator bool() const ;
		Metadata get_metadata() const ;
		unsigned int number_of_samples() const ;
		OptionalSnpCount total_number_of_snps() const ;
		std::string get_source_spec() const ;
		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const ;
		void set_field_mapping( std::string const& key, std::string const& value ) ;

		std::size_t get_index_of_first_data_line() const { return m_metadata_parser->get_number_of_lines() + 1 ; }
		std::size_t get_index_of_first_data_column() const { return 9 ; }

		void set_strict_mode( bool value ) ;

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
		
		void get_snp_identifying_data_impl( 
			std::istream* stream_ptr,
			IntegerSetter const& set_number_of_samples,
			StringSetter const& set_SNPID,
			StringSetter const& set_RSID,
			ChromosomeSetter const& set_chromosome,
			SNPPositionSetter const& set_SNP_position,
			AlleleSetter const& set_allele1,
			AlleleSetter const& set_allele2
		) ;
		

		VariantDataReader::UniquePtr read_variant_data_impl() ;

		void ignore_snp_probability_data_impl() ;
		void reset_to_start_impl() ;

	private:
		std::string const m_spec ;
		CompressionType m_compression_type ;
		std::auto_ptr< std::istream > m_stream_ptr ;
		vcf::MetadataParser::UniquePtr m_metadata_parser ;
		vcf::MetadataParser::Metadata m_metadata  ;
		typedef boost::ptr_map< std::string, vcf::VCFEntryType > EntryTypeMap ;
		EntryTypeMap m_info_types ;
		EntryTypeMap m_format_types ;
		typedef boost::bimap< std::string, std::string > FieldMapping ;
		FieldMapping m_field_mapping ;
		std::string m_genotype_field ;
		std::string m_intensity_field ;
		bool m_have_id_data ;
		bool m_strict_mode ;

		std::vector< std::string > const m_column_names ;
		std::size_t const m_number_of_samples ;
		OptionalSnpCount m_number_of_lines ;
		
		// We record the alleles per SNP, so that they can be used on subsequent SNPs.
		std::vector< std::string > m_variant_alleles ;
		
		std::string m_CHROM ;
		std::string m_POS ;
		std::string m_ID ;
		std::string m_REF ;
		std::string m_ALT ;
		std::string m_QUAL ;
		std::string m_FILTER ;
		std::string m_INFO ;
		
	private:
		void setup() ;
		void check_genotype_probability_field( std::string const& field ) const ;
		std::vector< std::string > read_column_names( std::istream& stream ) const ;
		std::string read_format() ;
		void reset_stream() ;
		void read_element( std::string& elt, char delim, std::size_t column ) const ;
	private:
		VCFFormatSNPDataSource( VCFFormatSNPDataSource const& other ) ;
	} ;
	
	
}

#endif
