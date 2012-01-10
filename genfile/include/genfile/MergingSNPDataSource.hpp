#ifndef GENFILE_MERGING_SNP_DATA_SOURCE_HPP
#define GENFILE_MERGING_SNP_DATA_SOURCE_HPP

#include <vector>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SNPDataSource.hpp"

namespace genfile {
	struct MergingSNPDataSource: public SNPDataSource {

		MergingSNPDataSource() ;
		~MergingSNPDataSource() ;
		
		void add_source( SNPDataSource::UniquePtr ) ;
	
		operator bool() const ;
		// Return the number of samples represented in the snps in this source.
		unsigned int number_of_samples() const ;
		// Return the total number of snps the source contains.
		unsigned int total_number_of_snps() const ;

		// Return a string identifying the source of the SNP data
		std::string get_source_spec() const ;
	
	private:
		void get_snp_identifying_data_impl( 
			IntegerSetter const& set_number_of_samples,
			StringSetter const& set_SNPID,
			StringSetter const& set_RSID,
			ChromosomeSetter const& set_chromosome,
			SNPPositionSetter const& set_SNP_position,
			AlleleSetter const& set_allele1,
			AlleleSetter const& set_allele2
		) ;	

		VariantDataReader::UniquePtr read_variant_data_impl()  ;

		void ignore_snp_probability_data_impl() ;
		void reset_to_start_impl() ;
	
	private:
		std::vector< SNPDataSource* > m_sources ;
		std::multimap< SNPIdentifyingData, std::size_t > m_current_snps ;
		
		void get_top_snp_in_source( std::size_t source_i ) ;
	} ;
}

#endif
