#ifndef GENFILE_THRESHOLDING_SNP_DATA_SOURCE_HPP
#define GENFILE_THRESHOLDING_SNP_DATA_SOURCE_HPP

#include "genfile/SNPDataSource.hpp"

namespace genfile {
	struct ThresholdingSNPDataSource: public SNPDataSource {
		ThresholdingSNPDataSource( SNPDataSource::UniquePtr source, double threshhold ) ;
		
		unsigned int number_of_samples() const ;
		unsigned int total_number_of_snps() const ;
		unsigned int number_of_sources() const ;
		unsigned int number_of_snps_in_source( std::size_t source_index ) const ;
		SNPDataSource const& get_source( std::size_t source_index ) const ;
		operator bool() const ;

	private:
		SNPDataSource::UniquePtr m_source ;
		double const m_threshhold ;

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

		void read_snp_probability_data_impl(
			GenotypeProbabilitySetter const& set_genotype_probabilities
		) ;
		
		void ignore_snp_probability_data_impl() ;
	} ;	
}

#endif
