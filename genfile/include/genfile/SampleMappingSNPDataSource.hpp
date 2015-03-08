
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_SAMPLE_MAPPING_SNP_DATA_SOURCE_HPP
#define GENFILE_SAMPLE_MAPPING_SNP_DATA_SOURCE_HPP

#include <vector>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/CohortIndividualSource.hpp"

namespace genfile {
	namespace impl {
		struct SampleMapping ;
	}

	struct SampleMappingSNPDataSource: public SNPDataSource {
		typedef std::auto_ptr< SampleMappingSNPDataSource > UniquePtr ;
		SampleMappingSNPDataSource(
			CohortIndividualSource const& reference_samples,
			std::string const& reference_sample_column,
			CohortIndividualSource const& samples,
			std::string const& mapped_sample_column,
			SNPDataSource::UniquePtr source
		) ;

		operator bool() const ;
		Metadata get_metadata() const ;
		unsigned int number_of_samples() const ;
		OptionalSnpCount total_number_of_snps() const ;
		std::string get_source_spec() const ;
		SNPDataSource const& get_parent_source() const ;

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
			std::auto_ptr< impl::SampleMapping > m_sample_mapping ;
			SNPDataSource::UniquePtr m_source ;
	} ;
}

#endif
