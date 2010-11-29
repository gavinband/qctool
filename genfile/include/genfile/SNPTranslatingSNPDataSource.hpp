#ifndef GENFILE_SNPTRANSLATING_SNP_DATA_SOURCE_HPP
#define GENFILE_SNPTRANSLATING_SNP_DATA_SOURCE_HPP

#include <map>
#include <string>
#include <memory>

#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPIdentifyingData.hpp"

namespace genfile {

	class SNPTranslatingSNPDataSource: public SNPDataSource
	{
	public:
		typedef std::auto_ptr< SNPTranslatingSNPDataSource > UniquePtr ;
		typedef std::map< SNPIdentifyingData, SNPIdentifyingData > Dictionary ;
		
		static UniquePtr create(
			SNPDataSource::UniquePtr source,
			Dictionary const& dict
		) ;

	public:
		SNPTranslatingSNPDataSource(
			SNPDataSource::UniquePtr source,
			Dictionary const& dict
		) ;
	
		operator bool() const {
			return  (*m_source) ;
		}

		unsigned int number_of_samples() const {
			return m_source->number_of_samples() ;
		}

		unsigned int total_number_of_snps() const {
			return m_source->total_number_of_snps() ;
		}

		std::string get_source_spec() const {
			return "snp-translated:" + m_source->get_source_spec() ;
		}

		SNPDataSource const& get_parent_source() const {
			return *m_source ;
		}

		SNPDataSource const& get_base_source() const {
			return m_source->get_base_source() ;
		}

	private:
		SNPDataSource::UniquePtr m_source ;
		Dictionary m_dictionary ;

	private:
		void reset_to_start_impl() {
			m_source->reset_to_start() ;
		}

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
