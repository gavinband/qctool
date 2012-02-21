#include "genfile/Error.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPIdentifyingDataTest.hpp"
#include "genfile/SNPIdentifyingDataFilteringSNPDataSource.hpp"
#include "genfile/get_set.hpp"

namespace genfile {
	// Create a SNPIdentifyingDataFilteringSNPDataSource from the given source and the given sample indices.
	std::auto_ptr< SNPIdentifyingDataFilteringSNPDataSource > SNPIdentifyingDataFilteringSNPDataSource::create(
		SNPDataSource::UniquePtr source,
		SNPIdentifyingDataTest::UniquePtr test
	) {
		return SNPIdentifyingDataFilteringSNPDataSource::UniquePtr(
			new SNPIdentifyingDataFilteringSNPDataSource( source, test )
		) ;
	}

	SNPIdentifyingDataFilteringSNPDataSource::SNPIdentifyingDataFilteringSNPDataSource(
		SNPDataSource::UniquePtr source,
		SNPIdentifyingDataTest::UniquePtr test
	):
	 	m_source( source ),
		m_test( test )
	{
		m_source->reset_to_start() ;
	}

	SNPIdentifyingDataFilteringSNPDataSource::operator bool() const {
		return (*m_source) ;
	}

	unsigned int SNPIdentifyingDataFilteringSNPDataSource::number_of_samples() const {
		return m_source->number_of_samples() ;
	}

	SNPDataSource::OptionalSnpCount SNPIdentifyingDataFilteringSNPDataSource::total_number_of_snps() const {
		return OptionalSnpCount() ;
	}

	std::string SNPIdentifyingDataFilteringSNPDataSource::get_source_spec() const {
		return "snp-id-data-filtered:" + m_source->get_source_spec() ;
	}
	
	SNPDataSource const& SNPIdentifyingDataFilteringSNPDataSource::get_parent_source() const {
		return (*m_source) ;
	}

	SNPDataSource const& SNPIdentifyingDataFilteringSNPDataSource::get_base_source() const {
		return m_source->get_base_source() ;
	}

	void SNPIdentifyingDataFilteringSNPDataSource::reset_to_start_impl() {
		m_source->reset_to_start() ;
	}

	void SNPIdentifyingDataFilteringSNPDataSource::get_snp_identifying_data_impl( 
		IntegerSetter const& set_number_of_samples,
		StringSetter const& set_SNPID,
		StringSetter const& set_RSID,
		ChromosomeSetter const& set_chromosome,
		SNPPositionSetter const& set_SNP_position,
		AlleleSetter const& set_allele1,
		AlleleSetter const& set_allele2
	) {
		SNPIdentifyingData snp ;
		while(
			m_source->get_snp_identifying_data( snp )
		) {
			if( (*m_test)( snp )) {
				break ;
			}
			else {
				m_source->ignore_snp_probability_data() ;
			}
		}
		if( *this ) {
			set_number_of_samples( number_of_samples() ) ;
			set_SNPID( snp.get_SNPID() ) ;
			set_RSID( snp.get_rsid() ) ;
			set_chromosome( snp.get_position().chromosome() ) ;
			set_SNP_position( snp.get_position().position() ) ;
			set_allele1( snp.get_first_allele() ) ;
			set_allele2( snp.get_second_allele() ) ;
		}
	}

	VariantDataReader::UniquePtr SNPIdentifyingDataFilteringSNPDataSource::read_variant_data_impl() {
		return m_source->read_variant_data() ;
	}

	void SNPIdentifyingDataFilteringSNPDataSource::ignore_snp_probability_data_impl() {
		m_source->ignore_snp_probability_data() ;
	}
}
