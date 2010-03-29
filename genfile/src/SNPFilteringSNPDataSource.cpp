#include "genfile/Error.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPIdentifyingDataTest.hpp"
#include "genfile/SNPFilteringSNPDataSource.hpp"

namespace genfile {
	// Create a SNPFilteringSNPDataSource from the given source and the given sample indices.
	std::auto_ptr< SNPFilteringSNPDataSource > SNPFilteringSNPDataSource::create(
		SNPDataSource& source,
		std::auto_ptr< SNPIdentifyingDataTest > snp_inclusion_test
	) {
		return std::auto_ptr< SNPFilteringSNPDataSource >(
			new SNPFilteringSNPDataSource( source, snp_inclusion_test )
		) ;
	}

	SNPFilteringSNPDataSource::SNPFilteringSNPDataSource(
		SNPDataSource& source,
		std::auto_ptr< SNPIdentifyingDataTest > snp_inclusion_test
	):
	 	m_source( source ),
		m_snp_inclusion_test( snp_inclusion_test )
	{
		m_indices_of_excluded_snps = get_indices_of_excluded_snps() ;
		m_source.reset_to_start() ;
	}

	std::set< std::size_t > SNPFilteringSNPDataSource::get_indices_of_excluded_snps() {
		std::set< std::size_t > indices_of_excluded_snps ;
		std::string SNPID, RSID ;
		GenomePosition position ;
		char allele1, allele2 ;
		while( m_source.get_snp_identifying_data(
				genfile::ignore(),
				genfile::set_value( SNPID ),
				genfile::set_value( RSID ),
				genfile::set_value( position.chromosome() ),
				genfile::set_value( position.position() ),
				genfile::set_value( allele1 ),
				genfile::set_value( allele2 )
			)
		) {
			if( !m_snp_inclusion_test->operator()( SNPID, RSID, position, allele1, allele2 )) {
				indices_of_excluded_snps.insert( indices_of_excluded_snps.end(), m_source.number_of_snps_read() ) ;
			}
			m_source.ignore_snp_probability_data() ;
		}
		return indices_of_excluded_snps ;
	}

	SNPFilteringSNPDataSource::operator bool() const {
		return m_source ;
	}

	unsigned int SNPFilteringSNPDataSource::number_of_samples() const {
		return m_source.number_of_samples() ;
	}

	unsigned int SNPFilteringSNPDataSource::total_number_of_snps() const {
		return m_source.total_number_of_snps() - m_indices_of_excluded_snps.size() ;
	}

	unsigned int SNPFilteringSNPDataSource::total_number_of_snps_before_filtering() const {
		return m_source.total_number_of_snps() ;
	}

	SNPIdentifyingDataTest const& SNPFilteringSNPDataSource::get_snp_inclusion_test() const {
		return *m_snp_inclusion_test ;
	}

	void SNPFilteringSNPDataSource::reset_to_start_impl() {
		m_source.reset_to_start() ;
	}

	void SNPFilteringSNPDataSource::get_snp_identifying_data_impl( 
		IntegerSetter const& set_number_of_samples,
		StringSetter const& set_SNPID,
		StringSetter const& set_RSID,
		ChromosomeSetter const& set_chromosome,
		SNPPositionSetter const& set_SNP_position,
		AlleleSetter const& set_allele1,
		AlleleSetter const& set_allele2
	) {
		while( m_source && m_indices_of_excluded_snps.find( m_source.number_of_snps_read() ) != m_indices_of_excluded_snps.end() ) {
			m_source.get_snp_identifying_data(
				ignore(),
				ignore(),
				ignore(),
				ignore(),
				ignore(),
				ignore(),
				ignore()
			) ;
			m_source.ignore_snp_probability_data() ;
		}
		m_source.get_snp_identifying_data(
			set_number_of_samples,
			set_SNPID,
			set_RSID,
			set_chromosome,
			set_SNP_position,
			set_allele1,
			set_allele2
		) ;
	}

	void SNPFilteringSNPDataSource::read_snp_probability_data_impl(
		GenotypeProbabilitySetter const& set_genotype_probabilities
	) {
		m_source.read_snp_probability_data( set_genotype_probabilities ) ;
	}

	void SNPFilteringSNPDataSource::ignore_snp_probability_data_impl() {
		m_source.ignore_snp_probability_data() ;
	}
}
