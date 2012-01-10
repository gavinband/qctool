#include <vector>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/MergingSNPDataSource.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/Error.hpp"

namespace genfile {
	MergingSNPDataSource::MergingSNPDataSource()
	{}

	MergingSNPDataSource::~MergingSNPDataSource() {
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			delete m_sources[i] ;
		}
	}

	void MergingSNPDataSource::add_source(
		SNPDataSource::UniquePtr source
	) {
		if( m_sources.size() > 0 && source->number_of_samples() != m_sources[0]->number_of_samples() ) {
			throw BadArgumentError(
				"genfile::MergingSNPDataSource::add_source()",
				"Cohort "
				+ string_utils::to_string( m_sources.size() + 1 )
				+ " has the wrong number of samples ("
				+ string_utils::to_string( source->number_of_samples() )
				+ " instead of "
				+ string_utils::to_string( m_sources[0]->number_of_samples() )
				+ ")"
				) ;
		}
		m_sources.push_back( 0 ) ;
		m_sources.back() = source.release() ;
	
		// read first SNP from this source.
		get_top_snp_in_source( m_sources.size() - 1 ) ;
	}

	void MergingSNPDataSource::get_top_snp_in_source( std::size_t source_i ) {
		assert( source_i < m_sources.size() ) ;
		SNPIdentifyingData snp ;
		if( m_sources[ source_i ]->get_snp_identifying_data( snp ) ) {
			m_current_snps.insert(
				std::make_pair( snp, source_i )
			) ;
		}
	}

	MergingSNPDataSource::operator bool() const {
		return !m_current_snps.empty() ;
	}

	// Return the number of samples represented in the snps in this source.
	unsigned int MergingSNPDataSource::number_of_samples() const {
		if( m_sources.empty() ) {
			return 0 ;
		}
		else {
			return m_sources[0]->number_of_samples() ;
		}
	}

	// Return the total number of snps the source contains.
	unsigned int MergingSNPDataSource::total_number_of_snps() const {
		unsigned int result = 0 ;
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			result += m_sources[i]->total_number_of_snps() ;
		}
		return result ;
	}

	// Return a string identifying the source of the SNP data
	std::string MergingSNPDataSource::get_source_spec() const {
		std::string result = "merge:";
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			if( i > 0 ) {
				result += "," ;
			}
			result += m_sources[i]->get_source_spec() ;
		}
		return result ;
	}

	void MergingSNPDataSource::get_snp_identifying_data_impl( 
		IntegerSetter const& set_number_of_samples,
		StringSetter const& set_SNPID,
		StringSetter const& set_RSID,
		ChromosomeSetter const& set_chromosome,
		SNPPositionSetter const& set_SNP_position,
		AlleleSetter const& set_allele1,
		AlleleSetter const& set_allele2
	) {
		if( m_current_snps.size() > 0 ) {
			SNPIdentifyingData const& snp = m_current_snps.begin()->first ;
			set_number_of_samples( number_of_samples() ) ;
			set_SNPID( snp.get_SNPID() ) ;
			set_RSID( snp.get_rsid() ) ;
			set_chromosome( snp.get_position().chromosome() ) ;
			set_SNP_position( snp.get_position().position() ) ;
			set_allele1( snp.get_first_allele() ) ;
			set_allele2( snp.get_second_allele() ) ;
		}
	}

	VariantDataReader::UniquePtr MergingSNPDataSource::read_variant_data_impl() {
		VariantDataReader::UniquePtr result ;
		assert( m_current_snps.size() > 0 ) ;
		std::size_t source_i = m_current_snps.begin()->second ;
		result = m_sources[ source_i ]->read_variant_data() ;
		m_current_snps.erase( m_current_snps.begin() ) ;
		get_top_snp_in_source( source_i ) ;
		return result ;
	}

	void MergingSNPDataSource::ignore_snp_probability_data_impl() {
		assert( m_current_snps.size() > 0 ) ;
		std::size_t source_i = m_current_snps.begin()->second ;
		m_sources[ source_i ]->ignore_snp_probability_data() ;
		m_current_snps.erase( m_current_snps.begin() ) ;
		get_top_snp_in_source( source_i ) ;
	}

	void MergingSNPDataSource::reset_to_start_impl() {
		m_current_snps.clear() ;
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			m_sources[ i ]->reset_to_start() ;
			get_top_snp_in_source( i ) ;
		}
	}
}
