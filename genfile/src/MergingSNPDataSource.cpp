
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/MergingSNPDataSource.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/Error.hpp"

namespace genfile {
	
	std::vector< std::string > MergingSNPDataSource::get_merge_strategies() {
		std::vector< std::string > result ;
		result.push_back( "keep-all" ) ;
		result.push_back( "drop-duplicates" ) ;
		return result ;
	}
	
	MergingSNPDataSource::UniquePtr MergingSNPDataSource::create(
		std::string const& merge_strategy,
		SNPIdentifyingData::CompareFields const& compare_fields
	) {
		MergingSNPDataSource::UniquePtr result ;
		if( merge_strategy == "keep-all" ) {
			result.reset( new KeepAllStrategyMergingSNPDataSource( compare_fields ) ) ;
		}
		else if( merge_strategy == "drop-duplicates" ) {
			result.reset( new DropDuplicatesStrategyMergingSNPDataSource( compare_fields ) ) ;
		}
		else {
			throw genfile::BadArgumentError( "genfile::MergingSNPDataSource::create()", "merge_strategy=\"" + merge_strategy + "\"" ) ;
		}
		return result ;
	}
	
	MergingSNPDataSource::MergingSNPDataSource( SNPIdentifyingData::CompareFields const& compare_fields ):
		m_current_snps( compare_fields )
	{}

	MergingSNPDataSource::~MergingSNPDataSource() {
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			delete m_sources[i] ;
		}
	}

	void MergingSNPDataSource::add_source(
		SNPDataSource::UniquePtr source,
		std::string const& id_prefix
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
		m_merge_id_prefixes.push_back( id_prefix ) ;

		// Get metadata.
		// Currently we just take the metadata from the first source.
		if( m_sources.size() == 1 ) {
			m_metadata = m_sources[0]->get_metadata() ;
		}

		// read first SNP from this source.
		get_top_snp_in_source( m_sources.size() - 1 ) ;
	}

	SNPDataSource::Metadata MergingSNPDataSource::get_metadata() const {
		return m_metadata ;
	}

	void MergingSNPDataSource::get_top_snp_in_source( std::size_t source_i ) {
		assert( source_i < m_sources.size() ) ;
		VariantIdentifyingData snp ;
		if( m_sources[ source_i ]->get_snp_identifying_data( &snp ) ) {
			m_current_snps.insert(
				std::make_pair( snp, source_i )
			) ;
		}
	}

	void MergingSNPDataSource::discard_top_snp_and_get_next_candidate() {
		assert( m_current_snps.size() > 0 ) ;
		std::size_t source_i = m_current_snps.begin()->second ;
		m_current_snps.erase( m_current_snps.begin() ) ;
		get_top_snp_in_source( source_i ) ;
	}

	KeepAllStrategyMergingSNPDataSource::KeepAllStrategyMergingSNPDataSource( SNPIdentifyingData::CompareFields const& compare_fields ):
		MergingSNPDataSource( compare_fields )
	{}

	void KeepAllStrategyMergingSNPDataSource::discard_top_snp_and_get_next() {
		discard_top_snp_and_get_next_candidate() ;
	}

	void DropDuplicatesStrategyMergingSNPDataSource::discard_top_snp_and_get_next() {
		assert( current_snps().size() > 0 ) ;
		
		genfile::SNPIdentifyingData const current_snp = current_snps().begin()->first ;
		discard_top_snp_and_get_next_candidate() ;

		while(
			current_snps().size() > 0
			&& current_snps().begin()->first.get_position() == current_snp.get_position()
		) {
			get_source( current_snps().begin()->second ).ignore_snp_probability_data() ;
			discard_top_snp_and_get_next_candidate() ;
		}
	}

	DropDuplicatesStrategyMergingSNPDataSource::DropDuplicatesStrategyMergingSNPDataSource( SNPIdentifyingData::CompareFields const& compare_fields ):
		MergingSNPDataSource( compare_fields )
	{}
		
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
	SNPDataSource::OptionalSnpCount MergingSNPDataSource::total_number_of_snps() const {
		unsigned int result = 0 ;
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			if( !m_sources[i]->total_number_of_snps() ) {
				return OptionalSnpCount() ;
			}
			result += *m_sources[i]->total_number_of_snps() ;
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
		VariantIdentifyingData* result
	) {
		if( m_current_snps.size() > 0 ) {
			SNPIdentifyingData const& snp = m_current_snps.begin()->first ;
			std::size_t source_index = m_current_snps.begin()->second ;
			*result = VariantIdentifyingData(
				m_merge_id_prefixes[ source_index ] + snp.get_SNPID(),
				snp.get_rsid(),
				snp.get_position(),
				snp.get_first_allele(),
				snp.get_second_allele()
			) ;
		}
	}

	VariantDataReader::UniquePtr MergingSNPDataSource::read_variant_data_impl() {
		VariantDataReader::UniquePtr result ;
		assert( m_current_snps.size() > 0 ) ;
		std::size_t source_i = m_current_snps.begin()->second ;
		result = m_sources[ source_i ]->read_variant_data() ;
		discard_top_snp_and_get_next() ;
		return result ;
	}

	void MergingSNPDataSource::ignore_snp_probability_data_impl() {
		assert( m_current_snps.size() > 0 ) ;
		std::size_t source_i = m_current_snps.begin()->second ;
		m_sources[ source_i ]->ignore_snp_probability_data() ;
		discard_top_snp_and_get_next() ;
	}

	void MergingSNPDataSource::reset_to_start_impl() {
		m_current_snps.clear() ;
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			m_sources[ i ]->reset_to_start() ;
			get_top_snp_in_source( i ) ;
		}
	}
}
