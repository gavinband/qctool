
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "genfile/Error.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPIdentifyingDataTest.hpp"
#include "genfile/SNPFilteringSNPDataSource.hpp"
#include "genfile/get_set.hpp"

namespace genfile {
	// Create a SNPFilteringSNPDataSource from the given source and the given sample indices.
	std::auto_ptr< SNPFilteringSNPDataSource > SNPFilteringSNPDataSource::create(
		SNPDataSource::UniquePtr source,
		IndexList const& indices_of_snps_to_include
	) {
		return SNPFilteringSNPDataSource::UniquePtr(
			new SNPFilteringSNPDataSource( source, indices_of_snps_to_include )
		) ;
	}

	SNPFilteringSNPDataSource::SNPFilteringSNPDataSource(
		SNPDataSource::UniquePtr source,
		IndexList indices_of_snps_to_include
	):
	 	m_source( source )
	{
		std::sort( indices_of_snps_to_include.begin(), indices_of_snps_to_include.end() ) ;
		for( std::size_t i = 0; i < m_source->total_number_of_snps(); ++i ) {
			if( !std::binary_search( indices_of_snps_to_include.begin(), indices_of_snps_to_include.end(), i )) {
				m_indices_of_excluded_snps.insert( i ) ;
			}
		}
		m_source->reset_to_start() ;
	}

	SNPFilteringSNPDataSource::operator bool() const {
		return (*m_source) ;
	}

	SNPDataSource::Metadata SNPFilteringSNPDataSource::get_metadata() const {
		return m_source->get_metadata() ;
	}

	unsigned int SNPFilteringSNPDataSource::number_of_samples() const {
		return m_source->number_of_samples() ;
	}

	void SNPFilteringSNPDataSource::get_sample_ids( GetSampleIds getter ) const {
		return m_source->get_sample_ids( getter ) ;
	}

	SNPDataSource::OptionalSnpCount SNPFilteringSNPDataSource::total_number_of_snps() const {
		if( m_source->total_number_of_snps() ) {
			return *m_source->total_number_of_snps() - m_indices_of_excluded_snps.size() ;
		}
		else {
			return OptionalSnpCount() ;
		}
	}

	SNPDataSource::OptionalSnpCount SNPFilteringSNPDataSource::total_number_of_snps_before_filtering() const {
		if( m_source->total_number_of_snps() ) {
			return *m_source->total_number_of_snps() ;
		}
		else {
			return OptionalSnpCount() ;
		}
	}
	
	std::string SNPFilteringSNPDataSource::get_source_spec() const {
		return "snp-filtered:" + m_source->get_source_spec() ;
	}
	
	SNPDataSource const& SNPFilteringSNPDataSource::get_parent_source() const {
		return (*m_source) ;
	}

	SNPDataSource const& SNPFilteringSNPDataSource::get_base_source() const {
		return m_source->get_base_source() ;
	}

	void SNPFilteringSNPDataSource::reset_to_start_impl() {
		m_source->reset_to_start() ;
	}

	void SNPFilteringSNPDataSource::get_snp_identifying_data_impl( 
		VariantIdentifyingData* result
	) {
		VariantIdentifyingData variant ;
		while( (*m_source) && m_indices_of_excluded_snps.find( m_source->number_of_snps_read() ) != m_indices_of_excluded_snps.end() ) {
			m_source->get_snp_identifying_data( &variant ) ;
			m_source->ignore_snp_probability_data() ;
		}
		m_source->get_snp_identifying_data( result ) ;
	}

	VariantDataReader::UniquePtr SNPFilteringSNPDataSource::read_variant_data_impl() {
		return m_source->read_variant_data() ;
	}

	void SNPFilteringSNPDataSource::ignore_snp_probability_data_impl() {
		m_source->ignore_snp_probability_data() ;
	}
}
