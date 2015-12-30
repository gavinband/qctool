
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "genfile/Error.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/VariantIdentifyingDataTest.hpp"
#include "genfile/SNPIdentifyingDataFilteringSNPDataSource.hpp"
#include "genfile/get_set.hpp"

namespace genfile {
	// Create a VariantIdentifyingDataFilteringSNPDataSource from the given source and the given sample indices.
	std::auto_ptr< VariantIdentifyingDataFilteringSNPDataSource > VariantIdentifyingDataFilteringSNPDataSource::create(
		SNPDataSource::UniquePtr source,
		VariantIdentifyingDataTest::UniquePtr test
	) {
		return VariantIdentifyingDataFilteringSNPDataSource::UniquePtr(
			new VariantIdentifyingDataFilteringSNPDataSource( source, test )
		) ;
	}

	std::auto_ptr< VariantIdentifyingDataFilteringSNPDataSource > VariantIdentifyingDataFilteringSNPDataSource::create(
		SNPDataSource::UniquePtr source,
		VariantIdentifyingDataTest const& test
	) {
		return VariantIdentifyingDataFilteringSNPDataSource::UniquePtr(
			new VariantIdentifyingDataFilteringSNPDataSource( source, test )
		) ;
	}

	VariantIdentifyingDataFilteringSNPDataSource::VariantIdentifyingDataFilteringSNPDataSource(
		SNPDataSource::UniquePtr source,
		VariantIdentifyingDataTest::UniquePtr test
	):
	 	m_source( source ),
		m_manage_test( true ),
		m_test( test.release() )
	{
		m_source->reset_to_start() ;
	}

	VariantIdentifyingDataFilteringSNPDataSource::VariantIdentifyingDataFilteringSNPDataSource(
		SNPDataSource::UniquePtr source,
		VariantIdentifyingDataTest const& test
	):
	 	m_source( source ),
		m_manage_test( false ),
		m_test( &test )
	{
		m_source->reset_to_start() ;
	}
	
	VariantIdentifyingDataFilteringSNPDataSource::~VariantIdentifyingDataFilteringSNPDataSource() {
		if( m_manage_test ) {
			delete m_test ;
		}
	}
	

	VariantIdentifyingDataFilteringSNPDataSource::operator bool() const {
		return (*m_source) ;
	}

	SNPDataSource::Metadata VariantIdentifyingDataFilteringSNPDataSource::get_metadata() const {
		return m_source->get_metadata() ;
	}

	unsigned int VariantIdentifyingDataFilteringSNPDataSource::number_of_samples() const {
		return m_source->number_of_samples() ;
	}

	void VariantIdentifyingDataFilteringSNPDataSource::get_sample_ids( GetSampleIds getter ) const {
		return m_source->get_sample_ids( getter ) ;
	}

	SNPDataSource::OptionalSnpCount VariantIdentifyingDataFilteringSNPDataSource::total_number_of_snps() const {
		return OptionalSnpCount() ;
	}

	std::string VariantIdentifyingDataFilteringSNPDataSource::get_source_spec() const {
		return "snp-id-data-filtered:" + m_source->get_source_spec() ;
	}
	
	SNPDataSource const& VariantIdentifyingDataFilteringSNPDataSource::get_parent_source() const {
		return (*m_source) ;
	}

	SNPDataSource const& VariantIdentifyingDataFilteringSNPDataSource::get_base_source() const {
		return m_source->get_base_source() ;
	}

	void VariantIdentifyingDataFilteringSNPDataSource::reset_to_start_impl() {
		m_source->reset_to_start() ;
	}

	void VariantIdentifyingDataFilteringSNPDataSource::get_snp_identifying_data_impl( VariantIdentifyingData* result ) {
		VariantIdentifyingData snp ;
		while(
			m_source->get_snp_identifying_data( &snp )
		) {
			if( (*m_test)( snp )) {
				break ;
			}
			else {
				m_filtered_out_snp_signal( snp ) ;
				m_source->ignore_snp_probability_data() ;
			}
		}
		if( *this ) {
			*result = snp ;
		}
	}

	VariantDataReader::UniquePtr VariantIdentifyingDataFilteringSNPDataSource::read_variant_data_impl() {
		return m_source->read_variant_data() ;
	}

	void VariantIdentifyingDataFilteringSNPDataSource::ignore_snp_probability_data_impl() {
		m_source->ignore_snp_probability_data() ;
	}
	
	void VariantIdentifyingDataFilteringSNPDataSource::send_filtered_out_SNPs_to( SNPSignal::slot_type callback ) {
		m_filtered_out_snp_signal.connect( callback ) ;
	}
}
