
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <fstream>
#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <Eigen/Core>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/get_set_eigen.hpp"
#include "genfile/vcf/get_set_eigen.hpp"
#include "components/SNPOutputComponent/SNPOutputComponent.hpp"

void SNPOutputComponent::declare_options( appcontext::OptionProcessor& options ) {
}


SNPOutputComponent::SNPOutputComponent(
	genfile::CohortIndividualSource const& samples,
	appcontext::OptionProcessor const& options,
	appcontext::UIContext& ui_context
):
	m_samples( samples ),
	m_options( options )
{}

void SNPOutputComponent::setup(
	genfile::SNPDataSink& sink,
	genfile::SNPDataSourceProcessor& processor
) {
	impl::SNPOutputter::UniquePtr outputter = impl::SNPOutputter::create( m_samples, sink ) ;
	processor.add_callback( genfile::SNPDataSourceProcessor::Callback::UniquePtr( outputter.release() ) ) ;
}

namespace {
	genfile::VariantEntry get_sample_entry( genfile::CohortIndividualSource const& samples, std::string const& name, std::size_t i ) {
		return samples.get_entry( i, name ) ;
	}
}

namespace impl {
	SNPOutputter::UniquePtr SNPOutputter::create(
		genfile::CohortIndividualSource const& samples,
		genfile::SNPDataSink& sink
	) {
		return SNPOutputter::UniquePtr( new SNPOutputter( samples, sink ) ) ;
	}

	SNPOutputter::SNPOutputter(
		genfile::CohortIndividualSource const& samples,
		genfile::SNPDataSink& sink
	):
		m_samples( samples ),
		m_manage( false ),
		m_sink( &sink )
	{
		assert( m_sink ) ;
	}
	
	void SNPOutputter::send_index_to( impl::SNPDataSourceIndex::UniquePtr index ) {
		m_index = index ;
	}
	
	SNPOutputter::~SNPOutputter() {
		if( m_manage ) {
			delete m_sink ;
		}
	}

	void SNPOutputter::begin_processing_snps( std::size_t number_of_samples, genfile::SNPDataSource::Metadata const& metadata ) {
		m_sink->set_metadata( metadata ) ;
		m_sink->set_sample_names( number_of_samples, boost::bind( get_sample_entry, boost::ref( m_samples ), "ID_1", _1 ) ) ;
	}
	
	void SNPOutputter::processed_snp( genfile::VariantIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
		if( m_index.get() ) {
			m_index->add_index_entry( snp, *m_sink ) ;
		}
		typedef std::map< std::string, std::vector< genfile::VariantEntry > > Info ;
		m_sink->write_variant_data( snp, data_reader, Info() ) ;
	}

	void SNPOutputter::end_processing_snps() {
		m_sink->finalise() ;
		if( m_index.get() ) {
			m_index->finalise() ;
		}
	}
}

