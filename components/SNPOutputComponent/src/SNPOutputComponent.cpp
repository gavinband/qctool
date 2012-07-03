
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <boost/bind.hpp>
#include <Eigen/Core>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/get_set_eigen.hpp"
#include "genfile/vcf/get_set_eigen.hpp"
#include "components/SNPOutputComponent/SNPOutputComponent.hpp"

namespace {
	genfile::VariantEntry get_sample_entry( genfile::CohortIndividualSource const& samples, std::string const& name, std::size_t i ) {
		return samples.get_entry( i, name ) ;
	}

	typedef std::map< std::string, std::vector< genfile::VariantEntry > > Info ;
	void send_results_to_sink(
		genfile::SNPDataSink& sink,
		genfile::SNPIdentifyingData const& snp,
		Eigen::MatrixXd const& genotypes,
		Info const& info = Info()
	) {
		sink.write_snp(
			genotypes.rows(),
			snp,
			genfile::GenotypeGetter< Eigen::MatrixXd >( genotypes, 0ul ),
			genfile::GenotypeGetter< Eigen::MatrixXd >( genotypes, 1ul ),
			genfile::GenotypeGetter< Eigen::MatrixXd >( genotypes, 2ul ),
			info
		) ;
	}
}

namespace impl {
	SNPOutputter::UniquePtr SNPOutputter::create( genfile::CohortIndividualSource const& samples, genfile::SNPDataSink::UniquePtr sink ) {
		return SNPOutputter::UniquePtr( new SNPOutputter( samples, sink ) ) ;
	}

	SNPOutputter::UniquePtr SNPOutputter::create( genfile::CohortIndividualSource const& samples, genfile::SNPDataSink& sink ) {
		return SNPOutputter::UniquePtr( new SNPOutputter( samples, sink ) ) ;
	}

	SNPOutputter::SNPOutputter( genfile::CohortIndividualSource const& samples, genfile::SNPDataSink::UniquePtr sink ):
		m_samples( samples ),
		m_manage( true ),
		m_sink( sink.release() )
	{
		assert( m_sink ) ;
	}

	SNPOutputter::SNPOutputter( genfile::CohortIndividualSource const& samples, genfile::SNPDataSink& sink ):
		m_samples( samples ),
		m_manage( false ),
		m_sink( &sink )
	{
		assert( m_sink ) ;
	}
	
	SNPOutputter::~SNPOutputter() {
		if( m_manage ) {
			delete m_sink ;
		}
	}

	void SNPOutputter::begin_processing_snps( std::size_t number_of_samples ) {
		m_sink->set_sample_names( number_of_samples, boost::bind( get_sample_entry, boost::ref( m_samples ), "ID_1", _1 ) ) ;
	}

	void SNPOutputter::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
		{
			genfile::vcf::GenotypeSetter< Eigen::MatrixBase< Eigen::MatrixXd > > setter( m_genotypes ) ;
			data_reader.get( "genotypes", setter ) ;
		}
		send_results_to_sink( *m_sink, snp, m_genotypes ) ;
	}

	void SNPOutputter::end_processing_snps() {
		// nothing to do.
	}
}

void SNPOutputComponent::setup( genfile::CohortIndividualSource const& samples, genfile::SNPDataSink::UniquePtr sink, genfile::SNPDataSourceProcessor& processor ) {
	processor.add_callback( impl::SNPOutputter::create( samples, sink )) ;
}

void SNPOutputComponent::setup( genfile::CohortIndividualSource const& samples, genfile::SNPDataSink& sink, genfile::SNPDataSourceProcessor& processor ) {
	processor.add_callback( impl::SNPOutputter::create( samples, sink )) ;
}

