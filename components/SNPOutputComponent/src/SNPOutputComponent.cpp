
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>
#include <Eigen/Core>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/get_set_eigen.hpp"
#include "genfile/vcf/get_set_eigen.hpp"
#include "components/SNPOutputComponent/SNPOutputComponent.hpp"

namespace impl {
	struct QCDBSNPDataSourceIndex: public SNPDataSourceIndex {
	public:
		typedef std::auto_ptr< QCDBSNPDataSourceIndex > UniquePtr ;
		typedef boost::shared_ptr< QCDBSNPDataSourceIndex > SharedPtr ;
		typedef std::map< std::string, std::pair< std::vector< std::string >, std::string > > Metadata ;
		static UniquePtr create(
			std::string const& filename,
			genfile::CohortIndividualSource const& samples,
			std::string const& analysis_name,
			Metadata const& metadata
		) ;

		static SharedPtr create_shared(
			std::string const& filename,
			genfile::CohortIndividualSource const& samples,
			std::string const& analysis_name,
			Metadata const& metadata
		) ;
	public:
		QCDBSNPDataSourceIndex(
			std::string const& filename,
			genfile::CohortIndividualSource const& samples,
			std::string const& analysis_name,
			Metadata const& metadata
		) ;
		
		void add_index_entry( genfile::SNPIdentifyingData2 const& snp, genfile::SNPDataSink& sink ) ;
		
	private:
		std::string const m_table_name ;
		qcdb::DBOutputter m_outputter ;
		genfile::CohortIndividualSource const& m_samples ;
		db::Connection::StatementPtr m_insert_stmt ;
	private:
		void setup() ;
	} ;

}

void SNPOutputComponent::declare_options( appcontext::OptionProcessor& options ) {
}


SNPOutputComponent::SNPOutputComponent(
	genfile::CohortIndividualSource const& samples,
	appcontext::OptionProcessor const& options,
	appcontext::UIContext& ui_context
):
	m_samples( samples ),
	m_options( options ),
	m_ui_context( ui_context )
{}

void SNPOutputComponent::setup( genfile::SNPDataSink& sink, genfile::SNPDataSourceProcessor& processor ) {
	impl::SNPOutputter::UniquePtr outputter = impl::SNPOutputter::create( m_samples, sink ) ;

	if( m_options.check( "-write-index" )) {
		std::string const filename = m_options.get< std::string >( "-write-index" ) ;
		impl::QCDBSNPDataSourceIndex::SharedPtr index
			= impl::QCDBSNPDataSourceIndex::create_shared(
				filename,
				m_samples,
				m_options.get< std::string >( "-analysis-name" ),
				m_options.get_values_as_map()
		) ;
		outputter->send_index_to(
			boost::bind(
				&impl::SNPDataSourceIndex::add_index_entry,
				index,
				_1, _2
			)
		) ;
	}

	processor.add_callback( genfile::SNPDataSourceProcessor::Callback::UniquePtr( outputter.release() ) ) ;
}

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
	SNPOutputter::UniquePtr SNPOutputter::create( genfile::CohortIndividualSource const& samples, genfile::SNPDataSink& sink ) {
		return SNPOutputter::UniquePtr( new SNPOutputter( samples, sink ) ) ;
	}

	SNPOutputter::SNPOutputter( genfile::CohortIndividualSource const& samples, genfile::SNPDataSink& sink ):
		m_samples( samples ),
		m_manage( false ),
		m_sink( &sink )
	{
		assert( m_sink ) ;
	}
	
	void SNPOutputter::send_index_to( IndexCallback callback ) {
		m_index_callback = callback ;
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
		if( m_index_callback ) {
			m_index_callback( snp, *m_sink ) ;
		}
		{
			genfile::vcf::GenotypeSetter< Eigen::MatrixBase< Eigen::MatrixXd > > setter( m_genotypes ) ;
			data_reader.get( "genotypes", setter ) ;
		}
		send_results_to_sink( *m_sink, snp, m_genotypes ) ;
	}

	void SNPOutputter::end_processing_snps() {
		// nothing to do.
	}
	
	QCDBSNPDataSourceIndex::UniquePtr QCDBSNPDataSourceIndex::create(
		std::string const& filename,
		genfile::CohortIndividualSource const& samples,
		std::string const& analysis_name,
		Metadata const& metadata
	) {
		return QCDBSNPDataSourceIndex::UniquePtr(
			new QCDBSNPDataSourceIndex(
				filename,
				samples,
				analysis_name,
				metadata
			)
		) ;
	}

	QCDBSNPDataSourceIndex::SharedPtr QCDBSNPDataSourceIndex::create_shared(
		std::string const& filename,
		genfile::CohortIndividualSource const& samples,
		std::string const& analysis_name,
		Metadata const& metadata
	) {
		return SharedPtr( create( filename, samples, analysis_name, metadata ).release() ) ;
	}
	
	QCDBSNPDataSourceIndex::QCDBSNPDataSourceIndex(
		std::string const& filename,
		genfile::CohortIndividualSource const& samples,
		std::string const& analysis_name,
		Metadata const& metadata
	):
		m_outputter(
			filename,
			analysis_name,
			metadata
		),
		m_samples( samples ),
		m_table_name( "DatafileIndex" )
	{
		setup() ;
	}
	
	void QCDBSNPDataSourceIndex::setup() {
		db::Connection& connection = m_outputter.connection() ;
		db::Connection::ScopedTransactionPtr transaction = connection.open_transaction( 240 ) ;
		connection.run_statement(
			"CREATE TABLE IF NOT EXISTS "
			+ m_table_name + " "
			+ " ("
			+ "   analysis_id INTEGER NOT NULL REFERENCES Entity( id ), "
			+ "   variant_id INTEGER NOT NULL REFERENCES Variant( id ), "
			+ "   file_id INTEGER NOT NULL REFERENCES Entity( id ), "
			+ "   file_offset INTEGER "
			+ ")"
		) ;

		connection.run_statement(
			"CREATE VIEW IF NOT EXISTS "
			+ m_table_name + "View "
			+ "AS"
			" SELECT A.id AS analysis_id, A.name AS analysis, V.id AS variant_id, V.chromosome AS chromosome, V.position AS position, V.rsid AS rsid, V.alleleA AS alleleA, V.alleleB AS alleleB,"
			" File.name AS filename, I.file_offset AS file_offset"
			" FROM " + m_table_name + " I"
			" INNER JOIN Entity A ON A.id = I.analysis_id"
			" INNER JOIN Entity File ON File.id = I.file_id"
			" INNER JOIN Variant V ON V.id = I.variant_id"
		) ;
		
		m_outputter.get_or_create_entity_data(
			m_outputter.analysis_id(),
			m_outputter.get_or_create_entity( "table", "Table holding results of an analysis" ),
			m_table_name
 		) ;
		
		m_insert_stmt = m_outputter.connection().get_statement(
			"INSERT INTO " + m_table_name + " VALUES( ?1, ?2, ?3, ?4 )"
		) ;
	}
	
	void QCDBSNPDataSourceIndex::add_index_entry( genfile::SNPIdentifyingData2 const& snp, genfile::SNPDataSink& sink ) {
		genfile::SNPDataSink::SinkPos const pos = sink.get_stream_pos() ;
		db::Connection::RowId const file_id = m_outputter.get_or_create_entity( pos.first->get_spec(), "File holding per-variant data." ) ;
		db::Connection::RowId const variant_id = m_outputter.get_or_create_variant( snp ) ;
		m_insert_stmt
			->bind( 1, m_outputter.analysis_id() )
			.bind( 2, variant_id )
			.bind( 3, file_id )
			.bind( 4, int64_t( pos.second ) )
			.step() ;

		m_insert_stmt->reset() ;
	}
}
