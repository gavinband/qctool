#include <string>
#include "qcdb/DBOutputter.hpp"
#include "components/HaplotypeFrequencyComponent/DBOutputter.hpp"

namespace haplotype_frequency_component {
	DBOutputter::UniquePtr DBOutputter::create( std::string const& filename, std::string const& analysis_name, DBOutputter::Metadata const& metadata ) {
		return UniquePtr( new DBOutputter( filename, analysis_name, metadata ) ) ;
	}

	DBOutputter::SharedPtr DBOutputter::create_shared( std::string const& filename, std::string const& analysis_name, DBOutputter::Metadata const& metadata ) {
		return DBOutputter::SharedPtr( new DBOutputter( filename, analysis_name, metadata ) ) ;
	}

	DBOutputter::DBOutputter( std::string const& filename, std::string const& analysis_name, DBOutputter::Metadata const& metadata ):
		qcdb::DBOutputter( filename, analysis_name, metadata ),
		m_max_transaction_count( 10000 )
	{
		db::Connection::ScopedTransactionPtr transaction = connection().open_transaction() ;
		connection().run_statement(
			"CREATE TABLE IF NOT EXISTS PairwiseSummaryData ( "
			"analysis_id INT NOT NULL, "
			"variant1_id INT NOT NULL, "
			"variant2_id INT NOT NULL, "
			"variable_id INT NOT NULL, "
			"value NULL, "
			"FOREIGN KEY( analysis_id ) REFERENCES Entity( id ), "
			"FOREIGN KEY( variant1_id ) REFERENCES Variant( id ), "
			"FOREIGN KEY( variant2_id ) REFERENCES Variant( id ), "
			"FOREIGN KEY( variable_id ) REFERENCES Entity( id ), "
			"UNIQUE( variant1_id, variant2_id, analysis_id, variable_id ) "
			")"
		) ;
		// Typical use is to look for all SNPs in r^2 with a given one, so add an index for this.
		connection().run_statement(
			"CREATE INDEX PairwiseSummaryDataIndex ON PairwiseSummaryData( analysis_id, variant1_id, variable_id )"
		) ;
		connection().run_statement(
			"CREATE VIEW IF NOT EXISTS PairwiseSummaryDataView AS "
			"SELECT "
			"V1.id AS variant1_id, V1.chromosome AS variant1_chromosome, V1.position AS variant1_position, V1.rsid AS variant1_rsid, "
			"V2.id AS variant2_id, V2.chromosome AS variant2_chromosome, V2.position AS variant2_position, V2.rsid AS variant2_rsid, "
			"Analysis.id AS analysis_id, Analysis.name AS analysis, "
			"Variable.id AS variable_id, Variable.name AS variable, "
			"PSD.value AS value "
			"FROM PairwiseSummaryData PSD "
			"INNER JOIN Variant V1 ON V1.id == PSD.variant1_id "
			"INNER JOIN Variant V2 ON V2.id == PSD.variant2_id "
			"INNER JOIN Entity Analysis ON Analysis.id == PSD.analysis_id "
			"INNER JOIN Entity Variable ON Variable.id == PSD.variable_id "
		) ;
		
		m_variable_class_id = get_or_create_entity( "per-variant pair variable", "per-variant pair variable values" ) ; 
		
		construct_statements() ;
	}

	DBOutputter::~DBOutputter() {
		write_data( m_data ) ;
	}

	void DBOutputter::operator()(
		std::string const& cohort,
		genfile::SNPIdentifyingData const& source_snp,
		genfile::SNPIdentifyingData const& target_snp,
		std::string const& variable,
		genfile::VariantEntry const& value
	) {
		m_data.resize( m_data.size() + 1 ) ;
		m_data.back().get<0>() = cohort ;
		m_data.back().get<1>() = source_snp ;
		m_data.back().get<2>() = target_snp ;
		m_data.back().get<3>() = variable ;
		m_data.back().get<4>() = value ;

		if( m_data.size() == m_max_transaction_count ) {
			write_data( m_data ) ;
			m_data.clear() ;
		}
	}

	void DBOutputter::construct_statements() {
		m_insert_summarydata_statement = connection().get_statement(
			"INSERT OR REPLACE INTO PairwiseSummaryData (  analysis_id, variant1_id, variant2_id, variable_id, value ) "
			"VALUES( ?1, ?2, ?3, ?4, ?5 )"
		) ;
	}

	void DBOutputter::write_data( Data const& data ) {
		db::Connection::ScopedTransactionPtr transaction = connection().open_transaction( 240 ) ; // wait 4 minutes if we have to.

		if( !transaction.get() ) {
			throw genfile::OperationFailedError( "HaplotypeFrequencyComponent::DBOutputter::write_data()", connection().get_spec(), "Opening transaction." ) ;
		}
		
		std::vector< db::Connection::RowId > variant1_ids( data.size() ) ;
		std::vector< db::Connection::RowId > variant2_ids( data.size() ) ;
		for( std::size_t i = 0; i < m_data.size(); ++i ) {
			if( i == 0 || data[i].get<1>() != data[i-1].get<1>() ) {
				variant1_ids[i] = get_or_create_variant( data[i].get<1>() ) ;
			} else {
				variant1_ids[i] = variant1_ids[i-1] ;
			}
			if( i == 0 || data[i].get<2>() != data[i-1].get<2>() ) {
				variant2_ids[i] = get_or_create_variant( data[i].get<2>() ) ;
			} else {
				variant2_ids[i] = variant2_ids[i-1] ;
			}
		}

		for( std::size_t i = 0; i < m_data.size(); ++i ) {
			db::Connection::RowId const variable_id = get_or_create_entity( data[i].get<3>(), data[i].get<3>(), m_variable_class_id ) ;

			store_comparison(
				variant1_ids[i],
				variant2_ids[i],
				analysis_id(),
				variable_id,
				data[i].get<4>()
			) ;
		}
	}

	void DBOutputter::store_comparison(
		db::Connection::RowId const variant1_id,
		db::Connection::RowId const variant2_id,
		db::Connection::RowId const analysis_id,
		db::Connection::RowId const variable_id,
		genfile::VariantEntry const& value
	) {
		m_insert_summarydata_statement
			->bind( 1, analysis_id )
			.bind( 2, variant1_id )
			.bind( 3, variant2_id )
			.bind( 4, variable_id )
			.bind( 5, value )
			.step() ;
		m_insert_summarydata_statement->reset() ;
	}
}
