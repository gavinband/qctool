#include <string>
#include <memory>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/thread/thread.hpp>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/Error.hpp"
#include "db/Connection.hpp"
#include "db/SQLStatement.hpp"
#include "components/SNPSummaryComponent/DBOutputter.hpp"
#include "../../qctool_version_autogenerated.hpp"

namespace impl {
	DBOutputter::UniquePtr DBOutputter::create( std::string const& filename, std::string const& cohort_name, std::string const& data_source_spec, std::string const& exclusions_name ) {
		return UniquePtr( new DBOutputter( filename, cohort_name, data_source_spec, exclusions_name ) ) ;
	}
	DBOutputter::SharedPtr DBOutputter::create_shared( std::string const& filename, std::string const& cohort_name, std::string const& data_source_spec, std::string const& exclusions_name ) {
		return SharedPtr( new DBOutputter( filename, cohort_name, data_source_spec, exclusions_name ) ) ;
	}

	DBOutputter::DBOutputter( std::string const& filename, std::string const& cohort_name, std::string const& data_source_spec, std::string const& exclusions_name ):
		m_connection( db::Connection::create( filename )),
		m_max_transaction_count( 10000 ),
		m_cohort_name( cohort_name ),
		m_source_spec( data_source_spec ),
		m_exclusions_name( exclusions_name )
	{
		db::Connection::ScopedTransactionPtr transaction = m_connection->open_transaction() ;
		m_connection->run_statement(
			"CREATE TABLE IF NOT EXISTS Variant ( id INTEGER PRIMARY KEY, snpid TEXT, rsid TEXT, chromosome TEXT, position INTEGER, alleleA TEXT, alleleB TEXT )"
		) ;
		m_connection->run_statement(
			"CREATE INDEX IF NOT EXISTS Variant_index ON Variant( chromosome, position, rsid )"
		) ;
		m_connection->run_statement(
			"CREATE TABLE IF NOT EXISTS Entity ( id INTEGER PRIMARY KEY, name TEXT, description TEXT )"
		) ;
		m_connection->run_statement(
			"CREATE TABLE IF NOT EXISTS EntityData ( "
			"entity_id INTEGER NOT NULL, "
			"variable_id INTEGER NOT NULL, "
			"value TEXT, "
			"FOREIGN KEY (entity_id) REFERENCES Entity( id ), "
			"FOREIGN KEY (variable_id) REFERENCES Entity( id ) "
			")"
		) ;
		m_connection->run_statement(
			"CREATE TABLE IF NOT EXISTS SummaryData ( "
			"variant_id INT, analysis_id INT, variable_id INT, value FLOAT, "
			"FOREIGN KEY( variant_id ) REFERENCES Variant( id ), "
			"FOREIGN KEY( analysis_id ) REFERENCES Entity( id ), "
			"FOREIGN KEY( variable_id ) REFERENCES Entity( id ), "
			"UNIQUE( analysis_id, variant_id, variable_id ) "
			")"
		) ;
		m_connection->run_statement(
			"CREATE INDEX IF NOT EXISTS SummaryDataIndex ON SummaryData( analysis_id, variant_id, variable_id )"
		) ;
		m_connection->run_statement(
			"CREATE INDEX IF NOT EXISTS EntityDataIndex ON EntityData( entity_id, variable_id )"
		) ;
		m_connection->run_statement(
			"CREATE VIEW IF NOT EXISTS SummaryDataView AS "
			"SELECT          V.id AS variant_id, V.chromosome, V.position, V.rsid, "
			"SD.analysis_id, Analysis.name, Variable.id AS variable_id, Variable.name AS variable, "
			"SD.value AS value "
			"FROM SummaryData SD "
			"INNER JOIN Variant V ON V.id == SD.variant_id "
			"INNER JOIN Entity Analysis ON Analysis.id = SD.analysis_id "
			"INNER JOIN Entity Variable ON Variable.id = SD.variable_id"
		) ;
		m_connection->run_statement(
			"CREATE VIEW IF NOT EXISTS SNPFilterView AS "
			"SELECT          Analysis.name, V.chromosome, V.position, V.rsid, "
			"MAF.value AS 'MAF', "
			"HWE.value AS 'minus_log10_exact_HW_p_value', "
			"Missing.value AS 'missing_proportion', "
			"Info.value AS 'impute_info' "
			"FROM            Variant V "
			"INNER JOIN      Entity Analysis "
			"LEFT OUTER JOIN      SummaryData Missing "
			"ON          Missing.variable_id IN ( SELECT id FROM Entity WHERE name = 'missing proportion' ) "
			"AND         Missing.analysis_id = Analysis.id "
			"AND         Missing.variant_id == V.id "
			"LEFT OUTER JOIN      SummaryData MAF "
			"ON          MAF.variable_id = ( SELECT id FROM Entity WHERE name = 'minor_allele_frequency' ) "
			"AND         MAF.analysis_id = Analysis.id "
			"AND         MAF.variant_id == V.id "
			"LEFT OUTER JOIN      SummaryData HWE "
			"ON          HWE.variant_id == V.id "
			"AND         HWE.analysis_id = Analysis.id "
			"AND         HWE.variable_id == ( SELECT id FROM Entity WHERE name == 'minus_log10_exact_HW_p_value' ) "
			"LEFT OUTER JOIN      SummaryData Info "
			"ON          Info.variant_id == V.id "
			"AND         Info.analysis_id = Analysis.id "
			"AND         Info.variable_id == ( SELECT id FROM Entity WHERE name == 'impute_info' ) "
			"WHERE       EXISTS( SELECT * FROM SummaryData WHERE analysis_id == Analysis.id ) "
		) ;
		
		m_connection->run_statement(
			"CREATE VIEW IF NOT EXISTS SNPTestView AS "
			"SELECT        Analysis.name, V.chromosome, V.position, V.rsid, "
			"Beta.value AS beta, "
			"Se.value AS se, "
			"Pvalue.value AS pvalue "
			"FROM          Variant V "
			"INNER JOIN    Entity Analysis "
			"LEFT OUTER JOIN SummaryData Beta "
			"ON       Beta.variable_id IN ( SELECT id FROM Entity WHERE name == 'beta_1' ) "
			"AND      Beta.analysis_id = Analysis.id "
			"AND      Beta.variant_id = V.id "
			"LEFT OUTER JOIN SummaryData SE "
			"ON       SE.variable_id IN ( SELECT id FROM Entity WHERE name == 'se_1' ) "
			"AND      SE.analysis_id = Analysis.id "
			"AND      SE.variant_id = V.id "
			"LEFT OUTER JOIN SummaryData Pvalue "
			"ON       Pvalue.variable_id IN ( SELECT id FROM Entity WHERE name == 'p_value' ) "
			"AND      Pvalue.analysis_id = Analysis.id "
			"AND      Pvalue.variant_id = V.id "
			"WHERE       EXISTS( SELECT * FROM SummaryData WHERE analysis_id == Analysis.id ) "
		) ;
		construct_statements() ;
		
		m_analysis_id = get_or_create_entity( m_cohort_name ) ;
		get_or_create_entity_data( m_analysis_id, get_or_create_entity( "data_source" ), m_source_spec ) ;
		if( m_exclusions_name != "" ) {
			get_or_create_entity_data( m_analysis_id, get_or_create_entity( "sample_exclusions" ), m_exclusions_name ) ;
		}
		
		reset_statements() ;
	}

	DBOutputter::~DBOutputter() {
		write_data( m_data ) ;
	}

	void DBOutputter::operator()(
		std::size_t index,
		genfile::SNPIdentifyingData const& snp,
		std::string const& computation_name,
		std::string const& variable,
		genfile::VariantEntry const& value
	) {
		m_data.resize( m_data.size() + 1 ) ;
		m_data.back().get<0>() = snp ;
		m_data.back().get<1>() = variable ;
		m_data.back().get<2>() = value ;

		if( m_data.size() == m_max_transaction_count ) {
			write_data( m_data ) ;
			m_data.clear() ;
		}
	}

	void DBOutputter::construct_statements() {
		m_find_variant_statement = m_connection->get_statement(
			"SELECT id FROM Variant WHERE rsid == ?1 AND chromosome == ?2 AND position == ?3"
		) ;
		m_insert_variant_statement = m_connection->get_statement(
			"INSERT INTO Variant ( snpid, rsid, chromosome, position, alleleA, alleleB) "
			"VALUES( ?1, ?2, ?3, ?4, ?5, ?6 )"
		) ;
		m_find_entity_statement = m_connection->get_statement( "SELECT * FROM Entity E WHERE name == ?1" ) ;
		m_insert_entity_statement = m_connection->get_statement( "INSERT INTO Entity ( name, description ) VALUES ( ?1, ?2 )" ) ;
		m_find_entity_data_statement = m_connection->get_statement( "SELECT * FROM EntityData WHERE entity_id == ?1 AND variable_id == ?2" ) ;
		m_insert_entity_data_statement = m_connection->get_statement( "INSERT OR REPLACE INTO EntityData ( entity_id, variable_id, value ) VALUES ( ?1, ?2, ?3 )" ) ;
		m_insert_summarydata_statement = m_connection->get_statement(
			"INSERT OR REPLACE INTO SummaryData ( variant_id, analysis_id, variable_id, value ) "
			"VALUES( ?1, ?2, ?3, ?4 )"
		) ;
	}

	void DBOutputter::reset_statements() {
		m_find_variant_statement->reset() ;
		m_insert_variant_statement->reset() ;
		m_find_entity_statement->reset() ;
		m_insert_entity_statement->reset() ;
		m_find_entity_data_statement->reset() ;
		m_insert_entity_data_statement->reset() ;
		m_insert_summarydata_statement->reset() ;
	}
	
	void DBOutputter::write_data( Data const& data ) {
		db::Connection::ScopedTransactionPtr transaction ;

		for( std::size_t i = 0; i < 1000; ++i ) {
			try {
				transaction = m_connection->open_transaction() ;
				break ;
			}
			catch( db::StatementStepError const& e ) {
				// wait a tenth of a second
				std::cerr << "SNPSummaryComponent::DBOutputter::write_data(): failed to open transaction, trying again in 0.1s...\n" ;
				boost::this_thread::sleep( boost::posix_time::milliseconds( 10 ) ) ;
			}
			catch( ... ) {
				std::cerr << "SNPSummaryComponent::write_data(): OMG, a strange exception was caught.\n" ;
				boost::this_thread::sleep( boost::posix_time::milliseconds( 10 ) ) ;
			}
		}
		if( !transaction.get() ) {
			throw genfile::OperationFailedError( "SNPSummaryComponent::write_data()", m_connection->get_spec(), "Opening transaction." ) ;
		}
		for( std::size_t i = 0; i < data.size(); ++i ) {
			store_data(
				data[i].get<0>(),
				data[i].get<1>(),
				data[i].get<2>()
			) ;
		}

		reset_statements() ;
	}

	db::Connection::RowId DBOutputter::get_or_create_snp( genfile::SNPIdentifyingData const& snp ) const {
		m_find_variant_statement->reset()
			.bind( 1, snp.get_rsid() )
			.bind( 2, std::string( snp.get_position().chromosome() ))
			.bind( 3, snp.get_position().position() )
			.step()
		;
		if( m_find_variant_statement->empty() ) {
			m_insert_variant_statement
				->reset()
				.bind( 1, snp.get_SNPID() )
				.bind( 2, snp.get_rsid() )
				.bind( 3, std::string( snp.get_position().chromosome() ) )
				.bind( 4, snp.get_position().position() )
				.bind( 5, snp.get_first_allele())
				.bind( 6, snp.get_second_allele())
				.step()
			;

			return m_connection->get_last_insert_row_id() ;
		} else {
			return m_find_variant_statement->get< db::Connection::RowId >( 0 ) ;
		}
	}

	db::Connection::RowId DBOutputter::get_or_create_variable( std::string const& name ) const {
		db::Connection::RowId result ;

		m_find_entity_statement
			->reset()
			.bind( 1, name )
			.step() ;

		if( m_find_entity_statement->empty() ) {
			m_insert_entity_statement
				->reset()
				.bind( 1, name )
				.step() ;
				
			result = m_connection->get_last_insert_row_id() ;
			
			m_insert_entity_data_statement
				->reset()
				.bind( 1, result )
				.bind( 2, get_or_create_entity( "tool" ))
				.bind( 3, "qctool revision " + std::string( globals::qctool_revision ) )
				.step()
			;
		} else {
			result = m_find_entity_statement->get< db::Connection::RowId >( 0 ) ;
		}
		return result ;
	}

	db::Connection::RowId DBOutputter::get_or_create_entity( std::string const& name ) const {
		db::Connection::RowId result ;

		m_find_entity_statement
			->reset()
			.bind( 1, name ).step() ;

		if( m_find_entity_statement->empty() ) {
			m_insert_entity_statement
				->reset()
				.bind( 1, name )
				.step() ;
				
			result = m_connection->get_last_insert_row_id() ;
		} else {
			result = m_find_entity_statement->get< db::Connection::RowId >( 0 ) ;
		}
		return result ;
	}

	db::Connection::RowId DBOutputter::get_or_create_entity( std::string const& name, std::string const& description ) const {
		db::Connection::RowId result ;

		m_find_entity_statement
			->reset()
			.bind( 1, name ).step() ;

		if( m_find_entity_statement->empty() ) {
			m_insert_entity_statement
				->reset()
				.bind( 1, name )
				.bind( 2, description )
				.step() ;
				
			result = m_connection->get_last_insert_row_id() ;
		} else {
			result = m_find_entity_statement->get< db::Connection::RowId >( 0 ) ;
		}
		return result ;
	}

	db::Connection::RowId DBOutputter::get_or_create_entity_data( db::Connection::RowId const entity_id, db::Connection::RowId const variable_id, genfile::VariantEntry const& value ) const {
		db::Connection::RowId result ;

		m_find_entity_data_statement
			->reset()
			.bind( 1, entity_id )
			.bind( 2, variable_id ).step() ;

		if( m_find_entity_data_statement->empty() ) {
			m_insert_entity_data_statement
				->reset()
				.bind( 1, entity_id )
				.bind( 2, variable_id )
				.bind( 3, value )
				.step() ;
			result = m_connection->get_last_insert_row_id() ;
		} else {
			result = m_find_entity_data_statement->get< db::Connection::RowId >( 0 ) ;
		}
		return result ;
	}

	void DBOutputter::store_data(
		genfile::SNPIdentifyingData const& snp,
		std::string const& variable,
		genfile::VariantEntry const& value
	) {

		db::Connection::RowId snp_id = get_or_create_snp( snp ) ;
		db::Connection::RowId variable_id = get_or_create_variable( variable );

		assert( m_insert_summarydata_statement.get() ) ;
		m_insert_summarydata_statement
			->reset()
			.bind( 1, snp_id )
			.bind( 2, m_analysis_id )
			.bind( 3, variable_id )
			.bind( 4, value )
			.step()
		;
	}
}
