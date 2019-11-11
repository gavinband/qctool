
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <memory>
#include <boost/optional.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/thread/thread.hpp>
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/Error.hpp"
#include "genfile/db/Connection.hpp"
#include "genfile/db/SQLStatement.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "qcdb/DBOutputter.hpp"

namespace qcdb {
	namespace impl {
		bool get_match_rsid( std::string spec ) {
			typedef genfile::VariantIdentifyingData::CompareFields CompareFields;
			CompareFields comparer( spec ) ;
			std::vector<int> const fields = comparer.get_compared_fields() ;

			if(
				std::find(
					fields.begin(),
					fields.end(),
					static_cast<int>(CompareFields::ePosition)
				) == fields.end()
				|| std::find(
					fields.begin(),
					fields.end(),
					static_cast<int>(CompareFields::eAlleles)
				) == fields.end()
			) {
				throw genfile::BadArgumentError(
					"qcdb::impl::get_match_rsid()",
					"spec=\"" + spec + "\"",
					"Spec must include position and alleles."
				) ;
			}
		
			bool match_rsid = (
				std::find(
					fields.begin(),
					fields.end(),
					static_cast<int>(CompareFields::eRSID)
				) != fields.end()
			) || (
				std::find(
					fields.begin(),
					fields.end(),
					static_cast<int>(CompareFields::eIDs)
				) != fields.end()
			) ;
			return match_rsid ;
		}
	}

	DBOutputter::UniquePtr DBOutputter::create(
		std::string const& filename,
		std::string const& analysis_name,
		std::string const& analysis_description,
		Metadata const& metadata,
		boost::optional< genfile::db::Connection::RowId > analysis_id,
		std::string const& snp_match_fields
	) {
		return UniquePtr( new DBOutputter( filename, analysis_name, analysis_description, metadata, analysis_id, snp_match_fields ) ) ;
	}
	DBOutputter::SharedPtr DBOutputter::create_shared(
		std::string const& filename,
		std::string const& analysis_name,
		std::string const& analysis_description,
		Metadata const& metadata,
		boost::optional< genfile::db::Connection::RowId > analysis_id,
		std::string const& snp_match_fields
	) {
		return SharedPtr( new DBOutputter( filename, analysis_name, analysis_description, metadata, analysis_id, snp_match_fields ) ) ;
	}

	DBOutputter::~DBOutputter() {}

	DBOutputter::DBOutputter(
		std::string const& filename,
		std::string const& analysis_name,
		std::string const& analysis_description,
		Metadata const& metadata,
		boost::optional< genfile::db::Connection::RowId > analysis_id,
		std::string const& snp_match_fields
	):
		m_connection( genfile::db::Connection::create( filename )),
		m_analysis_name( analysis_name ),
		m_analysis_chunk( analysis_description ),
		m_metadata( metadata ),
		m_create_indices( true ),
		m_match_rsid( impl::get_match_rsid( snp_match_fields )),
		m_analysis_id( analysis_id )
	{
		try {
			m_connection->run_statement( "PRAGMA journal_mode = OFF" ) ;
			m_connection->run_statement( "PRAGMA synchronous = OFF" ) ;
		}
		catch( genfile::db::Error const& ) {
			std::cerr << "qcdb::DBOutputter::DBOutputter(): unable to set PRAGMA synchronous=OFF, is another connection using this database?" ;
		}

		genfile::db::Connection::ScopedTransactionPtr transaction = m_connection->open_transaction( 7200 ) ;

		m_connection->run_statement(
			"CREATE TABLE IF NOT EXISTS Variant ( id INTEGER PRIMARY KEY, rsid TEXT, chromosome TEXT, position INTEGER, alleleA TEXT, alleleB TEXT )"
		) ;
		m_connection->run_statement(
			"CREATE INDEX IF NOT EXISTS Variant_position_index ON Variant( chromosome, position )"
		) ;
		m_connection->run_statement(
			"CREATE TABLE IF NOT EXISTS VariantIdentifier ( variant_id INTEGER NOT NULL, identifier TEXT, FOREIGN KEY( variant_id ) REFERENCES Variant( id ) ) "
		) ;
		m_connection->run_statement(
			"CREATE INDEX IF NOT EXISTS VariantIdentifierIdentifierIndex ON VariantIdentifier( identifier )"
		) ;
		m_connection->run_statement(
			"CREATE TABLE IF NOT EXISTS Analysis ( "
				"id INTEGER PRIMARY KEY, "
				"name TEXT, "
				"chunk TEXT"
			")"
		) ;
		m_connection->run_statement(
			"CREATE TABLE IF NOT EXISTS AnalysisProperty( "
			"analysis_id INTEGER NOT NULL REFERENCES Analysis( id ), "
			"property TEXT NOT NULL, "
			"value TEXT, "
			"source TEXT"
			")"
		) ;
		m_connection->run_statement(
			"CREATE VIEW IF NOT EXISTS AnalysisPropertyView AS "
			"SELECT analysis_id, name, property, value "
			"FROM AnalysisProperty AP "
			"INNER JOIN Analysis A "
			"ON A.id = AP.analysis_id"
		) ;
		m_connection->run_statement(
			"CREATE TABLE IF NOT EXISTS AnalysisStatus ( "
				"analysis_id INTEGER NOT NULL REFERENCES Analysis( id ), "
				"started TEXT NOT NULL, "
				"completed TEXT, "
				"status TEXT NOT NULL "
			")"
		) ;
		m_connection->run_statement(
			"CREATE VIEW IF NOT EXISTS AnalysisStatusView AS "
			"SELECT analysis_id, name AS analysis, chunk, started, completed, status "
			"FROM AnalysisStatus AST "
			"INNER JOIN Analysis A "
			"ON A.id == AST.analysis_id"
		) ;

		m_connection->run_statement(
			"CREATE TABLE IF NOT EXISTS Variable ( "
				"`analysis_id` INTEGER NOT NULL REFERENCES Analysis( id ), "
				"`table` TEXT NOT NULL, "
				"`name` TEXT NOT NULL, "
				"`description` TEXT"
			")"
		) ;

		m_connection->run_statement(
			"CREATE VIEW IF NOT EXISTS VariantView AS "
			"SELECT          V.id AS id, V.rsid AS rsid, V.chromosome AS chromosome, V.position AS position, V.alleleA AS alleleA, V.alleleB AS alleleB, "
			"GROUP_CONCAT( VI.identifier ) AS alternate_identifier "
			"FROM Variant V "
			"LEFT OUTER JOIN VariantIdentifier VI "
			"  ON VI.variant_id = V.id "
			"GROUP BY V.id"
		) ;
		construct_statements() ;
		store_metadata() ;
	}

	void DBOutputter::finalise( long options ) {
		if( options & eCreateIndices ) {
			genfile::db::Connection::ScopedTransactionPtr transaction = m_connection->open_transaction( 7200 ) ;
			m_connection->run_statement(
				"CREATE INDEX IF NOT EXISTS Variant_rsid_index ON Variant( rsid )"
			) ;
			m_connection->run_statement(
				"CREATE INDEX IF NOT EXISTS VariantIdentifierVariantIndex ON VariantIdentifier( variant_id )"
			) ;
		}
		genfile::db::Connection::ScopedTransactionPtr transaction = m_connection->open_transaction( 7200 ) ;
		end_analysis( m_analysis_id.get() ) ;
	}

	void DBOutputter::construct_statements() {
		m_insert_variable_statement = m_connection->get_statement( "INSERT INTO Variable ( analysis_id, `table`, `name`, `description` ) VALUES( ?, ?, ?, ? )" ) ;
		m_find_analysis_statement = m_connection->get_statement( "SELECT id FROM Analysis WHERE name == ?1 AND chunk == ?2" ) ;
		m_insert_analysis_statement = m_connection->get_statement( "INSERT INTO Analysis( name, chunk ) VALUES ( ?1, ?2 )" ) ;
		m_insert_analysis_property_statement = m_connection->get_statement( "INSERT OR REPLACE INTO AnalysisProperty ( analysis_id, property, value, source ) VALUES ( ?1, ?2, ?3, ?4 )" ) ;

		{
			std::string find_variant_statement_sql = "SELECT id, rsid FROM Variant WHERE chromosome == ?1 AND position == ?2 AND alleleA = ?3 AND alleleB = ?4" ;
			if( m_match_rsid ) {
				find_variant_statement_sql += " AND rsid == ?5" ;
			}
			m_find_variant_statement = connection().get_statement( find_variant_statement_sql ) ;
		}

		m_insert_variant_statement = connection().get_statement(
			"INSERT INTO Variant ( rsid, chromosome, position, alleleA, alleleB ) "
			"VALUES( ?1, ?2, ?3, ?4, ?5 )"
		) ;
		m_find_variant_identifier_statement = m_connection->get_statement( "SELECT * FROM VariantIdentifier WHERE variant_id == ?1 AND identifier == ?2" ) ;
		m_insert_variant_identifier_statement = m_connection->get_statement( "INSERT INTO VariantIdentifier( variant_id, identifier ) VALUES ( ?1, ?2 )" ) ;
	}

	void DBOutputter::store_metadata() {
		try {
			if( m_analysis_id ) {
				genfile::db::Connection::StatementPtr find_analysis_stmt = m_connection->get_statement(
					"SELECT name, chunk FROM Analysis WHERE id == ?1"
				) ;
				find_analysis_stmt
					->bind( 1, m_analysis_id.get() )
					.step() ;

				if( find_analysis_stmt->empty() ) {
					create_analysis(
						m_analysis_id.get(),
						m_analysis_name,
						m_analysis_chunk
					) ;
				} else {
					std::string const existing_analysis_name = find_analysis_stmt->get_column< std::string >(0) ;
					std::string const existing_analysis_chunk = find_analysis_stmt->get_column< std::string >(1) ;
					if( existing_analysis_name != m_analysis_name ) {
						throw genfile::BadArgumentError(
							"qcdb::DBOutputter::store_metadata()",
							"m_analysis_id=" + genfile::string_utils::to_string( m_analysis_id.get() ),
							"Existing analysis has mismatching name (\""
							+ existing_analysis_name
							+ "\", exopected \""
							+ m_analysis_name
							+ "\")"
						) ;
					}
					if( existing_analysis_chunk != m_analysis_chunk) {
						throw genfile::BadArgumentError(
							"qcdb::DBOutputter::store_metadata()",
							"m_analysis_id=" + genfile::string_utils::to_string( m_analysis_id.get() ),
							"Existing analysis has mismatching chunk (\""
							+ existing_analysis_chunk
							+ "\", exopected \""
							+ m_analysis_chunk
							+ "\")"
						) ;
					}
				}
			} else {
				m_analysis_id = create_analysis(
					m_analysis_name,
					m_analysis_chunk
				) ;
			}
		} catch( genfile::db::StatementStepError const& e ) {
			throw genfile::BadArgumentError( "qcdb::DBOutputter::store_metadata()", "analysis_name=\"" + m_analysis_name + "\"", "An analysis with name \"" + m_analysis_name + "\" and chunk \"" + m_analysis_chunk + "\" already exists" ) ;
		}

		start_analysis( m_analysis_id.get() ) ;

		for( Metadata::const_iterator i = m_metadata.begin(); i != m_metadata.end(); ++i ) {
			set_analysis_property(
				m_analysis_id.get(),
				i->first,
				genfile::string_utils::join( i->second.first, "," ),
				i->second.second
			) ;
		}
	}
	
	void DBOutputter::create_variable( std::string const& table, std::string const& name ) const {
		m_insert_variable_statement
			->bind( 1, m_analysis_id.get() )
			.bind( 2, table )
			.bind( 3, name )
			.step() ;
		m_insert_variable_statement->reset() ;
	}
	
	void DBOutputter::start_analysis( genfile::db::Connection::RowId const analysis_id ) const {
		genfile::db::Connection::StatementPtr stmnt = m_connection->get_statement( "INSERT INTO AnalysisStatus( analysis_id, started, status ) VALUES( ?, ?, ? )" ) ;
		stmnt->bind( 1, analysis_id ) ;
		stmnt->bind( 2, appcontext::get_current_time_as_string() ) ;
		stmnt->bind( 3, "incomplete" ) ;
		stmnt->step() ;
	}

	void DBOutputter::end_analysis( genfile::db::Connection::RowId const analysis_id ) const {
		genfile::db::Connection::StatementPtr stmnt = m_connection->get_statement( "UPDATE AnalysisStatus SET completed = ?, status = ? WHERE analysis_id == ?" ) ;
		stmnt->bind( 1, appcontext::get_current_time_as_string() ) ;
		stmnt->bind( 2, "success" ) ;
		stmnt->bind( 3, analysis_id ) ;
		stmnt->step() ;
	}

	void DBOutputter::create_analysis( genfile::db::Connection::RowId id, std::string const& name, std::string const& chunk ) const {
		genfile::db::Connection::StatementPtr insert_analysis_stmt = m_connection->get_statement(
			"INSERT INTO Analysis( id, name, chunk ) VALUES ( ?, ?, ? )"
		) ;
		insert_analysis_stmt
			->bind( 1, id )
			.bind( 2, name )
			.bind( 3, chunk )
			.step() ;
	}

	genfile::db::Connection::RowId DBOutputter::create_analysis( std::string const& name, std::string const& description ) const {
		genfile::db::Connection::RowId result ;
		m_insert_analysis_statement
			->bind( 1, name )
			.bind( 2, description )
			.step() ;
			
		result = m_connection->get_last_insert_row_id() ;
		m_insert_analysis_statement->reset() ;

		return result ;
	}

	genfile::db::Connection::RowId DBOutputter::set_analysis_property(
		genfile::db::Connection::RowId const analysis_id,
		std::string const& property,
		genfile::VariantEntry const& value,
		std::string const& aux
	) const {
		genfile::db::Connection::RowId result ;

		m_insert_analysis_property_statement
			->bind( 1, analysis_id )
			.bind( 2, property )
			.bind( 3, value )
			.bind( 4, aux )
			.step() ;
		result = m_connection->get_last_insert_row_id() ;
		m_insert_analysis_property_statement->reset() ;
		return result ;
	}

	void DBOutputter::add_alternative_variant_identifier(
		genfile::db::Connection::RowId const variant_id,
		std::string const& identifier,
		std::string const& rsid
	) const {
		if( identifier != rsid  && identifier != "---" && identifier != "." ) {
			add_variant_identifier( variant_id, identifier ) ;
		}
	}

	void DBOutputter::add_variant_identifier( genfile::db::Connection::RowId const variant_id, std::string const& identifier ) const {
		m_find_variant_identifier_statement
			->bind( 1, variant_id )
			.bind( 2, identifier )
			.step() ;
		if( m_find_variant_identifier_statement->empty() ) {
			m_insert_variant_identifier_statement
				->bind( 1, variant_id )
				.bind( 2, identifier )
				.step() ;
			m_insert_variant_identifier_statement->reset() ;
		}
		m_find_variant_identifier_statement->reset() ;
	}

	genfile::db::Connection::RowId DBOutputter::get_or_create_variant( genfile::VariantIdentifyingData const& snp ) const {
		genfile::db::Connection::RowId result ;
		if( snp.get_position().chromosome().is_missing() ) {
			m_find_variant_statement->bind_NULL( 1 ) ;
		} else {
			m_find_variant_statement
				->bind( 1, std::string( snp.get_position().chromosome() ) ) ;
		}
		m_find_variant_statement
			->bind( 2, snp.get_position().position() )
			.bind( 3, snp.get_allele(0) )
			.bind( 4, snp.get_allele(1) ) ;
		if( m_match_rsid ) {
			m_find_variant_statement->bind( 5, snp.get_primary_id() ) ;
		}
		m_find_variant_statement->step() ;

		if( m_find_variant_statement->empty() ) {
			if( snp.get_position().chromosome().is_missing() ) {
				m_insert_variant_statement
					->bind_NULL( 2 ) ;
			} else {
				m_insert_variant_statement
					->bind( 2, std::string( snp.get_position().chromosome() ) ) ;
			}

			m_insert_variant_statement
				->bind( 1, snp.get_primary_id() )
				.bind( 3, snp.get_position().position() )
				.bind( 4, snp.get_allele(0))
				.bind( 5, snp.get_allele(1))
				.step()
			;

			result = connection().get_last_insert_row_id() ;
			m_insert_variant_statement->reset() ;
			snp.get_identifiers(
				boost::bind(
					&DBOutputter::add_alternative_variant_identifier,
					this,
					result,
					_1,
					snp.get_primary_id()
				),
				1
			) ;
		} else {
			result = m_find_variant_statement->get< genfile::db::Connection::RowId >( 0 ) ;
			std::string const rsid = m_find_variant_statement->get< std::string >( 1 ) ;
			add_alternative_variant_identifier( result, snp.get_primary_id(), rsid ) ;
			snp.get_identifiers(
				boost::bind(
					&DBOutputter::add_alternative_variant_identifier,
					this,
					result,
					_1,
					rsid
				),
				1
			) ;
		}
		m_find_variant_statement->reset() ;
		return result ;
	}
	
}
