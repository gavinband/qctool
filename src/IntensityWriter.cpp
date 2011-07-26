#include <iostream>
#include <boost/bind.hpp>
#include "genfile/FileUtils.hpp"
#include "db/SQLite3Connection.hpp"
#include "db/SQLStatement.hpp"
#include "genfile/Error.hpp"
#include "genfile/zlib.hpp"
#include "genfile/endianness_utils.hpp"
#include "IntensityWriter.hpp"

IntensityWriter::IntensityWriter( std::string const& filename ):
	m_filename( filename ),
	m_connection( new db::SQLite3Connection( filename ))
{
	setup( *m_connection ) ;
}

void IntensityWriter::setup( db::Connection& connection ) {
	db::Connection::StatementPtr statement ;
	connection.run_statement(
		"CREATE TABLE IF NOT EXISTS FileInfo ( key TEXT NOT NULL UNIQUE, value TEXT );"
	) ;
	statement = connection.get_statement(
		"INSERT OR REPLACE INTO FileInfo VALUES( ?1, ?2 )"
	) ;
	statement->bind( 1, "format" ) ;
	statement->bind( 2, "VCDBv0.1" ) ;
	statement->step() ;

	connection.run_statement(
		"CREATE TABLE IF NOT EXISTS Meta ( id INTEGER PRIMARY KEY, name TEXT NOT NULL UNIQUE, description TEXT ); "
		"CREATE INDEX IF NOT EXISTS Meta_name ON Meta( name )"
	) ;

	statement = connection.get_statement(
		"INSERT OR REPLACE INTO Meta VALUES( ?1, ?2, ?3 ) ;"
	) ;

	statement->bind( 1, 1 ) ;
	statement->bind( 2, "zlib_compressed_serialized" ) ;
	statement->bind( 3, "zlib compressed data, serialized using qctool." ) ;
	statement->step() ;

	connection.run_statement(
		"CREATE TABLE IF NOT EXISTS SNP ( "
		"id INTEGER PRIMARY KEY, "
		"rsid TEXT NOT NULL, "
		"chromosome TEXT NOT NULL, "
		"position INTEGER NOT NULL, "
		"alleleA TEXT NOT NULL, "
		"alleleB TEXT NOT NULL, "
		"UNIQUE( rsid, chromosome, position ) "
		")"
	) ;

	connection.run_statement(
		"CREATE INDEX IF NOT EXISTS SNP_rsid ON SNP( rsid )"
	) ;
	statement = connection.get_statement(
		"CREATE INDEX IF NOT EXISTS SNP_position ON SNP( chromosome, position )"
	) ;
	statement->step() ;
	
	statement = connection.get_statement(
		"CREATE TABLE IF NOT EXISTS Data ( "
			"snp_id INTEGER NOT NULL, "
			"field_id INTEGER NOT NULL, "
			"storage_id INTEGER NOT NULL, "
			"uncompressed_size INTEGER NOT NULL, "
			"data BLOB, "
			"UNIQUE( snp_id, field_id, storage_id )"
			"FOREIGN KEY (field_id) REFERENCES Field( id ), "
			"FOREIGN KEY (snp_id) REFERENCES SNP( id ), "
			"FOREIGN KEY (storage_id) REFERENCES Storage( id ) "
			") ;"
	) ;
	statement->step() ;
	statement = connection.get_statement(
		"CREATE INDEX IF NOT EXISTS Data_snp ON Data( snp_id )"
	) ;
	statement->step() ;
	statement = connection.get_statement(
		"CREATE INDEX IF NOT EXISTS Data_field ON Data( field_id )"
	) ;
	statement->step() ;
	statement = connection.get_statement(
		"CREATE INDEX IF NOT EXISTS Data_storage ON Data( storage_id )"
	) ;

	statement->step() ;
}

void IntensityWriter::begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) {
	m_number_of_samples = number_of_samples ;
	m_number_of_snps = number_of_snps ;
	m_number_of_snps_written = 0 ;
}

void IntensityWriter::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
	std::vector< std::string > fields ;
	data_reader.get_supported_specs( boost::bind( &std::vector< std::string >::push_back, &fields, _1 )) ;
	
	try {
		db::Connection::StatementPtr statement ;
		db::Connection::StatementPtr transaction = m_connection->get_statement(
			"BEGIN"
		) ;
		transaction->step() ;


		db::Connection::RowId snp_row_id ;
		{
			statement  = m_connection->get_statement(
				"SELECT id, rsid, chromosome, position, alleleA, alleleB FROM SNP WHERE rsid = ?1 AND chromosome = ?2 AND position = ?3"
			) ;
			statement->bind( 1, snp.get_rsid() ) ;
			statement->bind( 2, genfile::string_utils::to_string( snp.get_position().chromosome() )) ;
			statement->bind( 3, snp.get_position().position() ) ;
			statement->step() ;
			if( statement->empty() ) {
				statement  = m_connection->get_statement(
					"INSERT INTO SNP( rsid, chromosome, position, alleleA, alleleB ) VALUES( ?1, ?2, ?3, ?4, ?5 )"
				) ;
				statement->bind( 1, snp.get_rsid() ) ;
				statement->bind( 2, genfile::string_utils::to_string( snp.get_position().chromosome() )) ;
				statement->bind( 3, snp.get_position().position() ) ;
				statement->bind( 4, std::string( 1, snp.get_first_allele() ) ) ;
				statement->bind( 5, std::string( 1, snp.get_second_allele() ) ) ;
				statement->step() ;
				snp_row_id = m_connection->get_last_insert_row_id() ;
			}
			else {
				snp_row_id = statement->get_column< db::Connection::RowId >( 0 ) ;
				std::string alleleA = statement->get_column< std::string >( 4 ) ;
				std::string alleleB = statement->get_column< std::string >( 5 ) ;
				if( alleleA != std::string( 1, snp.get_first_allele() ) || alleleB != std::string( 1, snp.get_second_allele() )) {
					throw genfile::MismatchError(
						"IntensityWriter::processed_snp()", m_filename + ":SNP",
						genfile::string_utils::to_string( snp ),
						"entry with id " + genfile::string_utils::to_string( snp_row_id )
					) ;
				}
			}
		}
		// If we get here, the right SNP is present and we have its snp_row_id.

		for( std::size_t field_i = 0; field_i < fields.size(); ++field_i ) {
			std::string const field = fields[ field_i ] ;
			// Make sure we've got these fields in Meta
			statement = m_connection->get_statement(
				"SELECT id FROM Meta WHERE name == ?1"
			) ;
			statement->bind( 1, field ) ;
			statement->step() ;
			if( statement->empty() ) {
				statement = m_connection->get_statement(
					"INSERT INTO Meta( name ) VALUES( ?1 )"
				) ;
				statement->bind( 1, field ) ;
				statement->step() ;

				statement = m_connection->get_statement(
					"SELECT id FROM Meta WHERE name == ?1"
				) ;
				statement->bind( 1, field ) ;
				statement->step() ;
			}

			db::Connection::RowId field_id = statement->get_column< int >( 0 ) ;
			statement->step() ;
			if( !statement->empty() ) {
				throw genfile::DuplicateKeyError( m_filename + ":Meta", "name=\"" + field + "\"" ) ;
			}
			
			// Look to see if the data is there already.
			statement = m_connection->get_statement(
				"SELECT snp_id FROM Data WHERE snp_id == ?1 AND field_id == ?2"
			) ;
			statement->bind( 1, snp_row_id ) ;
			statement->bind( 2, field_id ) ;
			statement->step() ;
			if( statement->empty() ) {
				// std::cerr << "SNP " << snp_row_id << ", field_id " << field_id << ", statement is empty.\n" ;
				// No data already.
				// Compress the data and store it.
				std::vector< std::vector< genfile::VariantEntry > > data( m_number_of_samples ) ;
				data_reader.get( field, genfile::VariantDataReader::set( data ) ) ;
				assert( data.size() == m_number_of_samples ) ;

				// count the data
				std::size_t count = 0 ;
				for( std::size_t i = 0; i < data.size(); ++i ) {
					count += data[i].size() ;
				}

				std::vector< char > buffer( count * 8 ) ;
				{
					char* begin = &buffer[0] ;
					char* end = begin + buffer.size() ;
					begin = genfile::write_small_integer( begin, end, uint64_t( m_number_of_samples ) ) ;
					for( std::size_t i = 0; i < data.size(); ++i ) {
						begin = genfile::write_small_integer( begin, end, data[i].size() ) ;
						for( std::size_t j = 0; j < data[i].size(); ++j ) {
							if( end < begin + 100 ) {
								std::size_t index = begin - &buffer[0] ;
								buffer.resize( buffer.size() + 1000 ) ;
								begin = &buffer[0] + index ;
								end = &buffer[0] + buffer.size() ;
							}
							begin = data[i][j].serialize( begin, end ) ;
						}
					}
					buffer.resize( begin - &buffer[0] ) ;
				}
			
				std::vector< char > compressed_buffer ;
				genfile::zlib_compress( buffer, &compressed_buffer ) ;
				db::Connection::StatementPtr statement = m_connection->get_statement(
					"INSERT OR REPLACE INTO Data VALUES( ?1, ?2, ?3, ?4, ?5 )"
				) ;
				statement->bind( 1, snp_row_id ) ;
				statement->bind( 2, field_id ) ;
				statement->bind( 3, 1 ) ; // zlib-compressed, serialised.
				statement->bind( 4, uint64_t( buffer.size() ) ) ;
				statement->bind( 5, &compressed_buffer[0], compressed_buffer.size() ) ;
				statement->step() ;
			}
		}
		transaction = m_connection->get_statement(
			"COMMIT"
		) ;
		transaction->step() ;
	}
	catch( db::SQLite3Error const& e ) {
		std::cerr << "(" << e.what() << "): " << e.description() << ".\n" ;
		throw genfile::OperationFailedError( "IntensityWriter::processed_snp()", "file://" + m_filename, "SQL insert" ) ;
	}
	++m_number_of_snps_written ;
}

void IntensityWriter::end_processing_snps() {
	
}



