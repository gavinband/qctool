#include <iostream>
#include <boost/bind.hpp>
#include "genfile/FileUtils.hpp"
#include "db/SQLite3Connection.hpp"
#include "db/SQLStatement.hpp"
#include "genfile/Error.hpp"
#include "genfile/zlib.hpp"
#include "IntensityWriter.hpp"

IntensityWriter::IntensityWriter( std::string const& filename ):
	m_filename( filename ),
	m_connection( new db::SQLite3Connection( filename ))
{
	setup( *m_connection ) ;
}

void IntensityWriter::setup( db::Connection& connection ) {
	db::Connection::StatementPtr statement ;
	statement = connection.get_statement(
		"CREATE TABLE IF NOT EXISTS Meta ( id INTEGER PRIMARY KEY, name TEXT NOT NULL UNIQUE, description TEXT )"
	) ;
	statement->step() ;
	statement = connection.get_statement(
		"INSERT OR REPLACE INTO Meta VALUES( ?1, ?2, ?3 ) ;"
	) ;
	statement->bind( 1, 1 ) ;
	statement->bind( 2, "no compression" ) ;
	statement->bind( 3, "Indicates that the data is not compressed" ) ;
	statement->step() ;

	statement->reset() ;
	statement->bind( 1, 2 ) ;
	statement->bind( 2, "zlib compression" ) ;
	statement->bind( 3, "Indicates that the data is compressed with zlib" ) ;
	statement->step() ;

	statement->reset() ;
	statement->bind( 1, 3 ) ;
	statement->bind( 2, "bzip2 compression" ) ;
	statement->bind( 3, "Indicates that the data is compressed with libbzip2" ) ;
	statement->step() ;

	statement = connection.get_statement(
		"CREATE TABLE IF NOT EXISTS SNP ( id INTEGER PRIMARY KEY, rsid TEXT NOT NULL, chromosome TEXT NOT NULL, position INTEGER NOT NULL, alleleA TEXT NOT NULL, alleleB TEXT NOT NULL )"
	) ;
	statement->step() ;
	statement = connection.get_statement(
		"DELETE FROM SNP"
	) ;
	statement->step() ;
	statement = connection.get_statement(
		"CREATE TABLE IF NOT EXISTS Data ( "
			"snp_id INTEGER NOT NULL, "
			"field_id INTEGER NOT NULL, "
			"storage_id INTEGER NOT NULL, "
			"uncompressed_size INTEGER NOT NULL, "
			"data BLOB, "
			"FOREIGN KEY (field_id) REFERENCES Field( id ), "
			"FOREIGN KEY (snp_id) REFERENCES SNP( id ), "
			"FOREIGN KEY (storage_id) REFERENCES Storage( id ) "
			") ;"
	) ;
	statement->step() ;
	statement = connection.get_statement(
		"DELETE FROM Data"
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

		{
			statement  = m_connection->get_statement(
				"INSERT INTO SNP VALUES( ?1, ?2, ?3, ?4, ?5, ?6 )"
			) ;
			statement->bind( 1, m_number_of_snps_written ) ;
			statement->bind( 2, snp.get_rsid() ) ;
			statement->bind( 3, genfile::string_utils::to_string( snp.get_position().chromosome() )) ;
			statement->bind( 4, snp.get_position().position() ) ;
			statement->bind( 5, std::string( 1, snp.get_first_allele() ) ) ;
			statement->bind( 6, std::string( 1, snp.get_second_allele() ) ) ;
			statement->step() ;
		}
		
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

			std::size_t meta_id = statement->get_column< int >( 0 ) ;
			statement->step() ;
			if( !statement->empty() ) {
				throw genfile::DuplicateKeyError( m_filename + ":Meta", "name=\"" + field + "\"" ) ;
			}
			
			std::vector< std::vector< genfile::VariantEntry > > data( m_number_of_samples ) ;
			data_reader.get( field, genfile::VariantDataReader::set( data ) ) ;
			assert( data.size() == m_number_of_samples ) ;

			// count the data.
			std::vector< double > values ;
			values.reserve( m_number_of_samples * 3 ) ;
			for( std::size_t i = 0; i < m_number_of_samples; ++i ) {
				for( std::size_t j = 0; j < data[i].size(); ++j ) {
					if( data[i][j].is_missing() ) {
						values.push_back( std::numeric_limits< double >::quiet_NaN() ) ;
					}
					else {
						values.push_back( data[i][j].as< double >() ) ;
					}
				}
			}

			std::vector< char > buffer ;
			genfile::zlib_compress( values, &buffer ) ;
			db::Connection::StatementPtr statement = m_connection->get_statement(
				"INSERT INTO Data VALUES( ?1, ?2, ?3, ?4, ?5 )"
			) ;
			statement->bind( 1, m_number_of_snps_written ) ;
			statement->bind( 2, meta_id ) ;
			statement->bind( 3, 2 ) ; // gzip compressed data.
			statement->bind( 4, buffer.size() ) ;
			statement->bind( 5, &buffer[0], buffer.size() ) ;
			statement->step() ;
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



