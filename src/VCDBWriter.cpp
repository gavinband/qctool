#include <iostream>
#include <boost/bind.hpp>
#include "genfile/FileUtils.hpp"
#include "db/SQLite3Connection.hpp"
#include "db/SQLStatement.hpp"
#include "genfile/Error.hpp"
#include "genfile/zlib.hpp"
#include "genfile/endianness_utils.hpp"
#include "DataStore.hpp"
#include "VCDBWriter.hpp"

VCDBWriter::UniquePtr VCDBWriter::create( std::string const& filename ) {
	VCDBWriter::UniquePtr result(
		new VCDBWriter(
			DataStore::create(
				"sqlite3://" + filename
			)
		)
	) ;
	return result ;
}

VCDBWriter::VCDBWriter( DataStore::UniquePtr store ):
	m_store( store )
{
	setup() ;
}

void VCDBWriter::setup() {
	m_storage_id = m_store->get_or_create_entity( "qctool_serialized", "serialized using qctool." ) ;
	m_store->set_relationship( "qctool_serialized", "is_a", "storage_type" ) ;
	m_store->get_or_create_entity( "per_variant_per_sample_data", "Indicates data pertains to each sample at a variant" ) ;	
	m_cohort_id = m_store->get_or_create_entity( "unnamed cohort", "An unnamed cohort" ) ;	
}

db::Connection::RowId VCDBWriter::get_or_create_field( std::string const& field, std::string const& type ) {
	EntityCache::const_iterator where = m_entity_cache.find( field ) ;
	if( where != m_entity_cache.end() ) {
		return where->second ;
	}
	else {
		DataStore::EntityId field_id = m_store->get_or_create_entity( field, "" ) ;
		m_store->get_or_create_entity( type, "" ) ;
		m_store->set_relationship( field, "has_type", type ) ;
		m_store->set_relationship( field, "stored_as", "qctool_serialized" ) ;
		m_store->set_relationship( field, "is_a", "per_variant_per_sample_data" ) ;
		m_entity_cache[ field ] = field_id ;
		return field_id ;
	}
}

void VCDBWriter::begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) {
	m_number_of_samples = number_of_samples ;
	m_number_of_snps = number_of_snps ;
	m_number_of_snps_written = 0 ;
}

void VCDBWriter::set_SNPs( std::vector< genfile::SNPIdentifyingData > const& snps ) {
	std::size_t const CHUNK_SIZE = 1000 ;
	for( std::size_t chunk = 0; chunk < ( snps.size() / CHUNK_SIZE ); ++chunk ) {
		DataStore::TransactionPtr transaction = m_store->open_transaction() ;
		for( std::size_t i = chunk * CHUNK_SIZE; i < std::min( (chunk+1) * CHUNK_SIZE, snps.size() ); ++i ) {
			m_store->get_or_create_SNP( snps[i] ) ;
		}
	}
}

namespace impl {
	void add_spec_to_map( std::map< std::string, std::string >* map, std::string const& name, std::string const& type ) {
		(*map)[ name ] = type ;
	}
}

void VCDBWriter::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
	std::map< std::string, std::string > fields ;
	data_reader.get_supported_specs( boost::bind( &impl::add_spec_to_map, &fields, _1, _2 )) ;
	
	try {
		DataStore::TransactionPtr transaction = m_store->open_transaction() ;
		DataStore::EntityId snp_id  = m_store->get_or_create_SNP( snp ) ;

		for(
			std::map< std::string, std::string >::const_iterator field_i = fields.begin();
			field_i != fields.end();
		 	++field_i
 		) {
			std::string const field = field_i->first ;
			std::string const type = field_i->second ;
			DataStore::EntityId const field_id = get_or_create_field( field, type ) ;

			std::vector< std::vector< genfile::VariantEntry > >& data = m_data ;
			data.resize( m_number_of_samples ) ;
			data_reader.get( field, data ) ;
			assert( data.size() == m_number_of_samples ) ;

			// count the data
			std::size_t count = 0 ;
			for( std::size_t i = 0; i < data.size(); ++i ) {
				count += data[i].size() ;
			}

			std::vector< char >& buffer = m_buffer ;
			m_buffer.resize( 100 + (count * 8) ) ;
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
		
			m_store->store_per_variant_data(
				snp_id,
				field_id,
				m_cohort_id,
				m_storage_id,
				&buffer[0],
				&buffer[0] + buffer.size()
			) ;
		}
	}
	catch( db::StatementPreparationError const& e ) {
		std::cerr << "(" << e.what() << "): " << e.description() << ": preparing statement \"" + e.SQL() + "\".\n" ;
		throw genfile::OperationFailedError( "VCDBWriter::processed_snp()", m_store->get_spec(), "SQL insert" ) ;
	}
	catch( db::Error const& e ) {
		std::cerr << "(" << e.what() << "): " << e.description() << ".\n" ;
		throw genfile::OperationFailedError( "VCDBWriter::processed_snp()", m_store->get_spec(), "SQL insert" ) ;
	}
	++m_number_of_snps_written ;
}

void VCDBWriter::end_processing_snps() {
	
}



