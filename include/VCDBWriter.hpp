#ifndef QCTOOL_VCDB_WRITER_HPP
#define QCTOOL_VCDB_WRITER_HPP

#include <vector>
#include <set>
#include <map>
#include <string>
#include <boost/noncopyable.hpp>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "db/SQLite3Connection.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "DataStore.hpp"

class VCDBWriter: public genfile::SNPDataSourceProcessor::Callback, public boost::noncopyable {
public:
	static void declare_options( appcontext::OptionProcessor& options ) ;

	typedef std::auto_ptr< VCDBWriter > UniquePtr ;
	static UniquePtr create( appcontext::OptionProcessor const& options ) ;
	
public:
	
	VCDBWriter( DataStore::UniquePtr store, std::vector< std::string > const& fields ) ;
	void begin_processing_snps( std::size_t number_of_samples ) ;
	void processed_snp( genfile::SNPIdentifyingData const& , genfile::VariantDataReader& data_reader ) ;
	void end_processing_snps() ;
	void set_SNPs( std::vector< genfile::SNPIdentifyingData > const& snps ) ;

private:
	DataStore::UniquePtr m_store ;
	std::set< std::string > const m_fields ;
	bool m_write_all_fields ;
	std::size_t m_number_of_samples ;
	std::size_t m_number_of_snps_written ;
	std::size_t m_transaction_count ;
	std::size_t m_transaction_limit ;
	DataStore::TransactionPtr m_transaction ;
	
	DataStore::EntityId m_cohort_id ;
	DataStore::EntityId m_storage_id ;

	typedef std::map< std::string, db::Connection::RowId > EntityCache ;
	EntityCache m_entity_cache ;

	std::vector< std::vector< genfile::VariantEntry > > m_data ;
	std::vector< char > m_buffer ;

private:
	void setup() ;
	db::Connection::RowId get_or_create_field( std::string const& field, std::string const& type ) ;
	void store_field( DataStore::EntityId const snp_id, std::string const& field, std::string const& type, genfile::VariantDataReader& data_reader ) ;
	void unsafe_store_field( DataStore::EntityId const snp_id, std::string const& field, std::string const& type, genfile::VariantDataReader& data_reader ) ;
} ;

#endif

