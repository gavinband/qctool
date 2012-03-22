#ifndef QCTOOL_QCDB_DB_OUTPUTTER_HPP
#define QCTOOL_QCDB_DB_OUTPUTTER_HPP

#include <string>
#include <memory>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/VariantEntry.hpp"
#include "db/Connection.hpp"
#include "db/SQLStatement.hpp"

namespace qcdb {
	struct DBOutputter {
		typedef std::auto_ptr< DBOutputter > UniquePtr ;
		typedef boost::shared_ptr< DBOutputter > SharedPtr ;
		typedef std::map< std::string, std::pair< std::vector< std::string >, std::string > > Metadata ;

		static UniquePtr create( std::string const& filename, std::string const& cohort_name, Metadata const& metadata ) ;
		static SharedPtr create_shared( std::string const& filename, std::string const& cohort_name, Metadata const& metadata ) ;

		DBOutputter( std::string const& filename, std::string const& cohort_name, Metadata const& metadata ) ;
		~DBOutputter() ;

		db::Connection::RowId get_or_create_entity( std::string const& name, std::string const& description ) const ;
		db::Connection::RowId get_or_create_entity_data( db::Connection::RowId const entity_id, db::Connection::RowId const variable_id, genfile::VariantEntry const& value ) const ;
		db::Connection::RowId get_or_create_variant( genfile::SNPIdentifyingData const& snp ) const ;
		void insert_summary_data( db::Connection::RowId snp_id, db::Connection::RowId variable_id, genfile::VariantEntry const& value ) const ;

	private:
		db::Connection::UniquePtr m_connection ;
		std::string const m_analysis_name ;
		Metadata const m_metadata ;

		db::Connection::StatementPtr m_find_entity_statement ;
		db::Connection::StatementPtr m_find_entity_data_statement ;
		db::Connection::StatementPtr m_insert_entity_statement ;
		db::Connection::StatementPtr m_insert_entity_data_statement ;
		db::Connection::StatementPtr m_find_variant_statement ;
		db::Connection::StatementPtr m_insert_variant_statement ;
		db::Connection::StatementPtr m_insert_summarydata_statement ;
		db::Connection::RowId m_analysis_id ;

	protected:
		db::Connection& connection() const { return *m_connection ; }
		db::Connection::RowId analysis_id() const { return m_analysis_id ; }
	private:
		void construct_statements() ;
		void store_metadata() ;
	} ;
}

#endif
