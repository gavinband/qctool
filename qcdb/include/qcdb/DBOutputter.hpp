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

		static UniquePtr create( std::string const& filename, std::string const& cohort_name, std::string const& data_source_spec, std::string const& exclusions_name ) ;
		static SharedPtr create_shared( std::string const& filename, std::string const& cohort_name, std::string const& data_source_spec, std::string const& exclusions_name ) ;

		DBOutputter( std::string const& filename, std::string const& cohort_name, std::string const& data_source_spec, std::string const& exclusions_name ) ;
		~DBOutputter() ;

		db::Connection::RowId get_or_create_entity( std::string const& name, std::string const& description ) const ;
		db::Connection::RowId get_or_create_entity_data( db::Connection::RowId const entity_id, db::Connection::RowId const variable_id, genfile::VariantEntry const& value ) const ;
		void insert_summary_data( db::Connection::RowId snp_id, db::Connection::RowId variable_id, genfile::VariantEntry const& value ) const ;

	private:
		db::Connection::UniquePtr m_connection ;
		std::string const m_analysis_name ;
		std::string const m_source_spec ;
		std::string const m_exclusions_name ;

		db::Connection::StatementPtr m_find_entity_statement ;
		db::Connection::StatementPtr m_find_entity_data_statement ;
		db::Connection::StatementPtr m_insert_entity_statement ;
		db::Connection::StatementPtr m_insert_entity_data_statement ;
		db::Connection::StatementPtr m_insert_summarydata_statement ;
		db::Connection::RowId m_analysis_id ;

	protected:
		db::Connection& connection() const { return *m_connection ; }

	private:
		void construct_statements() ;

	} ;
}

#endif
