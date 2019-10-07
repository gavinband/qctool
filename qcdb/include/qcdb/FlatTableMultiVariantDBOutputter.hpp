
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_QCDB_FLAT_TABLE_MULTI_VARIANT_OUTPUTTER_HPP
#define QCTOOL_QCDB_FLAT_TABLE_MULTI_VARIANT_OUTPUTTER_HPP

#include <string>
#include <memory>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/bimap.hpp>
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/db/Connection.hpp"
#include "genfile/db/SQLStatement.hpp"
#include "qcdb/MultiVariantStorage.hpp"
#include "qcdb/DBOutputter.hpp"

namespace qcdb {
	struct FlatTableMultiVariantDBOutputter: public MultiVariantStorage {
		typedef std::auto_ptr< FlatTableMultiVariantDBOutputter > UniquePtr ;
		typedef boost::shared_ptr< FlatTableMultiVariantDBOutputter > SharedPtr ;
		typedef MultiVariantStorage::Metadata Metadata ;
		typedef MultiVariantStorage::AnalysisId AnalysisId ;
		typedef std::vector< genfile::VariantIdentifyingData > Key ;
			
		static UniquePtr create(
			std::string const& filename,
			std::size_t const number_of_key_fields,
			std::string const& analysis_name,
			std::string const& analysis_description,
			Metadata const& metadata,
			std::string const& snp_match_fields = "position,alleles",
			boost::optional< genfile::db::Connection::RowId > analysis_id = boost::optional< genfile::db::Connection::RowId >()
		) ;

		FlatTableMultiVariantDBOutputter(
			std::string const& filename,
			std::size_t number_of_key_fields,
			std::string const& analysis_name,
			std::string const& analysis_description,
			Metadata const& metadata,
			std::string const& snp_match_fields = "position,alleles",
			boost::optional< genfile::db::Connection::RowId > analysis_id = boost::optional< genfile::db::Connection::RowId >()
		) ;

		~FlatTableMultiVariantDBOutputter() ;

		void set_table_name( std::string const& table_name ) ;
		void set_without_rowid() ;

		void set_variant_names( std::vector< std::string > const& names ) ;
		void add_variable( std::string const& ) ;

		void create_new_key( Key const& key ) ;

		void store_data_for_key(
			Key const& key,
			std::string const& variable,
			genfile::VariantEntry const& value
		) ;
		
		void finalise( long options = qcdb::eCreateIndices ) ;

		AnalysisId analysis_id() const ;

	private:
		qcdb::DBOutputter m_outputter ;
		std::vector< std::string > m_key_entry_names ;
		std::string m_table_name ;
		bool m_without_rowid ;
		std::size_t const m_max_variants_per_block ;
		genfile::db::Connection::StatementPtr m_insert_data_sql ;
		std::vector< Key > m_keys ;
		typedef boost::bimap< std::string, std::size_t > VariableMap ;
		VariableMap m_variables ;
		typedef std::map< std::pair< std::size_t, std::size_t >, genfile::VariantEntry > ValueMap ;
		ValueMap m_values ;

	private:
		std::string get_table_name() const ;
		void store_block() ;
		void create_schema() ;
		void store_data_for_variants(
			std::size_t const key_i,
			genfile::db::Connection::RowId const analysis_id,
			std::vector< genfile::db::Connection::RowId > const& variant_ids
		) ;
	} ;
}
#endif
