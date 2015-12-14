
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_SNP_SUMMARY_COMPONENT_DB_OUTPUTTER_HPP
#define QCTOOL_SNP_SUMMARY_COMPONENT_DB_OUTPUTTER_HPP

#include <string>
#include <memory>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/VariantEntry.hpp"
#include "db/Connection.hpp"
#include "db/SQLStatement.hpp"
#include "qcdb/DBOutputter.hpp"
#include "qcdb/Storage.hpp"

namespace snp_summary_component {
	struct DBOutputter: public qcdb::Storage {
		typedef qcdb::DBOutputter::Metadata Metadata ;
		typedef qcdb::Storage Storage ;
		typedef boost::shared_ptr< DBOutputter > SharedPtr ;

		static UniquePtr create( std::string const& filename, std::string const& analysis_name, std::string const& analysis_description, Metadata const& metadata ) ;
		static SharedPtr create_shared(
			std::string const& filename,
			std::string const& analysis_name,
			std::string const& analysis_description,
			Metadata const& metadata
		) ;

		DBOutputter(
			std::string const& filename,
			std::string const& analysis_name,
			std::string const& analysis_description,
			Metadata const& metadata
		) ;

		~DBOutputter() ;

		// Set the name of the table where results are stored
		void set_table_name( std::string const& table_name ) ;

		void add_variable(
			std::string const& 
		) ;

		void create_new_variant( genfile::VariantIdentifyingData const& ) ;
		void store_per_variant_data(
			genfile::VariantIdentifyingData const& snp,
			std::string const& variable,
			genfile::VariantEntry const& value
		) ;
		
		void finalise( long options ) ;

		db::Connection::RowId analysis_id() const ;
	private:
		qcdb::DBOutputter m_outputter ;
		std::size_t const m_max_transaction_count ;
		db::Connection::RowId m_variable_id ;
		typedef std::vector< boost::tuple< genfile::VariantIdentifyingData, std::string, genfile::VariantEntry > > Data ;
		Data m_data ;
		
		std::string m_table_name ;
		db::Connection::StatementPtr m_insert_data_statement ;

	private:
		void create_schema() ;
		void write_data( Data const& data ) ;

		void store_data(
			db::Connection::RowId const snp_id,
			std::string const& variable,
			genfile::VariantEntry const& value
		) ;
		
		void insert_summary_data(
			db::Connection::RowId snp_id,
			db::Connection::RowId variable_id,
			genfile::VariantEntry const& value
		) const ;
	} ;
}

#endif
