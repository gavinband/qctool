
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_SNP_SUMMARY_COMPONENT_FLAT_TABLE_DB_OUTPUTTER_HPP
#define QCTOOL_SNP_SUMMARY_COMPONENT_FLAT_TABLE_DB_OUTPUTTER_HPP

#include <string>
#include <memory>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/bimap.hpp>
#include "genfile/SNPIdentifyingData2.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/VariantEntry.hpp"
#include "db/Connection.hpp"
#include "db/SQLStatement.hpp"
#include "qcdb/DBOutputter.hpp"
#include "components/SNPSummaryComponent/Storage.hpp"

namespace snp_summary_component {
	struct FlatTableDBOutputter: public Storage {
		typedef std::auto_ptr< FlatTableDBOutputter > UniquePtr ;
		typedef boost::shared_ptr< FlatTableDBOutputter > SharedPtr ;
		typedef qcdb::DBOutputter::Metadata Metadata ;
		static UniquePtr create( std::string const& filename, std::string const& cohort_name, Metadata const& metadata ) ;
		static SharedPtr create_shared( std::string const& filename, std::string const& cohort_name, Metadata const& metadata ) ;

		FlatTableDBOutputter(
			std::string const& filename,
			std::string const& cohort_name,
			Metadata const& metadata
		) ;

		~FlatTableDBOutputter() ;

		void set_table_name( std::string const& table_name ) ;
		
		void store_per_variant_data(
			genfile::SNPIdentifyingData2 const& snp,
			std::string const& variable,
			genfile::VariantEntry const& value
		) ;
		
		void finalise() ;

	private:
		qcdb::DBOutputter m_outputter ;
		std::string m_table_name ;
		std::size_t const m_max_snps_per_block ;
		db::Connection::StatementPtr m_insert_data_sql ;
		std::vector< genfile::SNPIdentifyingData2 > m_snps ;
		typedef boost::bimap< std::string, std::size_t > VariableMap ;
		VariableMap m_variables ;
		typedef std::map< std::pair< std::size_t, std::size_t >, genfile::VariantEntry > ValueMap ;
		ValueMap m_values ;

	private:
		std::string get_table_name() const ;
		void store_block() ;
		void create_schema() ;
		void store_data_for_variant(
			std::size_t const,
			db::Connection::RowId const,
			db::Connection::RowId const
		) ;
		void create_variables() ;
	} ;
}

#endif
