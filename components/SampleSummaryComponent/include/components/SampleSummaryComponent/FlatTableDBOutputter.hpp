
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_SAMPLE_SUMMARY_COMPONENT_FLAT_TABLE_DB_OUTPUTTER_HPP
#define QCTOOL_SAMPLE_SUMMARY_COMPONENT_FLAT_TABLE_DB_OUTPUTTER_HPP

#include <string>
#include <memory>
#include <map>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/bimap.hpp>
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "db/Connection.hpp"
#include "db/SQLStatement.hpp"
#include "qcdb/DBOutputter.hpp"
#include "components/SampleSummaryComponent/SampleStorage.hpp"
#include "qcdb/StorageOptions.hpp"

namespace sample_stats {
	struct FlatTableDBOutputter: public SampleStorage {
		typedef std::auto_ptr< FlatTableDBOutputter > UniquePtr ;
		typedef boost::shared_ptr< FlatTableDBOutputter > SharedPtr ;

		static UniquePtr create(
			std::string const& filename,
			std::string const& analysis_name,
			std::string const& analysis_description,
			qcdb::DBOutputter::Metadata const& metadata,
			genfile::CohortIndividualSource const& samples,
			boost::optional< db::Connection::RowId > analysis_id
		) ;
		static SharedPtr create_shared(
			std::string const& filename,
			std::string const& analysis_name,
			std::string const& analysis_description,
			qcdb::DBOutputter::Metadata const& metadata,
			genfile::CohortIndividualSource const& samples,
			boost::optional< db::Connection::RowId > analysis_id
		) ;

		FlatTableDBOutputter(
			std::string const& filename,
			std::string const& analysis_name,
			std::string const& analysis_description,
			qcdb::DBOutputter::Metadata const& metadata,
			genfile::CohortIndividualSource const& samples,
			boost::optional< db::Connection::RowId > analysis_id
		) ;
		~FlatTableDBOutputter() ;

		void store_per_sample_data(
			std::string const& computation_name,
			std::size_t sample,
			std::string const& variable,
			std::string const& description,
			genfile::VariantEntry const& value
		) ;
		
		void set_table_name( std::string const& ) ;

		void finalise( long options = qcdb::eCreateIndices ) ;
		
		void add_variable( std::string const& ) ;
		
		AnalysisId analysis_id() const ;

	private:
		qcdb::DBOutputter m_outputter ;
		std::string m_table_name ;
		genfile::CohortIndividualSource const& m_samples ;
		std::size_t const m_max_transaction_count ;
		db::Connection::RowId m_variable_id ;
		db::Connection::StatementPtr m_find_sample_statement ;
		db::Connection::StatementPtr m_insert_sample_statement ;
		db::Connection::StatementPtr m_insert_data_sql ;

		std::vector< genfile::VariantEntry > m_sample_ids ;
		typedef std::vector< boost::tuple< std::string, std::size_t, std::string, std::string, genfile::VariantEntry > > Data ;
		Data m_data ;

		typedef boost::bimap< std::string, std::size_t > VariableMap ;
		VariableMap m_variables ;
		typedef std::map< std::pair< std::size_t, std::size_t >, genfile::VariantEntry > ValueMap ;
		ValueMap m_values ;

	private:
		VariableMap::left_const_iterator add_variable_impl( std::string const& ) ;
		void create_schema() ;
		void create_variables() ;
		void construct_statements() ;
		void store_samples( genfile::CohortIndividualSource const& samples ) ;
		db::Connection::RowId get_or_create_sample( genfile::VariantEntry const& identifier, std::size_t index ) const ;
		void store_block() ;
		void store_data_for_sample( db::Connection::RowId analysis_id, std::size_t sample_index, genfile::VariantEntry const& sample_id ) ;
	} ;
}

#endif
