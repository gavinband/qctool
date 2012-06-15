
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef CALL_COMPARER_COMPONENT_CALL_COMPARER_DB_OUTPUTTER_HPP
#define CALL_COMPARER_COMPONENT_CALL_COMPARER_DB_OUTPUTTER_HPP

#include <string>
#include <memory>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <boost/thread/thread_time.hpp>
#include <boost/thread/thread.hpp>
#include "db/Connection.hpp"
#include "db/SQLStatement.hpp"
#include "qcdb/DBOutputter.hpp"
#include "components/CallComparerComponent/PairwiseCallComparerManager.hpp"

struct CallComparerDBOutputter: public qcdb::DBOutputter, public PairwiseCallComparerManager::ComparisonClient, public PairwiseCallComparerManager::MergeClient {
	typedef std::auto_ptr< CallComparerDBOutputter > UniquePtr ;
	typedef boost::shared_ptr< CallComparerDBOutputter > SharedPtr ;
	
	static UniquePtr create( std::string const& filename, std::string const& analysis, Metadata const& metadata = Metadata() ) ;
	static SharedPtr create_shared( std::string const& filename, std::string const& analysis, Metadata const& metadata = Metadata() ) ;

	CallComparerDBOutputter( std::string const& filename, std::string const& analysis, Metadata const& metadata = Metadata() ) ;
	~CallComparerDBOutputter() ;
	void begin_comparisons( genfile::SNPIdentifyingData const& snp ) ;
	void set_result(
		std::string const& callset1,
		std::string const& callset2,
		std::string const& comparison_method,
		std::string const& comparison_variable,
		genfile::VariantEntry const& value
	) ;
	void end_comparisons() ;

	void set_result(
		std::string const& comparison_method,
		std::string const& comparison_variable,
		genfile::VariantEntry const& value
	) ;

private:
	std::size_t const m_max_transaction_count ;
	db::Connection::RowId const m_callset_id ;
	db::Connection::StatementPtr m_insert_comparison_statement ;
	
	typedef std::vector< boost::tuple< genfile::SNPIdentifyingData, std::string, std::string, std::string, std::string, genfile::VariantEntry > > Data ;
	Data m_data ;

	genfile::SNPIdentifyingData m_snp ;

private:
	void construct_statements() ;
	void write_data( Data const& data ) ;
	void store_comparison(
		db::Connection::RowId const snp_id,
		std::string const& callset1,
		std::string const& callset2,
		std::string const& comparison_method,
		std::string const& comparison_variable,
		genfile::VariantEntry const& value
	) ;
} ;

#endif
