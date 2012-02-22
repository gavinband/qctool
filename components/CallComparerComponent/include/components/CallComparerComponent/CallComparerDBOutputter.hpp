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
#include "components/CallComparerComponent/PairwiseCallComparerManager.hpp"

struct CallComparerDBOutputter: public PairwiseCallComparerManager::ComparisonClient, public PairwiseCallComparerManager::MergeClient {
	typedef std::auto_ptr< CallComparerDBOutputter > UniquePtr ;
	typedef boost::shared_ptr< CallComparerDBOutputter > SharedPtr ;
	
	static UniquePtr create( std::string const& filename ) ;
	static SharedPtr create_shared( std::string const& filename ) ;

	CallComparerDBOutputter( std::string const& filename ) ;

	
	~CallComparerDBOutputter() ;
	void begin_comparisons( genfile::SNPIdentifyingData const& snp ) ;
	void end_comparisons() ;
	void set_result(
		std::string const& callset1,
		std::string const& callset2,
		std::string const& comparison_method,
		std::string const& comparison_variable,
		genfile::VariantEntry const& value
	) ;

	void set_result(
		std::string const& comparison_method,
		std::string const& comparison_variable,
		genfile::VariantEntry const& value
	) ;

private:
	db::Connection::UniquePtr m_connection ;
	std::size_t const m_max_transaction_count ;

	db::Connection::StatementPtr m_find_variant_statement ;
	db::Connection::StatementPtr m_insert_variant_statement ;
	db::Connection::StatementPtr m_find_entity_statement ;
	db::Connection::StatementPtr m_insert_entity_statement ;
	db::Connection::StatementPtr m_insert_comparison_statement ;
	db::Connection::StatementPtr m_insert_summarydata_statement ;
	
	typedef std::vector< boost::tuple< genfile::SNPIdentifyingData, std::string, std::string, std::string, std::string, genfile::VariantEntry > > Data ;
	Data m_data ;

	genfile::SNPIdentifyingData m_snp ;

private:
	void construct_statements() ;
	void write_data( Data const& data ) ;
	void store_comparison(
		genfile::SNPIdentifyingData const& snp,
		std::string const& callset1,
		std::string const& callset2,
		std::string const& comparison_method,
		std::string const& comparison_variable,
		genfile::VariantEntry const& value
	) ;
} ;

#endif
