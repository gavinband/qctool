#ifndef CALL_COMPARER_COMPONENT_CALL_COMPARER_FILE_OUTPUTTER_HPP
#define CALL_COMPARER_COMPONENT_CALL_COMPARER_FILE_OUTPUTTER_HPP

#include <memory>
#include <boost/shared_ptr.hpp>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "components/CallComparerComponent/PairwiseCallComparerManager.hpp"

struct CallComparerFileOutputter: public PairwiseCallComparerManager::ComparisonClient, public PairwiseCallComparerManager::MergeClient {
	typedef std::auto_ptr< CallComparerFileOutputter > UniquePtr ;
	typedef boost::shared_ptr< CallComparerFileOutputter > SharedPtr ;
	
	static UniquePtr create( std::string const& filename, std::string const& analysis ) ;
	static SharedPtr create_shared( std::string const& filename, std::string const& analysis ) ;

	CallComparerFileOutputter( std::string const& filename ) ;

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
	std::string const m_filename ;
	statfile::BuiltInTypeStatSink::UniquePtr m_sink ;
	genfile::SNPIdentifyingData m_snp ;
} ;

#endif
