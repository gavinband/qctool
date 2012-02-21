#ifndef QCTOOL_CALL_COMPARER_COMPONENT_HPP
#define QCTOOL_CALL_COMPARER_COMPONENT_HPP

#include <boost/noncopyable.hpp>
#include "PairwiseCallComparerManager.hpp"

struct CallComparerComponent: public genfile::SNPDataSourceProcessor::Callback, public virtual boost::noncopyable
{
public:
	typedef std::auto_ptr< CallComparerComponent > UniquePtr ;
	static UniquePtr create( appcontext::OptionProcessor const& options ) ;
	static UniquePtr create( PairwiseCallComparerManager::UniquePtr call_comparer, std::vector< std::string > const& call_fields ) ;
	static void declare_options( appcontext::OptionProcessor& options ) ;

public:
	CallComparerComponent( PairwiseCallComparerManager::UniquePtr call_comparer, std::vector< std::string > const& call_fields )  ;

	PairwiseCallComparerManager& get_comparer() { return *m_call_comparer ; }

public:
	void add_callbacks( genfile::SNPDataSourceProcessor& processor ) ;
	void begin_processing_snps( std::size_t number_of_samples ) ;
	void processed_snp( genfile::SNPIdentifyingData const&, genfile::VariantDataReader& data_reader ) ;
	void end_processing_snps() ;
	
private:
	PairwiseCallComparerManager::UniquePtr m_call_comparer ;
	std::vector< std::string > m_call_fields ;
	PairwiseCallComparerManager::MergeClient::UniquePtr m_consensus_caller ;
} ;


#endif
