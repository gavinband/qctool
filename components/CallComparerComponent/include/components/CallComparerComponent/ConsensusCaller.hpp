#ifndef QCTOOL_CALLCOMPARER_COMPONENT_CONSENSUS_CALLER_HPP
#define QCTOOL_CALLCOMPARER_COMPONENT_CONSENSUS_CALLER_HPP

#include <map>
#include <string>
#include <Eigen/Core>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "components/CallComparerComponent/PairwiseCallComparerManager.hpp"

struct ConsensusCaller: PairwiseCallComparerManager::Client
{
	ConsensusCaller( std::string const& comparison_method, double threshhold ):
		m_comparison_method( comparison_method ),
		m_pvalue_threshhold( 0.0001 )
	{}
	
	void begin_comparisons( genfile::SNPIdentifyingData const& snp ) {
		m_comparison_values.clear() ;
	}

	void end_comparisons() {
		// Divide the calls into groups of calls all having p-value greater than the threshhold.
		// We assume comparisons are transitive: in reality this isn't the case but it's simpler to
		// assume transitivity here.

		Eigen::MatrixXd mismatches( m_comparison_values.size(), m_comparison_values.size() ) ;
		for( )
		
		std::vector< std::set< std::string > > groups ;
		std::map< std::string, std::size_t > group_indices ;
		ComparisonValues::const_iterator
			i = m_comparison_values.begin(),
			end_i = m_comparison_values.end() ;
		for( ; i != end_i; ++i ) {
			std::map< std::string, std::size_t >::iterator
				callset1_i = group_indices.find( i->first.first ),
				callset2_i = group_indices.find( i->first.second ) ;
			
		}
	}

	void send_results(
		std::string const& callset1,
		std::string const& callset2,
		std::string const& comparison_method,
		std::string const& comparison_variable,
		genfile::VariantEntry const& value
	) {
		if( comparison_method == m_comparison_method && comparison_variable == "pvalue" )) {
			m_comparison_values[ std::make_pair( callset1, callset2 ) ] = value ;
		}
	}
	
private:
	std::string m_comparison_method ;
	typedef std::map< std::pair< std::string, std::string >, genfile::VariantEntry > ComparisonValues ;
	ComparisonValues m_comparison_values ;
} ;

#endif

