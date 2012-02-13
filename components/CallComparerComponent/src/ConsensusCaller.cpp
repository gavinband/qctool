#include <map>
#include <string>
#include <boost/bimap.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <Eigen/Core>

#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "components/CallComparerComponent/PairwiseCallComparerManager.hpp"
#include "components/CallComparerComponent/ConsensusCaller.hpp"

#define CONSENSUS_CALLER_DEBUG 1

ConsensusCaller::SharedPtr ConsensusCaller::create_shared( std::string const& comparison_method, double threshhold ) {
	return SharedPtr( new ConsensusCaller( comparison_method, threshhold )) ;
}

ConsensusCaller::ConsensusCaller( std::string const& comparison_method, double threshhold ):
	m_comparison_method( comparison_method ),
	m_range( std::make_pair( threshhold, 100.0 ))
{}

std::string ConsensusCaller::get_spec() const {
	return m_spec ;
}

std::string ConsensusCaller::get_result_as_string() const {
	std::string result ;
	foreach( std::string const& call, m_concordant_calls ) {
		if( result.size() > 0 ) {
			result += "," ;
		}
		result += call ;
	}
	return result ;
}

void ConsensusCaller::begin_comparisons( genfile::SNPIdentifyingData const& snp ) {
	m_snp = snp ;
	m_comparison_values.clear() ;
	m_concordant_calls.clear() ;
}

void ConsensusCaller::set_result(
	std::string const& callset1,
	std::string const& callset2,
	std::string const& comparison_method,
	std::string const& comparison_variable,
	genfile::VariantEntry const& value
) {
	if( comparison_method == m_comparison_method && comparison_variable == "pvalue" ) {
		m_comparison_values[ std::make_pair( callset1, callset2 ) ] = value ;
	}
}

void ConsensusCaller::end_comparisons() {
	// Divide the calls into groups of calls all having p-value greater than the threshhold.
	// We assume comparisons are transitive: in reality this isn't the case but it's simpler to
	// assume transitivity here.

	boost::bimap< std::string, int > call_names ;
	{
		int index = 0 ;
		foreach( ComparisonValues::value_type const& key, m_comparison_values ) {
			if( call_names.left.find( key.first.first ) == call_names.left.end() ) {
					call_names.left.insert( std::make_pair( key.first.first, index++ ) );
			}
			if( call_names.left.find( key.first.second ) == call_names.left.end() ) {
				call_names.left.insert( std::make_pair( key.first.second, index++ ) );
			}
		}
	}
	
	Eigen::MatrixXd mismatches = Eigen::MatrixXd::Zero( call_names.size(), call_names.size() ) ;
	foreach( ComparisonValues::value_type value, m_comparison_values ) {
		double pvalue = value.second.as< double >() ;
		if( pvalue < m_range.first || pvalue > m_range.second ) {
			std::size_t const i = call_names.left.at( value.first.first ) ;
			std::size_t const j = call_names.left.at( value.first.second ) ;
			mismatches( i, j )++ ;
			mismatches( j, i )++ ;
		}
	}

	// At a really good SNP there are no mismatches (=> all rows are zero.)
	// At a quite good SNP there is only one mismatch (=> only one row is nonzero.)
	Eigen::VectorXd counts = mismatches.rowwise().sum() ;
	for( int i = 0; i < counts.size(); ++i ) {
		if( counts[i] == 0 ) {
			m_concordant_calls.insert( call_names.right.at( i ) ) ;
		}
	}
	
#if CONSENSUS_CALLER_DEBUG
	std::cerr << "############################\n" ;
	std::cerr << "SNP: " << m_snp << ":\n" ;

	foreach( ComparisonValues::value_type value, m_comparison_values ) {
		std::cerr << value.first.first << " / " << value.first.second << ": P-value = " << value.second << ".\n" ;
	}

	std::cerr << " mismatches are:\n" ;
	typedef boost::bimap< std::string, int >::right_map::value_type V ;
	foreach( V const& value, call_names.right ) {
		std::cerr << value.first << ":" << value.second << " " ;
	}
	std::cerr << "\n" ;
	std::cerr << mismatches << "\n";
	std::cerr << "\n" ;
	std::cerr << "Consensus calls are: " ;
	std::copy( m_concordant_calls.begin(), m_concordant_calls.end(), std::ostream_iterator< std::string >( std::cerr, " " ) ) ;
	std::cerr << ".\n" ;
#endif
}
