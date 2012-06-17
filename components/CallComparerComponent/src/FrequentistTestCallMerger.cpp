
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <map>
#include <string>
#include <boost/bimap.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <Eigen/Core>

#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "components/CallComparerComponent/PairwiseCallComparerManager.hpp"
#include "components/CallComparerComponent/FrequentistTestCallMerger.hpp"

//#define CONSENSUS_CALLER_DEBUG 1

FrequentistTestCallMerger::SharedPtr FrequentistTestCallMerger::create_shared( std::string const& comparison_method, double threshhold ) {
	return SharedPtr( new FrequentistTestCallMerger( comparison_method, threshhold )) ;
}

FrequentistTestCallMerger::FrequentistTestCallMerger( std::string const& comparison_method, double threshhold ):
	m_comparison_method( comparison_method ),
	m_pvalue_range( std::make_pair( threshhold, 100.0 ))
{}

std::string FrequentistTestCallMerger::get_spec() const {
	return m_spec ;
}

std::string FrequentistTestCallMerger::get_result_as_string() const {
	std::string result ;
	foreach( std::string const& call, m_concordant_calls ) {
		if( result.size() > 0 ) {
			result += "," ;
		}
		result += call ;
	}
	return result ;
}

void FrequentistTestCallMerger::begin_comparisons( genfile::SNPIdentifyingData const& snp ) {
	m_snp = snp ;
	m_comparison_values.clear() ;
	m_concordant_calls.clear() ;
}

void FrequentistTestCallMerger::set_result(
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

void FrequentistTestCallMerger::end_comparisons() {
	// Divide the calls into groups of calls all having p-value greater than the threshhold.
	// We assume comparisons are transitive: in reality this isn't the case but it's simpler to
	// assume transitivity here.

	typedef boost::bimap< std::string, int > CallNames ;
	CallNames call_names ;
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
		if( pvalue < m_pvalue_range.first || pvalue > m_pvalue_range.second ) {
			std::size_t const i = call_names.left.at( value.first.first ) ;
			std::size_t const j = call_names.left.at( value.first.second ) ;
			mismatches( i, j )++ ;
			mismatches( j, i )++ ;
		}
	}

#if CONSENSUS_CALLER_DEBUG
	std::cerr << std::resetiosflags( std::ios::floatfield ) ;
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
#endif

	for( CallNames::left_map::const_iterator i = call_names.left.begin(); i != call_names.left.end(); ++i ) {
		m_concordant_calls.insert( i->first ) ;
	}
	
	// At a really good SNP there are no mismatches (=> all rows are zero.)
	// At a quite good SNP we will have a (N-1)x(N-1) sub-table of zeroes.
	if( mismatches.array().maxCoeff() > 0 ) {
		// try knocking out one call.
		for( int i = 0; i < mismatches.rows(); ++i ) {
			Eigen::RowVectorXd saved_row = mismatches.row(i) ;
			Eigen::VectorXd saved_col = mismatches.col(i) ;
			mismatches.row(i).setZero() ;
			mismatches.col(i).setZero() ;
			if( mismatches.array().maxCoeff() == 0 ) {
				m_concordant_calls.erase( call_names.right.at( i ) ) ;
				mismatches.col(i) = saved_col ;
				mismatches.row(i) = saved_row ;
				break ;
			}
			mismatches.col(i) = saved_col ;
			mismatches.row(i) = saved_row ;
		}
		if( m_concordant_calls.size() == call_names.size() ) {
			m_concordant_calls.clear() ;
		}
	}
	

#if CONSENSUS_CALLER_DEBUG
	std::cerr << "Consensus calls are: " ;
	std::copy( m_concordant_calls.begin(), m_concordant_calls.end(), std::ostream_iterator< std::string >( std::cerr, " " ) ) ;
	std::cerr << ".\n" ;
#endif

	if( m_concordant_calls.size() == 1 ) {
		m_concordant_calls.clear() ;
	}
}
