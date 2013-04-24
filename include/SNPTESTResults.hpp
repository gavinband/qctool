
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef SNPTEST_RESULTS_HPP
#define SNPTEST_RESULTS_HPP

#include "FlatFileFrequentistGenomeWideAssociationResults.hpp"

struct SNPTESTResults: public FlatFileFrequentistGenomeWideAssociationResults {
	SNPTESTResults(
		genfile::SNPIdentifyingDataTest::UniquePtr test
	) ;

	void add_variable( std::string const& variable ) ;
	std::string get_summary( std::string const& prefix, std::size_t target_column ) const ;

private:
	genfile::SNPIdentifyingDataTest::UniquePtr m_exclusion_test ;
	std::set< std::string > m_variables ;
	
	std::set< std::string > get_desired_columns() const ;
	bool read_snp( statfile::BuiltInTypeStatSource& source, genfile::SNPIdentifyingData& snp ) const ;
	bool check_if_snp_accepted( std::size_t snp_i ) const ;
	void store_value(
		int snp_index,
		std::string const& variable,
		double value
	) ;
} ;



#endif

