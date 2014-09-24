
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

	void set_effect_size_column_regex( std::string const& beta_column_regex ) ;
	EffectParameterNamePack get_effect_parameter_names() const ;
	void add_variable( std::string const& variable ) ;
	
	std::string get_summary( std::string const& prefix, std::size_t target_column ) const ;

private:
	genfile::SNPIdentifyingDataTest::UniquePtr m_exclusion_test ;
	std::string m_beta_column_regex ;
	std::vector< std::string > m_beta_columns ;
	std::vector< std::string > m_se_columns ;
	std::vector< std::string > m_cov_columns ;
	std::string m_pvalue_column ;
	std::string m_info_column ;
	std::set< std::string > m_variables ;

private:
	DesiredColumns setup_columns( std::vector< std::string > const& column_names ) ;
	bool read_snp( statfile::BuiltInTypeStatSource& source, genfile::SNPIdentifyingData& snp ) const ;
	bool check_if_snp_accepted( std::size_t snp_i ) const ;
	void store_value(
		int snp_index,
		std::string const& variable,
		double value
	) ;
} ;



#endif

