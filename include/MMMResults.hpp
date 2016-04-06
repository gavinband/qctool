
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef MMM_RESULTS_HPP
#define MMM_RESULTS_HPP

#include <boost/optional.hpp>
#include "genfile/Chromosome.hpp"
#include "genfile/VariantIdentifyingDataTest.hpp"
#include "FlatFileFrequentistGenomeWideAssociationResults.hpp"

struct MMMResults: public FlatFileFrequentistGenomeWideAssociationResults {
	MMMResults(
		genfile::VariantIdentifyingDataTest::UniquePtr test,
		boost::optional< genfile::Chromosome > const chromosome_hint = boost::optional< genfile::Chromosome >()
	) ;

	void set_effect_size_column_regex( std::string const& beta_column_regex ) ;
	void add_variable( std::string const& variable ) ;
	std::string get_summary( std::string const& prefix, std::size_t target_column ) const ;
	int get_number_of_effect_parameters() const ;
	EffectParameterNamePack get_effect_parameter_names() const ;

private:
	genfile::VariantIdentifyingDataTest::UniquePtr m_exclusion_test ;
	std::string m_effect_column_regex ;
	std::string m_se_column_regex ;
	boost::optional< genfile::Chromosome > m_chromosome_hint ;
	std::set< std::string > m_variables ;
	DesiredColumns m_required_columns ;
	
private:
	DesiredColumns setup_columns( std::vector< std::string > const& column_names ) ;
	bool read_snp( statfile::BuiltInTypeStatSource& source, genfile::SNPIdentifyingData& snp ) const ;
	bool check_if_snp_accepted( std::size_t snp_i ) const ;
	void store_value(
		int snp_index,
		std::string const& variable,
		std::string const& value
	) ;
} ;



#endif
