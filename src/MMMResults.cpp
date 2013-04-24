
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "MMMResults.hpp"

MMMResults::MMMResults(
	genfile::SNPIdentifyingDataTest::UniquePtr test
):
 	m_exclusion_test( test )
{}

void MMMResults::add_variable( std::string const& variable ) {
	m_variables.insert( variable ) ;
}

std::string MMMResults::get_summary( std::string const& prefix, std::size_t target_column ) const {
	return prefix + "mmm " + FlatFileFrequentistGenomeWideAssociationResults::get_summary( "", target_column ) ;
}

bool MMMResults::read_snp( statfile::BuiltInTypeStatSource& source, genfile::SNPIdentifyingData& snp ) const {
	return( source >> snp.position().chromosome() >> snp.SNPID() >> snp.rsid() >> snp.position().position() >> snp.first_allele() >> snp.second_allele() ) ;
}

bool MMMResults::check_if_snp_accepted( std::size_t snp_i ) const {
	return
		!m_exclusion_test.get()
		|| m_exclusion_test->operator()( m_snps[ snp_i ] )
	;
}

std::set< std::string > MMMResults::get_desired_columns() const {
	std::set< std::string > required_columns ;
	required_columns.insert( "z_est" ) ;
	required_columns.insert( "z_se" ) ;
	required_columns.insert( "pval" ) ;
	required_columns.insert( "var_info_all" ) ;
	required_columns.insert( "freq_1" ) ;
	required_columns.insert( "gen_00" ) ;
	required_columns.insert( "gen_01" ) ;
	required_columns.insert( "gen_11" ) ;
	required_columns.insert( "gen_NULL" ) ;
	required_columns.insert( m_variables.begin(), m_variables.end() ) ;
	return required_columns ;
}

void MMMResults::store_value(
	int snp_index,
	std::string const& variable,
	double const value
) {
	using genfile::string_utils::to_repr ;
	
	if( variable == "z_est" ) {
		m_betas( snp_index, 0 ) = value ;
	}
	else if( variable == "z_se" ) {
		m_ses( snp_index, 0 ) = value ;
	}
	else if( variable == "pval" ) {
		m_pvalues( snp_index ) = value ;
	}
	else if( variable == "var_info_all" ) {
		m_info( snp_index ) = value ;
	}
	else if( variable == "freq_1" ) {
		m_maf( snp_index ) = value ;
	}
	else if( variable == "gen_00" ) {
		m_sample_counts( snp_index, 0 ) = value ;
	}
	else if( variable == "gen_01" ) {
		m_sample_counts( snp_index, 1 ) = value ;
	}
	else if( variable == "gen_11" ) {
		m_sample_counts( snp_index, 2 ) = value ;
	}
	else if( variable == "gen_NULL" ) {
		m_sample_counts( snp_index, 3 ) = value ;
	}
}
