
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "MMMResults.hpp"

MMMResults::MMMResults(
	genfile::SNPIdentifyingDataTest::UniquePtr test
):
 	m_exclusion_test( test ),
	m_effect_column_regex( "z_est" ),
	m_se_column_regex( "z_se" )
{}

void MMMResults::add_variable( std::string const& variable ) {
	m_variables.insert( variable ) ;
}

void MMMResults::set_effect_size_column_regex( std::string const& effect_column_regex ) {
	m_effect_column_regex = effect_column_regex ;
	m_se_column_regex = genfile::string_utils::replace_all( effect_column_regex, "_est", "_se" ) ;
}

std::string MMMResults::get_summary( std::string const& prefix, std::size_t target_column ) const {
	return prefix + "mmm " + FlatFileFrequentistGenomeWideAssociationResults::get_summary( "", target_column ) ;
}

int const MMMResults::get_number_of_effect_parameters() const {
	return 1 ;
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

std::set< std::pair< std::string, bool > > MMMResults::get_desired_columns() const {
	std::set< std::pair< std::string, bool > > required_columns ;
	required_columns.insert( std::make_pair( m_effect_column_regex, true ) ) ;
	required_columns.insert( std::make_pair( m_se_column_regex, true ) ) ;
	required_columns.insert( std::make_pair( "pval", true ) ) ;
	required_columns.insert( std::make_pair( "var_info_all", true ) ) ;
	required_columns.insert( std::make_pair( "freq_1", true ) ) ;
	required_columns.insert( std::make_pair( "gen_00", true ) ) ;
	required_columns.insert( std::make_pair( "gen_01", true ) ) ;
	required_columns.insert( std::make_pair( "gen_11", true ) ) ;
	required_columns.insert( std::make_pair( "gen_NULL", true ) ) ;
	{
		std::set< std::string >::const_iterator
			i = m_variables.begin(),
			end_i = m_variables.end() ;
		for( ; i != end_i; ++i ) {
			required_columns.insert( std::make_pair( *i, true ) ) ;
		}
	}
	return required_columns ;
}

void MMMResults::store_value(
	int snp_index,
	std::string const& variable,
	double const value
) {
	using genfile::string_utils::to_repr ;
	
	if( variable == m_effect_column_regex ) {
		m_betas( snp_index, 0 ) = value ;
	}
	else if( variable == m_se_column_regex ) {
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
