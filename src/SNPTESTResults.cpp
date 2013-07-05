
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "SNPTESTResults.hpp"

SNPTESTResults::SNPTESTResults(
	genfile::SNPIdentifyingDataTest::UniquePtr test
):
 	m_exclusion_test( test ),
	m_beta_column_regex( ".*beta_1.*" ),
	m_se_column_regex( ".*se_1.*" )
{}

void SNPTESTResults::add_variable( std::string const& variable ) {
	m_variables.insert( variable ) ;
	/* Make sure we prepare storage. */
	m_extra_variables[ variable ] ;
}

std::string SNPTESTResults::get_summary( std::string const& prefix, std::size_t target_column ) const {
	return prefix + "SNPTEST " + FlatFileFrequentistGenomeWideAssociationResults::get_summary( "", target_column ) ;
}

void SNPTESTResults::set_effect_size_column_regex( std::string const& beta_column_regex ) {
	m_beta_column_regex = beta_column_regex ;
	m_se_column_regex = genfile::string_utils::replace_all( beta_column_regex, "beta", "se" ) ;
}

std::set< std::pair< std::string, bool > > SNPTESTResults::get_desired_columns() const {
	std::set< std::pair< std::string, bool > > desired_columns ;
	desired_columns.insert( std::make_pair( m_beta_column_regex, true ) ) ;
	desired_columns.insert( std::make_pair( m_se_column_regex, true ) ) ;
	desired_columns.insert( std::make_pair( ".*_pvalue", true ) ) ;
	desired_columns.insert( std::make_pair( "(all_)?info", true ) ) ;
	desired_columns.insert( std::make_pair( "all_maf", true ) ) ;
	desired_columns.insert( std::make_pair( "all_AA", true ) ) ;
	desired_columns.insert( std::make_pair( "all_AB", true ) ) ;
	desired_columns.insert( std::make_pair( "all_BB", true ) ) ;
	desired_columns.insert( std::make_pair( "all_A", false ) ) ;
	desired_columns.insert( std::make_pair( "all_B", false ) ) ;
	desired_columns.insert( std::make_pair( "all_NULL", true ) ) ;
	{
		std::set< std::string >::const_iterator
			i = m_variables.begin(),
			end_i = m_variables.end() ;
		for( ; i != end_i; ++i ) {
			desired_columns.insert( std::make_pair( *i, true ) ) ;
		}
	}
	return desired_columns ;
}

bool SNPTESTResults::read_snp( statfile::BuiltInTypeStatSource& source, genfile::SNPIdentifyingData& snp ) const {
	return( source >> snp.SNPID() >> snp.rsid() >> snp.position().chromosome() >> snp.position().position() >> snp.first_allele() >> snp.second_allele() ) ;
}

bool SNPTESTResults::check_if_snp_accepted( std::size_t snp_i ) const {
	return
		! m_exclusion_test.get()
		|| m_exclusion_test->operator()( m_snps[ snp_i ] )
	;
}

void SNPTESTResults::store_value(
	int snp_index,
	std::string const& variable,
	double value
) {
	using genfile::string_utils::to_repr ;
	
	if( variable == m_beta_column_regex ) {
		m_betas( snp_index, 0 ) = value ;
	}
	else if( variable == m_se_column_regex ) {
		m_ses( snp_index, 0 ) = value ;
	}
	else if( variable == "*_pvalue" ) {
		m_pvalues( snp_index ) = value ;
	}
	else if( variable == "info" ) {
		m_info( snp_index ) = value ;
	}
	else if( variable == "all_maf" ) {
		m_maf( snp_index ) = value ;
	}
	else if( variable == "all_AA" ) {
		m_sample_counts( snp_index, 0 ) = value ;
	}
	else if( variable == "all_AB" ) {
		m_sample_counts( snp_index, 1 ) = value ;
	}
	else if( variable == "all_BB" ) {
		m_sample_counts( snp_index, 2 ) = value ;
	}
	else if( variable == "all_NULL" ) {
		m_sample_counts( snp_index, 3 ) = value ;
	}
	else if( m_variables.find( variable ) != m_variables.end() ) {
		m_extra_variables[ variable ][ snp_index ] = value ;
	}
}
