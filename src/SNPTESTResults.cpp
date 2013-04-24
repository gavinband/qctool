
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "SNPTESTResults.hpp"

SNPTESTResults::SNPTESTResults(
	genfile::SNPIdentifyingDataTest::UniquePtr test
):
 	m_exclusion_test( test )
{}

void SNPTESTResults::add_variable( std::string const& variable ) {
	m_variables.insert( variable ) ;
	/* Make sure we prepare storage. */
	m_extra_variables[ variable ] ;
}

std::string SNPTESTResults::get_summary( std::string const& prefix, std::size_t target_column ) const {
	return prefix + "SNPTEST " + FlatFileFrequentistGenomeWideAssociationResults::get_summary( "", target_column ) ;
}

std::set< std::string > SNPTESTResults::get_desired_columns() const {
	std::set< std::string > desired_columns ;
	desired_columns.insert( "*_beta_1" ) ;
	desired_columns.insert( "*_beta_2" ) ;
	desired_columns.insert( "*_se_1" ) ;
	desired_columns.insert( "*_se_2" ) ;
	desired_columns.insert( "*_pvalue" ) ;
	desired_columns.insert( "info" ) ;
	desired_columns.insert( "all_maf" ) ;
	desired_columns.insert( "all_AA" ) ;
	desired_columns.insert( "all_AB" ) ;
	desired_columns.insert( "all_BB" ) ;
	desired_columns.insert( "all_NULL" ) ;
	desired_columns.insert( m_variables.begin(), m_variables.end() ) ;

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
	
	if( variable == "*_beta_1" ) {
		m_betas( snp_index, 0 ) = value ;
	}
	else if( variable == "*_se_1" ) {
		m_ses( snp_index, 0 ) = value ;
	}
	if( variable == "*_beta_2" ) {
		m_betas( snp_index, 1 ) = value ;
	}
	else if( variable == "*_se_2" ) {
		m_ses( snp_index, 1 ) = value ;
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
