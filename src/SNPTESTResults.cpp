
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <boost/regex.hpp>
#include "genfile/string_utils/string_utils.hpp"
#include "SNPTESTResults.hpp"

SNPTESTResults::SNPTESTResults(
	genfile::SNPIdentifyingDataTest::UniquePtr test
):
 	m_exclusion_test( test ),
	m_beta_column_regex( ".*_beta_<i>($|[^0-9].*)" )
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
//	m_beta_column_regex = beta_column_regex ;
}

int SNPTESTResults::get_number_of_effect_parameters() const {
	return m_beta_columns.size() ;
}

namespace impl {
	bool matches_a_name( std::vector< std::string > names, boost::regex const& regex ) {
		std::size_t i = 0 ;
		for( i = 0; i < names.size(); ++i ) {
			if( boost::regex_match( names[ i ], regex ) ) {
				break ;
			}
		}
		return i < names.size() ;
	}
}


void SNPTESTResults::setup_columns( std::vector< std::string > const& column_names ) {
	using namespace genfile::string_utils ;
	
	DesiredColumns desired_columns ;
	
	for( std::size_t i = 0; i < 100; ++i ) {
		std::string const beta_regex = replace_all( m_beta_column_regex, "<i>", to_string(i) ) ;
		std::string const se_regex = replace_all( replace_all( m_beta_column_regex, "<i>", to_string(i) ), "_beta", "_se" ) ;

		std::size_t beta_column_i = 0 ;
		for( beta_column_i = 0; beta_column_i < column_names.size(); ++beta_column_i ) {
			boost::regex regex( beta_regex ) ;
			if( boost::regex_match( column_names[ beta_column_i ], regex ) ) {
				break ;
			}
		}
		if( impl::matches_a_name( column_names, boost::regex( beta_regex ) ) ) {
			desired_columns.insert( std::make_pair( beta_regex, true )) ;
			desired_columns.insert( std::make_pair( se_regex, true )) ;
			m_beta_columns.push_back( beta_regex ) ;
			m_se_columns.push_back( se_regex ) ;
			for( std::size_t j = i+1; j < 100; ++j ) {
				std::string const other_beta_regex = replace_all( m_beta_column_regex, "<i>", to_string(j) ) ;
				if( impl::matches_a_name( column_names, boost::regex( other_beta_regex ) ) ) {
					std::string const cov_regex = replace_all(
						replace_all( m_beta_column_regex, "<i>", to_string(i) + "," + to_string(j) ),
						"_beta",
						"_cov"
					) ;
					desired_columns.insert( std::make_pair( cov_regex, true )) ;
					m_cov_columns.push_back( cov_regex ) ;
				}
			}
		}
	}
	
	desired_columns.insert( std::make_pair( replace_all( m_beta_column_regex, "_beta_<i>", "(ml|em|score|threshhold|lrt)_pvalue" ), true ) ) ;
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
			std::set< std::pair< std::string, bool > >::iterator where = desired_columns.find( std::make_pair( *i, false ) ) ;
			if( where != desired_columns.end() ) {
				desired_columns.erase( where ) ;
			}
			desired_columns.insert( std::make_pair( *i, false ) ) ;
		}
	}
	m_desired_columns = desired_columns ;
}

std::set< std::pair< std::string, bool > > SNPTESTResults::get_desired_columns() const {
	return m_desired_columns ;
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
	bool matched = false ;
	
	for( std::size_t i = 0; !matched && i < m_beta_columns.size(); ++i ) {
		if( variable == m_beta_columns[i] ) {
			m_betas( snp_index, i ) = value ;
			matched = true ;
		} else if( variable == m_se_columns[i] ) {
			m_ses( snp_index, i ) = value ;
			matched = true ;
		}
	}
	for( std::size_t i = 0; !matched && i < m_cov_columns.size(); ++i ) {
		if( variable == m_cov_columns[i] ) {
			m_covariance( snp_index, i ) = value ;
			matched = true ;
		}
	}

	if( matched ) {
		return ;
	}
	if( variable == ".*_pvalue" ) {
		m_pvalues( snp_index ) = value ;
	}
	else if( variable == "(all_)?info" ) {
		m_info( snp_index ) = value ;
	}
	else if( variable == "all_maf" ) {
		m_maf( snp_index ) = value ;
	}
	else if( variable == "all_A" ) {
		m_sample_counts( snp_index, 0 ) = value ;
	}
	else if( variable == "all_B" ) {
		m_sample_counts( snp_index, 1 ) = value ;
	}
	else if( variable == "all_AA" ) {
		m_sample_counts( snp_index, 2 ) = value ;
	}
	else if( variable == "all_AB" ) {
		m_sample_counts( snp_index, 3 ) = value ;
	}
	else if( variable == "all_BB" ) {
		m_sample_counts( snp_index, 4 ) = value ;
	}
	else if( variable == "all_NULL" ) {
		m_sample_counts( snp_index, 5 ) = value ;
	}
	else if( m_variables.find( variable ) != m_variables.end() ) {
		m_extra_variables[ variable ][ snp_index ] = value ;
	}
}
