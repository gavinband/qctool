
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <boost/regex.hpp>
#include <boost/format.hpp>
#include "genfile/Error.hpp"
#include "genfile/string_utils/string_utils.hpp"
#include "SNPTESTResults.hpp"
#include "EffectParameterNamePack.hpp"

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

EffectParameterNamePack SNPTESTResults::get_effect_parameter_names() const {
	return EffectParameterNamePack(
		m_beta_columns,
		m_se_columns,
		m_cov_columns
	) ;
}

namespace impl {
	boost::optional< SNPTESTResults::ColumnSpec > get_matching_name(
		std::vector< std::string > names,
		boost::regex const& regex,
		bool required,
		std::string const& type
	) {
		boost::optional< SNPTESTResults::ColumnSpec > result ;
		std::size_t i = 0 ;
		for( i = 0; i < names.size(); ++i ) {
			if( boost::regex_match( names[ i ], regex ) ) {
				break ;
			}
		}
		if( i < names.size() ) {
			result = SNPTESTResults::ColumnSpec( names[i], regex, i, required, type ) ;
		}
		return result ;
	}
	
	void insert_if_matched(
		std::vector< std::string > names,
		boost::regex const& regex,
		std::string const& type,
		std::vector< SNPTESTResults::ColumnSpec >* result
	) {
		assert( result ) ;
		boost::optional< SNPTESTResults::ColumnSpec > matched = get_matching_name( names, regex, false, type ) ;
		if( matched ) {
			result->push_back( matched.get() ) ;
		}
	}

	void insert_matched(
		std::vector< std::string > names,
		boost::regex const& regex,
		std::string const& type,
		std::vector< SNPTESTResults::ColumnSpec >* result
	) {
		assert( result ) ;
		boost::optional< SNPTESTResults::ColumnSpec > matched = get_matching_name( names, regex, true, type ) ;
		if( matched ) {
			result->push_back( matched.get() ) ;
		} else {
			throw genfile::BadArgumentError(
				"impl::insert_matched()",
				"names",
				( boost::format( "No name matched \"%s\"" ) % regex.str() ).str()
			) ;
		}
	}
}


SNPTESTResults::DesiredColumns SNPTESTResults::setup_columns( std::vector< std::string > const& column_names ) {
	using namespace genfile::string_utils ;
	using boost::regex ;
	DesiredColumns result ;
	
	{
		regex unwanted_column_details( "^.*frequentist_(add|dom|het|gen|rec)(_newml|ml|score|em|threshhold|expected)*_" ) ;

		std::vector< std::size_t > beta_indices ;
		for( std::size_t i = 0; i < 1000; ++i ) {
			regex const beta_regex = regex( replace_all( m_beta_column_regex, "<i>", to_string(i) ) ) ;
			regex const se_regex = regex( replace_all( replace_all( m_beta_column_regex, "<i>", to_string(i) ), "_beta", "_se" ) ) ;

			boost::optional< ColumnSpec > matched_beta_column = impl::get_matching_name( column_names, beta_regex, true, "beta" ) ;
			if( matched_beta_column ) {
				boost::optional< ColumnSpec > matched_se_column = impl::get_matching_name( column_names, se_regex, true, "se" ) ;
				if( !matched_se_column ) {
					throw genfile::BadArgumentError(
						"SNPTESTResults::setup_columns()",
						"column_names",
						( boost::format( "No column matching \"%s\" although one matches \"%s\"." ) % se_regex % beta_regex ).str()
					) ;
				}
				
				matched_beta_column->set_simplified_name(
					boost::regex_replace( matched_beta_column->name(), unwanted_column_details, "" )
				) ;
				matched_se_column->set_simplified_name(
					boost::regex_replace( matched_se_column->name(), unwanted_column_details, "" )
				) ;

				result.push_back( *matched_beta_column ) ;
				result.push_back( *matched_se_column ) ;
				m_beta_columns.push_back( matched_beta_column->simplified_name() ) ;
				m_se_columns.push_back( matched_se_column->simplified_name() ) ;
				beta_indices.push_back( i ) ;
			}
		}
	
		for( std::size_t i = 0; i < beta_indices.size(); ++i ) {
			for( std::size_t j = i+1; j < beta_indices.size(); ++j ) {
				std::string const cov_regex = replace_all(
					replace_all( m_beta_column_regex, "<i>", to_string( beta_indices[i] ) + "," + to_string( beta_indices[j] ) ),
					"_beta",
					"_cov"
				) ;
				boost::optional< ColumnSpec > matched_cov_column = impl::get_matching_name( column_names, regex( cov_regex ), true, "cov" ) ;
				if( !matched_cov_column ) {
					throw genfile::BadArgumentError(
						"SNPTESTResults::setup_columns()",
						"column_names",
						( boost::format( "No column matching \"%s\"." ) % cov_regex ).str()
					) ;
				}
				matched_cov_column->set_simplified_name(
					boost::regex_replace( matched_cov_column->name(), unwanted_column_details, "" )
				) ;
				result.push_back( *matched_cov_column ) ;
				m_cov_columns.push_back( matched_cov_column->simplified_name() ) ;
			}
		}
		
		if( m_beta_columns.size() == 0 ) {
			throw genfile::BadArgumentError(
				"SNPTESTResults::setup_columns()",
				"column_names",
				( boost::format( "No column matching beta regex (\"%s\") could be found." ) % m_beta_column_regex ).str()
			) ;
		}
	}

	{
		regex const pvalue_regex( replace_all( m_beta_column_regex, "_beta_<i>", "(ml|em|score|threshhold|lrt)_pvalue" ) ) ;
		boost::optional< ColumnSpec > matched_pvalue_column = impl::get_matching_name( column_names, pvalue_regex, true, "pvalue" ) ;
		if( !matched_pvalue_column ) {
			throw genfile::BadArgumentError(
				"SNPTESTResults::setup_columns()",
				"column_names",
				( boost::format( "No column matching overall pvalue regex (\"%s\") could be found." ) % pvalue_regex.str() ).str()
			) ;
		}
		result.push_back( matched_pvalue_column.get() ) ;
		m_pvalue_column = matched_pvalue_column.get().name() ;
	}

	{
		regex const info_regex( "(all_)?info" ) ;
		boost::optional< ColumnSpec > matched_info_column = impl::get_matching_name( column_names, info_regex, true, "info" ) ;
		if( !matched_info_column ) {
			throw genfile::BadArgumentError(
				"SNPTESTResults::setup_columns()",
				"column_names",
				( boost::format( "No column matching info regex (\"%s\") could be found." ) % info_regex.str() ).str()
			) ;
		}
		result.push_back( matched_info_column.get() ) ;
		m_info_column = matched_info_column.get().name() ;
	}
	impl::insert_matched( column_names, regex( "all_maf" ), "maf", &result ) ;
	impl::insert_matched( column_names, regex( "all_AA" ), "counts", &result ) ;
	impl::insert_matched( column_names, regex( "all_AB" ), "counts", &result ) ;
	impl::insert_matched( column_names, regex( "all_BB" ), "counts", &result ) ;
	impl::insert_if_matched( column_names, regex( "all_A" ), "counts", &result ) ;
	impl::insert_if_matched( column_names, regex( "all_B" ), "counts", &result ) ;
	impl::insert_matched( column_names, regex( "all_NULL" ), "counts", &result ) ;
	impl::insert_if_matched( column_names, regex( "cases_AA" ), "extra", &result ) ;
	impl::insert_if_matched( column_names, regex( "cases_AB" ), "extra", &result ) ;
	impl::insert_if_matched( column_names, regex( "cases_BB" ), "extra", &result ) ;
	impl::insert_if_matched( column_names, regex( "cases_A" ), "counts", &result ) ;
	impl::insert_if_matched( column_names, regex( "cases_B" ), "counts", &result ) ;
	impl::insert_if_matched( column_names, regex( "cases_NULL" ), "extra", &result ) ;
	impl::insert_if_matched( column_names, regex( "controls_AA" ), "extra", &result ) ;
	impl::insert_if_matched( column_names, regex( "controls_AB" ), "extra", &result ) ;
	impl::insert_if_matched( column_names, regex( "controls_BB" ), "extra", &result ) ;
	impl::insert_if_matched( column_names, regex( "controls_A" ), "extra", &result ) ;
	impl::insert_if_matched( column_names, regex( "controls_B" ), "extra", &result ) ;
	impl::insert_if_matched( column_names, regex( "controls_NULL" ), "extra", &result ) ;
	
	return result ;
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
	if( variable == m_pvalue_column ) {
		m_pvalues( snp_index ) = value ;
	}
	else if( variable == m_info_column ) {
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
