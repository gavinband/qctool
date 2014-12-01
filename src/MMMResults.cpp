
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <boost/regex.hpp>
#include <boost/optional.hpp>
#include "genfile/Chromosome.hpp"
#include "MMMResults.hpp"
#include "EffectParameterNamePack.hpp"

// #define DEBUG_MMMRESULTS 1

namespace {
	namespace impl {	
		boost::optional< MMMResults::ColumnSpec > get_matching_name(
			std::vector< std::string > names,
			boost::regex const& regex,
			bool required,
			std::string const& type
		) {
			boost::optional< MMMResults::ColumnSpec > result ;
			std::size_t i = 0 ;
			for( i = 0; i < names.size(); ++i ) {
#if DEBUG_MMMRESULTS
				std::cerr << "Looking for regex " << regex.str() << "in " << names[i] << "...\n" ;
#endif
				if( boost::regex_match( names[ i ], regex ) ) {
#if DEBUG_MMMRESULTS
					std::cerr << "Found!\n" ;
#endif
					break ;
				}
			}
			if( i < names.size() ) {
				result = MMMResults::ColumnSpec( names[i], regex, i, required, type ) ;
			}
			return result ;
		}
	}
}
MMMResults::MMMResults(
	genfile::SNPIdentifyingDataTest::UniquePtr test,
	boost::optional< genfile::Chromosome > const chromosome_hint
):
 	m_exclusion_test( test ),
	m_effect_column_regex( "z_est" ),
	m_se_column_regex( "z_se" ),
	m_chromosome_hint( chromosome_hint )
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

int MMMResults::get_number_of_effect_parameters() const {
	return 1 ;
}

EffectParameterNamePack MMMResults::get_effect_parameter_names() const {
	return EffectParameterNamePack(
		std::vector< std::string >( 1, m_effect_column_regex ),
		std::vector< std::string >( 1, m_se_column_regex ),
		std::vector< std::string >()
	) ;
}


bool MMMResults::read_snp( statfile::BuiltInTypeStatSource& source, genfile::SNPIdentifyingData& snp ) const {
	if( m_chromosome_hint ) {
		snp.position().chromosome() = *m_chromosome_hint ;
	}
	return( source >> snp.SNPID() >> snp.rsid() >> snp.position().position() >> snp.first_allele() >> snp.second_allele() ) ;
}

bool MMMResults::check_if_snp_accepted( std::size_t snp_i ) const {
	return
		!m_exclusion_test.get()
		|| m_exclusion_test->operator()( m_snps[ snp_i ] )
	;
}


MMMResults::DesiredColumns MMMResults::setup_columns( std::vector< std::string > const& columns ) {
	DesiredColumns result ;
	using boost::regex ;
	result.push_back( *impl::get_matching_name( columns, regex( m_effect_column_regex ), true, "beta" ) ) ;
	result.push_back( *impl::get_matching_name( columns, regex( m_se_column_regex ), true, "se" ) ) ;
	result.push_back( *impl::get_matching_name( columns, regex( "pval" ), true, "pvalue" ) ) ;
	result.push_back( *impl::get_matching_name( columns, regex( "var_info_all" ), true, "info" ) ) ;
	result.push_back( *impl::get_matching_name( columns, regex( "freq_1" ), true, "info" ) ) ;
	result.push_back( *impl::get_matching_name( columns, regex( "gen_00" ), true, "counts" ) ) ;
	result.push_back( *impl::get_matching_name( columns, regex( "gen_01" ), true, "counts" ) ) ;
	result.push_back( *impl::get_matching_name( columns, regex( "gen_11" ), true, "counts" ) ) ;
	result.push_back( *impl::get_matching_name( columns, regex( "gen_NULL" ), true, "counts" ) ) ;
	
	return result ;	
}

void MMMResults::store_value(
	int snp_index,
	std::string const& variable,
	std::string const& value
) {
	using genfile::string_utils::to_repr ;
    double const NA = std::numeric_limits< double >::quiet_NaN() ;
	if( variable == m_effect_column_regex ) {
		m_betas( snp_index, 0 ) = ( value == "NA" ? NA: to_repr< double >( value ) ) ;
	}
	else if( variable == m_se_column_regex ) {
		m_ses( snp_index, 0 ) = ( value == "NA" ? NA: to_repr< double >( value ) ) ;
	}
	else if( variable == "pval" ) {
		m_pvalues( snp_index ) = ( value == "NA" ? NA: to_repr< double >( value ) ) ;
	}
	else if( variable == "var_info_all" ) {
		m_info( snp_index ) = ( value == "NA" ? NA: to_repr< double >( value ) ) ;
	}
	else if( variable == "freq_1" ) {
		m_maf( snp_index ) = ( value == "NA" ? NA: to_repr< double >( value ) ) ;
	}
	else if( variable == "gen_00" ) {
		m_sample_counts( snp_index, 2 ) = ( value == "NA" ? NA: to_repr< double >( value ) ) ;
	}
	else if( variable == "gen_01" ) {
		m_sample_counts( snp_index, 3 ) = ( value == "NA" ? NA: to_repr< double >( value ) ) ;
	}
	else if( variable == "gen_11" ) {
		m_sample_counts( snp_index, 4 ) = ( value == "NA" ? NA: to_repr< double >( value ) ) ;
	}
	else if( variable == "gen_NULL" ) {
		m_sample_counts( snp_index, 5 ) = ( value == "NA" ? NA: to_repr< double >( value ) ) ;
	}
}
