
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <map>
#include "genfile/CohortIndividualSource.hpp"
#include "components/SNPSummaryComponent/StratifyingSNPSummaryComputation.hpp"
                                                       
StratifyingSNPSummaryComputation::StratifyingSNPSummaryComputation( SNPSummaryComputation::UniquePtr computation, std::string const& stratification_name, StrataMembers const& strata_members ):
 	m_computation( computation ),
	m_stratification_name( stratification_name ),
	m_strata_members( strata_members )
{}

namespace impl {
	struct NameMungingSNPComputationCallback {
		NameMungingSNPComputationCallback( SNPSummaryComputation::ResultCallback callback, std::string const& prefix, std::string const& suffix ):
			m_callback( callback ),
			m_prefix( prefix ),
			m_suffix( suffix )
		{
			assert( m_callback ) ;
		}

		NameMungingSNPComputationCallback( NameMungingSNPComputationCallback const& other ):
			m_callback( other.m_callback ),
			m_prefix( other.m_prefix ),
			m_suffix( other.m_suffix )
		{}

		NameMungingSNPComputationCallback& operator=( NameMungingSNPComputationCallback const& other ) {
			m_callback = other.m_callback ;
			m_prefix = other.m_prefix ;
			m_suffix = other.m_suffix ;
			return *this ;
		}
		
		void operator()( std::string const& value_name, genfile::VariantEntry const& value ) const {
			m_callback( m_prefix + value_name + m_suffix, value ) ;
		}
		
	private:
		SNPSummaryComputation::ResultCallback m_callback ;
		std::string m_prefix ;
		std::string m_suffix ;
	} ;
}

void StratifyingSNPSummaryComputation::operator()(
	SNPIdentifyingData const& snp,
	Genotypes const& genotypes,
	SampleSexes const& sexes,
	genfile::VariantDataReader& data_reader,
	ResultCallback callback
) {
	for( StrataMembers::const_iterator strata_i = m_strata_members.begin(); strata_i != m_strata_members.end(); ++strata_i ) {
		std::vector< int > const& members = strata_i->second ;
		Genotypes strata_genotypes = Genotypes::Constant( members.size(), genotypes.cols(), 0 ) ;
		SampleSexes strata_sexes( members.size(), '.' ) ;
		for( std::size_t i = 0; i < members.size(); ++i ) {
			strata_genotypes.row( i ) = genotypes.row( members[i] ) ;
			strata_sexes[ i ] = sexes[ members[i] ] ;
		}
		
		m_computation->operator()(
			snp,
			strata_genotypes,
			strata_sexes,
			data_reader,
			impl::NameMungingSNPComputationCallback( callback, "", " [" + m_stratification_name + "=" + genfile::string_utils::to_string( strata_i->first ) + "] " )
		) ;
	}
}

std::string StratifyingSNPSummaryComputation::get_summary( std::string const& prefix, std::size_t column_width ) const {
	return m_computation->get_summary( prefix, column_width ) + ", stratified on \"" + m_stratification_name + "\"." ;
}

