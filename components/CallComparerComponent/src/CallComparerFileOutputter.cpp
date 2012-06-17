
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "components/CallComparerComponent/CallComparerFileOutputter.hpp"

CallComparerFileOutputter::UniquePtr CallComparerFileOutputter::create( std::string const& filename, std::string const& ) {
	return UniquePtr( new CallComparerFileOutputter( filename ) ) ;
}

CallComparerFileOutputter::SharedPtr CallComparerFileOutputter::create_shared( std::string const& filename, std::string const& ) {
	return SharedPtr( new CallComparerFileOutputter( filename ) ) ;
}

CallComparerFileOutputter::CallComparerFileOutputter( std::string const& filename ):
	m_filename( filename ),
	m_sink( statfile::BuiltInTypeStatSink::open( filename ))
{
	(*m_sink) | "SNPID" | "rsid" | "chromosome" | "position" | "alleleA" | "alleleB" | "callset_1" | "callset_2" | "comparison_method" | "comparison_variable" | "value" ;
}

void CallComparerFileOutputter::begin_comparisons( genfile::SNPIdentifyingData const& snp ) {
	m_snp = snp ;
}

void CallComparerFileOutputter::end_comparisons() {}

void CallComparerFileOutputter::set_result(
	std::string const& callset1,
	std::string const& callset2,
	std::string const& comparison_method,
	std::string const& comparison_variable,
	genfile::VariantEntry const& value
) {
	(*m_sink)
		<< m_snp.get_SNPID()
		<< m_snp.get_rsid()
		<< std::string( m_snp.get_position().chromosome() )
		<< m_snp.get_position().position()
		<< m_snp.get_first_allele()
		<< m_snp.get_second_allele()
		<< callset1
		<< callset2
		<< comparison_method
		<< comparison_variable ;
	(*m_sink)
		<< value.as< double >()
		<< statfile::end_row() ;
	;
}

void CallComparerFileOutputter::set_result(
	std::string const& comparison_method,
	std::string const& accepted_calls,
	PairwiseCallComparerManager::Calls const&
) {
	(*m_sink)
		<< m_snp.get_SNPID()
		<< m_snp.get_rsid()
		<< std::string( m_snp.get_position().chromosome() )
		<< m_snp.get_position().position()
		<< m_snp.get_first_allele()
		<< m_snp.get_second_allele()
		<< "NA"
		<< "NA"
		<< comparison_method
		<< "accepted_calls"
		<< accepted_calls
		<< statfile::end_row() ;
	;
}
