
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef SNP_ID_SINK_HPP
#define SNP_ID_SINK_HPP

#include <string>
#include <iostream>
#include <memory>
#include <stddef.h>

#include "SNPInListCondition.hpp"
#include "genfile/SNPDataSink.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"

struct SNPIDSink: public genfile::SNPDataSink
{
	SNPIDSink( std::string const& filename ) {
		setup( filename ) ;
	}
	
	operator bool() const { return *m_sink ; }
	
	void write_snp_impl(
		uint32_t,
		std::string SNPID,
		std::string RSID,
		genfile::Chromosome chromosome,
		uint32_t SNP_position,
		std::string first_allele,
		std::string second_allele,
		GenotypeProbabilityGetter const&,
		GenotypeProbabilityGetter const&,
		GenotypeProbabilityGetter const&,
		Info const& info
	) {
		(*m_sink)
			<< SNPID
			<< RSID
			<< chromosome
			<< SNP_position
			<< first_allele
			<< second_allele
			<< statfile::end_row() ;
	} ;
	
	std::string get_spec() const { return "(SNPIDSink)" ; }
	
private:
	
	void setup( std::string const& filename ) {
		m_sink = statfile::BuiltInTypeStatSink::open( filename ) ;
		(*m_sink ) | "SNPID" | "rsid" | "chromosome" | "position" | "alleleA" | "alleleB" ;
	}
	
	statfile::BuiltInTypeStatSink::UniquePtr m_sink ;
} ;

#endif
