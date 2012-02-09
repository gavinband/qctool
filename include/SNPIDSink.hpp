#ifndef SNP_ID_SINK_HPP
#define SNP_ID_SINK_HPP

#include <string>
#include <iostream>
#include <memory>
#include <stddef.h>

#include "SNPInListCondition.hpp"
#include "genfile/SNPDataSink.hpp"

struct SNPIDSink: public genfile::SNPDataSink
{
	SNPIDSink( std::string const& filename ) {
		setup( filename ) ;
	}
	
	operator bool() const { return *m_stream_ptr ; }
	
	void write_snp_impl(
		uint32_t,
		std::string SNPID,
		std::string RSID,
		genfile::Chromosome chromosome,
		uint32_t SNP_position,
		std::string,
		std::string,
		GenotypeProbabilityGetter const&,
		GenotypeProbabilityGetter const&,
		GenotypeProbabilityGetter const&
	) {
		stream() << SNPID << " " << RSID << " " << chromosome << " " << SNP_position << "\n" ;
	} ;
	
	std::string get_spec() const { return "(SNPIDSink)" ; }
	
private:
	
	std::ostream& stream() { return *m_stream_ptr ; }
	
	void setup( std::string const& filename ) {
		m_stream_ptr = open_file_for_output( filename ) ;
	}
	
	std::auto_ptr< std::ostream > m_stream_ptr ;
} ;

#endif
