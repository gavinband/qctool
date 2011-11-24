#ifndef GENFILE_VCF_FORMAT_SNP_DATA_SINK_HPP
#define GENFILE_VCF_FORMAT_SNP_DATA_SINK_HPP

#include <string>
#include <iostream>
#include <boost/function.hpp>
#include "genfile/SNPDataSink.hpp"

namespace genfile {
	struct VCFFormatSNPDataSink: public SNPDataSink {
		VCFFormatSNPDataSink( std::string const& filename ) ;

	private:
		void write_header( std::size_t number_of_samples ) const ;
		
		// Methods required by SNPDataSink
		operator bool() const { return *m_stream_ptr ; }
		
		void write_snp_impl(
			uint32_t number_of_samples,
			std::string SNPID,
			std::string RSID,
			Chromosome chromosome,
			uint32_t SNP_position,
			std::string first_allele,
			std::string second_allele,
			GenotypeProbabilityGetter const& get_AA_probability,
			GenotypeProbabilityGetter const& get_AB_probability,
			GenotypeProbabilityGetter const& get_BB_probability
		) ;
		
	private:
		std::string m_filename ;
		std::auto_ptr< std::ostream > m_stream_ptr ;
		bool m_have_written_header ;
		std::size_t m_number_of_samples ;
		double const m_call_threshhold ;
		boost::function< std::string ( std::size_t ) > m_sample_name_getter ;
	} ;
}

#endif
