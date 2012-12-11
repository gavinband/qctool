
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

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
		void write_header( std::size_t number_of_samples, SampleNameGetter sample_name_getter ) const ;
		
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
			GenotypeProbabilityGetter const& get_BB_probability,
			Info const& info
		) ;
		
		std::string get_spec() const ;

		void set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter ) ;
		std::ostream::streampos get_stream_pos() const ;
		
	private:
		std::string const m_filename ;
		CompressionType const m_compression_type ;
		std::auto_ptr< std::ostream > m_stream_ptr ;
		bool m_have_written_header ;
		std::size_t m_number_of_samples ;
		double const m_call_threshhold ;
		std::vector< std::string > m_output_fields ;
	private:
		void write_info( Info const& info ) ;
	} ;
}

#endif
