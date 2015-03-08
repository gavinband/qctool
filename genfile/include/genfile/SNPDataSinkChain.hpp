
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef SNPDATAsinkCHAIN_HPP
#define SNPDATAsinkCHAIN_HPP

#include <iostream>
#include <string>
#include "genfile/snp_data_utils.hpp"
#include "genfile/SNPDataSink.hpp"

namespace genfile {
	// class SNPDataSinkChain represents a SNPDataSink
	// which outputs its data sequentially to a collection of other SNPDataSinks
	class SNPDataSinkChain: public SNPDataSink
	{
	public:
		SNPDataSinkChain() ;

		~SNPDataSinkChain() ;

		void add_sink( SNPDataSink::UniquePtr sink ) ;
		operator bool() const ;
		void move_to_next_sink() ;
		std::size_t number_of_sinks() const ;
		SNPDataSink const& sink( std::size_t i ) const ;
		std::size_t index_of_current_sink() const ;
		std::string get_spec() const ;
		SinkPos get_stream_pos() const ;

	private:
		void set_metadata_impl( Metadata const& metadata ) ;
		void set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter name_getter ) ;
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

		void write_variant_data_impl(
			SNPIdentifyingData const& id_data,
			VariantDataReader& data_reader,
			Info const& info
		) ;

		void finalise_impl() ;

	private:

		std::vector< SNPDataSink* > m_sinks ;
		std::size_t m_current_sink ;
	} ;
}

#endif
