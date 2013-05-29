
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <string>
#include "genfile/snp_data_utils.hpp"
#include "genfile/SNPDataSink.hpp"
#include "genfile/SNPDataSinkChain.hpp"

namespace genfile {
	SNPDataSinkChain::SNPDataSinkChain(): m_current_sink(0) {}

	SNPDataSinkChain::~SNPDataSinkChain() {
		for( std::size_t i = 0; i < m_sinks.size(); ++i ) {
			delete m_sinks[i] ;
		}
	}

	void SNPDataSinkChain::add_sink( SNPDataSink::UniquePtr sink ) {
		if( m_sinks.empty() ) {
			m_current_sink = 0 ;
		}
		m_sinks.push_back( sink.release() ) ;
		
	}

	SNPDataSinkChain::operator bool() const {
		if( m_current_sink < m_sinks.size() ) {
			bool result = static_cast< bool >( *m_sinks[ m_current_sink ] ) ;
			return result ;
		}
		else {
			return false ;
		}
	}

	void SNPDataSinkChain::move_to_next_sink() {
		assert( m_current_sink < m_sinks.size() ) ;
		++m_current_sink ;
	}
	
	std::size_t SNPDataSinkChain::number_of_sinks() const { return m_sinks.size() ; }

	SNPDataSink const& SNPDataSinkChain::sink( std::size_t i ) const { return *m_sinks[ i ] ; }
	std::size_t SNPDataSinkChain::index_of_current_sink() const { return m_current_sink ; }

	std::string SNPDataSinkChain::get_spec() const {
		std::string result = "chain:" ;
		for( std::size_t i = 0; i < m_sinks.size(); ++i ) {
			if( i > 0 ) {
				result += "," ;
			}
			result += m_sinks[i]->get_spec() ;
		}
		return result ;
	}

	SNPDataSinkChain::SinkPos SNPDataSinkChain::get_stream_pos() const {
		return sink( m_current_sink ).get_stream_pos() ;
	}

	void SNPDataSinkChain::set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter name_getter ) {
		for( std::size_t i = 0; i < m_sinks.size(); ++i ) {
			m_sinks[i]->set_sample_names( number_of_samples, name_getter ) ;
		}
	}
	
	void SNPDataSinkChain::write_snp_impl(
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
	) {
		assert( m_current_sink < m_sinks.size() ) ;
		m_sinks[m_current_sink]->write_snp(
			number_of_samples,
			SNPID,
			RSID,
			chromosome,
			SNP_position,
			first_allele,
			second_allele,
			get_AA_probability,
			get_AB_probability,
			get_BB_probability,
			info
		) ;
	}

	void SNPDataSinkChain::write_variant_data_impl(
		SNPIdentifyingData const& id_data,
		VariantDataReader& data_reader,
		Info const& info
	) {
		assert( m_current_sink < m_sinks.size() ) ;
		m_sinks[ m_current_sink ]->write_variant_data(
			id_data,
			data_reader,
			info
		) ;
	}
}
