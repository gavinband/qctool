
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef SNPDATAsinkCHAIN_HPP
#define SNPDATAsinkCHAIN_HPP

#include <iostream>
#include <string>
#include "snp_data_utils.hpp"
#include "SNPDataSink.hpp"

namespace genfile {
	// class SNPDataSinkChain represents a SNPDataSink
	// which outputs its data sequentially to a collection of other SNPDataSinks
	class SNPDataSinkChain: public SNPDataSink
	{
	public:
		SNPDataSinkChain(): m_current_sink(0) {}

		~SNPDataSinkChain() {
			for( std::size_t i = 0; i < m_sinks.size(); ++i ) {
				delete m_sinks[i] ;
			}
		} ;

		void add_sink( SNPDataSink::UniquePtr sink ) {
			if( m_sinks.empty() ) {
				m_current_sink = 0 ;
			}
			m_sinks.push_back( sink.release() ) ;
			
		}

		operator bool() const {
			if( m_current_sink < m_sinks.size() ) {
				bool result = static_cast< bool >( *m_sinks[ m_current_sink ] ) ;
				return result ;
			}
			else {
				return false ;
			}
		}

		void move_to_next_sink() {
			assert( m_current_sink < m_sinks.size() ) ;
			++m_current_sink ;
		}
		
		std::size_t number_of_sinks() const { return m_sinks.size() ; }

		SNPDataSink const& sink( std::size_t i ) const { return *m_sinks[ i ] ; }
		std::size_t index_of_current_sink() const { return m_current_sink ; }

		std::string get_spec() const {
			std::string result = "chain:" ;
			for( std::size_t i = 0; i < m_sinks.size(); ++i ) {
				if( i > 0 ) {
					result += "," ;
				}
				result += m_sinks[i]->get_spec() ;
			}
			return result ;
		}

	private:
		void set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter name_getter ) {
			for( std::size_t i = 0; i < m_sinks.size(); ++i ) {
				m_sinks[i]->set_sample_names( number_of_samples, name_getter ) ;
			}
		}
		
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

	private:

		std::vector< SNPDataSink* > m_sinks ;
		std::size_t m_current_sink ;
	} ;
}

#endif