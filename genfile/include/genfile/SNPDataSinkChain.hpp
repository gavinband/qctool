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

		void add_sink( std::auto_ptr< SNPDataSink > sink ) {
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

		SNPDataSink const& sink( std::size_t i ) const { return *m_sinks[ i ] ; }
		std::size_t index_of_current_sink() const { return m_current_sink ; }

	private:
		void write_snp_impl(
			uint32_t number_of_samples,
			std::string SNPID,
			std::string RSID,
			Chromosome chromosome,
			uint32_t SNP_position,
			char first_allele,
			char second_allele,
			GenotypeProbabilityGetter const& get_AA_probability,
			GenotypeProbabilityGetter const& get_AB_probability,
			GenotypeProbabilityGetter const& get_BB_probability
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
				get_BB_probability
			) ;
		}

	private:

		std::vector< SNPDataSink* > m_sinks ;
		std::size_t m_current_sink ;
	} ;
}

#endif