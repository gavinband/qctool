#ifndef SNPDATAsourceCHAIN_HPP
#define SNPDATAsourceCHAIN_HPP

#include <iostream>
#include <string>
#include "snp_data_utils.hpp"
#include "SNPDataSource.hpp"

namespace genfile {
	// class SNPDataSourceChain represnets a SNPDataSource
	// which gets it data sequentially from a collection of other SNPDataSources
	class SNPDataSourceChain: public SNPDataSource
	{
	public:
		SNPDataSourceChain(): m_current_source(0), m_number_of_samples(0) {}

		~SNPDataSourceChain() {
			for( std::size_t i = 0; i < m_sources.size(); ++i ) {
				delete m_sources[i] ;
			}
		} ;

		void add_source( std::auto_ptr< SNPDataSource > source ) {
			if( m_sources.empty() ) {
				m_number_of_samples = source->number_of_samples() ;
				m_current_source = 0 ;
			}
			else if( source->number_of_samples() != m_number_of_samples ) {
				throw FileContainsSNPsOfDifferentSizes() ;
			}
			m_sources.push_back( source.release() ) ;
			m_number_of_snps_read.push_back( 0u ) ;
		}

		unsigned int number_of_samples() const { return m_number_of_samples ; }
		unsigned int total_number_of_snps() const {
			unsigned int total_number_of_snps = 0 ;
			for( std::size_t i = 0; i < m_sources.size(); ++i ) {
				total_number_of_snps += m_sources[i]->total_number_of_snps() ;
			}
			return total_number_of_snps ;
		}

		FormatType format() const {
			assert( m_current_source < m_sources.size() ) ;
			return m_sources[ m_current_source ]->format() ;
		}

		operator bool() const {
			if( m_current_source < m_sources.size() ) {
				bool result = static_cast< bool >( *m_sources[ m_current_source ] ) ;
				return result ;
			}
			else {
				return false ;
			}
		}

		std::istream const& stream() const {
			assert( m_current_source < m_sources.size() ) ;
			return m_sources[ m_current_source ]->stream() ;
		}

		std::istream& stream() {
			assert( m_current_source < m_sources.size() ) ;
			return m_sources[ m_current_source ]->stream() ;
		}


	private:

		void pre_read_snp() {
			// Make sure we switch to the next source when necessary.
			while(( m_current_source < m_sources.size())
				&& (m_number_of_snps_read[ m_current_source ] >= m_sources[m_current_source]->total_number_of_snps())) {
				++m_current_source ;
			}
		}

		void post_read_snp() {
			++m_number_of_snps_read[ m_current_source ] ;
		}
	
		std::vector< SNPDataSource* > m_sources ;
		std::vector< std::size_t > m_number_of_snps_read ;
		std::size_t m_current_source ;
		unsigned int m_number_of_samples ;
	} ;
}

#endif