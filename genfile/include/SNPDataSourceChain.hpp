#ifndef SNPDATAPROVIDERCHAIN_HPP
#define SNPDATAPROVIDERCHAIN_HPP

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
		SNPDataSourceChain(): m_current_provider(0), m_number_of_samples(0) {}

		~SNPDataSourceChain() {
			for( std::size_t i = 0; i < m_providers.size(); ++i ) {
				delete m_providers[i] ;
			}
		} ;

		void add_provider( std::auto_ptr< SNPDataSource > provider ) {
			if( m_providers.empty() ) {
				m_number_of_samples = provider->number_of_samples() ;
				m_current_provider = 0 ;
			}
			else if( provider->number_of_samples() != m_number_of_samples ) {
				throw FileContainsSNPsOfDifferentSizes() ;
			}
			m_providers.push_back( provider.release() ) ;
			m_numbers_of_snps_read.push_back( 0u ) ;
		}

		unsigned int number_of_samples() const { return m_number_of_samples ; }
		unsigned int total_number_of_snps() const {
			unsigned int total_number_of_snps = 0 ;
			for( std::size_t i = 0; i < m_providers.size(); ++i ) {
				total_number_of_snps += m_providers[i]->total_number_of_snps() ;
			}
			return total_number_of_snps ;
		}

		FormatType format() const {
			assert( m_current_provider < m_providers.size() ) ;
			return m_providers[ m_current_provider ]->format() ;
		}

		operator bool() const {
			if( m_current_provider < m_providers.size() ) {
				bool result = static_cast< bool >( *m_providers[ m_current_provider ] ) ;
				return result ;
			}
			else {
				return false ;
			}
		}

		std::istream const& stream() const {
			assert( m_current_provider < m_providers.size() ) ;
			return m_providers[ m_current_provider ]->stream() ;
		}

		std::istream& stream() {
			assert( m_current_provider < m_providers.size() ) ;
			return m_providers[ m_current_provider ]->stream() ;
		}


	private:

		void pre_read_snp() {
			// Make sure we switch to the next source when necessary.
			while(( m_current_provider < m_providers.size())
				&& (m_numbers_of_snps_read[ m_current_provider ] < m_providers[m_current_provider]->total_number_of_snps())) {
				++m_current_provider ;
			}
		}
	
		std::vector< SNPDataSource* > m_providers ;
		std::vector< std::size_t > m_numbers_of_snps_read ;
		std::size_t m_current_provider ;
		unsigned int m_number_of_samples ;
	} ;
}

#endif