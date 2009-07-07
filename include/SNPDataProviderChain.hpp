#ifndef SNPDATAPROVIDERCHAIN_HPP
#define SNPDATAPROVIDERCHAIN_HPP

#include <iostream>
#include <string>
#include "snp_data_utils.hpp"
#include "SNPDataProvider.hpp"

// class SNPDataProviderChain represnets a SNPDataProvider
// which gets it data sequentially from a collection of other SNPDataProviders
class SNPDataProviderChain: public SNPDataProvider
{
public:
	SNPDataProviderChain(): m_current_provider(0), m_number_of_samples(0) {}

	~SNPDataProviderChain() {
		for( std::size_t i = 0; i < m_providers.size(); ++i ) {
			delete m_providers[i] ;
		}
	} ;

	void add_provider( std::auto_ptr< SNPDataProvider > provider ) {
		if( m_providers.empty() ) {
			m_number_of_samples = provider->number_of_samples() ;
			m_current_provider = 0 ;
		}
		else if( provider->number_of_samples() != m_number_of_samples ) {
			throw FileContainsSNPsOfDifferentSizes() ;
		}
		m_providers.push_back( provider.release() ) ;
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
			std::istream const& str = stream() ;
			return ( str ? true : false ) ;
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

	void prepare_to_read_snp() {
		// Make sure we switch to the next source when necessary.
		for (
			m_providers[ m_current_provider ]->stream().peek() ;
			(!(*m_providers[m_current_provider])) && (m_current_provider < m_providers.size()) ;
			++m_current_provider
		) {
			m_providers[ m_current_provider ]->stream().peek() ;
		}
	}
	
	std::vector< SNPDataProvider* > m_providers ;
	std::size_t m_current_provider ;
	unsigned int m_number_of_samples ;
} ;

#endif