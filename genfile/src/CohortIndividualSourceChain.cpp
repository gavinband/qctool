#include <boost/ptr_container/ptr_vector.hpp>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/CohortIndividualSourceChain.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/Error.hpp"

namespace genfile {
	CohortIndividualSourceChain::CohortIndividualSourceChain() {	
	}

	CohortIndividualSourceChain::~CohortIndividualSourceChain() {}
	
	void CohortIndividualSourceChain::add_source( CohortIndividualSource::UniquePtr source ) {
		if( m_sources.size() > 0 && source->get_column_spec() != m_sources[0].get_column_spec() ) {
			throw BadArgumentError( "CohortIndividualSourceChain::add_source()", "Column name mismatch in source." ) ;
		}
		m_sources.push_back( source.release() ) ;
	}
	
	std::size_t CohortIndividualSourceChain::get_number_of_individuals() const {
		std::size_t count = 0 ;
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			count += m_sources[i].get_number_of_individuals() ;
		}
		return count ;
	}

	CohortIndividualSource::ColumnSpec CohortIndividualSourceChain::get_column_spec() const {
		if( m_sources.size() == 0 ) {
			return ColumnSpec() ;
		}
		else {
			return m_sources[0].get_column_spec() ;
		}
	}

	CohortIndividualSource::Entry CohortIndividualSourceChain::get_entry( std::size_t sample_i, std::string const& column_name ) const {
		Entry result ;
		std::size_t i = 0 ;
		for( ; i < m_sources.size(); ++i ) {
			if( sample_i >= m_sources[i].get_number_of_individuals() ) {
				sample_i -= m_sources[i].get_number_of_individuals() ;
			}
			else {
				result = m_sources[i].get_entry( sample_i, column_name ) ;
				break ;
			}
		}
		assert( i < m_sources.size() ) ;
		return result ;
	}

	std::string CohortIndividualSourceChain::get_source_spec() const {
		return "(chain of " + string_utils::to_string( m_sources.size() ) + " sources)" ;
	}
}

