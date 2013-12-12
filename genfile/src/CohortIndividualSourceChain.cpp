
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <boost/ptr_container/ptr_vector.hpp>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/CohortIndividualSourceChain.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/Error.hpp"

namespace genfile {
	CohortIndividualSourceChain::CohortIndividualSourceChain() {	
	}

	CohortIndividualSourceChain::~CohortIndividualSourceChain() {}
	
	void CohortIndividualSourceChain::add_source( CohortIndividualSource::UniquePtr source, std::string const& source_name ) {
		m_sources.push_back( source.release() ) ;
        m_source_names.push_back( source_name ) ;
		if( m_sources.size() == 1 ) {
			m_column_spec = m_sources[0].get_column_spec() ;
		}
		else {
			ColumnSpec new_column_spec = m_sources.back().get_column_spec() ;
            if( !new_column_spec.find( "qctool:cohort_name" )) {
                new_column_spec.add_column( "qctool:cohort_name", e_DISCRETE_COVARIATE ) ;
            }
			std::vector< std::string > column_names = new_column_spec.get_names() ;
			for( std::size_t i = 0; i < column_names.size(); ++i ) {
				// we remove any columns which match name but mismatch type, and any columns which match completely.
				// we extend our set of columns with any genuinely new columns.
				if( m_column_spec.check_for_column( column_names[i] ) ) {
					if( m_column_spec[ column_names[i] ] != new_column_spec[ column_names[i] ] ) {
						m_column_spec.remove( column_names[i] ) ;
						new_column_spec.remove( column_names[i] ) ;
					}
					else {
						new_column_spec.remove( column_names[i] ) ;
					}
				}
			}
			m_column_spec = m_column_spec + new_column_spec ;
		}
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
			return m_column_spec ;
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
                if( column_name == "qctool:cohort_name" ) {
                    return m_source_names[i] ;
                }
				else if( m_sources[i].check_for_column( column_name )) {
					result = m_sources[i].get_entry( sample_i, column_name ) ;
				}
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

