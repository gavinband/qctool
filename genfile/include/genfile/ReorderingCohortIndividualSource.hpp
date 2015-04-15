
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_REORDERING_COHORTINDIVIDUALSOURCE_HPP

#include <string>
#include <memory>
#include <iosfwd>
#include <set>
#include <algorithm>

#include <boost/variant/variant.hpp>
#include <boost/variant/get.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include "genfile/Error.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/MissingValue.hpp"
#include "genfile/VariantEntry.hpp"

namespace genfile {
	// Base class for classes which provide random-access view
	// to a set of samples, 
	class ReorderingCohortIndividualSource: public CohortIndividualSource
	{
	public:
		ReorderingCohortIndividualSource(
			CohortIndividualSource::UniquePtr source,
			std::vector< std::size_t > const& order,
			std::vector< std::string > const& ids
		): 
			m_source( source ),
			m_order( order ),
			m_ids( ids )
		{
			assert( m_source.get() ) ;
			assert( m_order.size() == m_source->size() ) ;
			assert( m_ids.size() == m_source->size() ) ;
			assert( *std::min_element( m_order.begin(), m_order.end() ) == 0 ) ;
			assert( *std::max_element( m_order.begin(), m_order.end() ) == ( m_order.size() - 1 ) ) ;
			assert( std::set< std::size_t > ( m_order.begin(), m_order.end() ).size() == m_order.size() ) ;
		}
		
		std::size_t get_number_of_individuals() const {
			return m_source->get_number_of_individuals() ;
		}

		ColumnSpec get_column_spec() const {
			ColumnSpec spec = m_source->get_column_spec() ;
			spec.add_column( "original_id", e_DISCRETE_COVARIATE ) ;
			spec.add_column( "original_index", e_DISCRETE_COVARIATE ) ;
			if( spec.check_for_column( "ID_2" )) {
				spec.remove( "ID_2" ) ;
			}
			return spec ;
		}

		bool check_for_column( std::string const& column_name ) const {
			return( column_name == "original_id" || column_name == "original_index" || m_source->check_for_column( column_name )) ;
		}

		Entry get_entry( std::size_t sample_i, std::string const& column_name ) const {
			std::size_t original_sample_i = m_order[ sample_i ] ;
			if( column_name == "original_id" ) {
				return m_source->get_entry( original_sample_i, "ID_1" ) ;
			} else if( column_name == "original_index" ) {
				return Entry::Integer( original_sample_i ) ;
			} else if( column_name == "ID_1" ) {
				return m_ids[ sample_i ] ;
			} else {
				return m_source->get_entry( original_sample_i, column_name ) ;
			}
		}

		CohortIndividualSource const& get_base_source() const {
			return m_source->get_base_source() ;
		}
		CohortIndividualSource const& get_parent_source() const {
			return *m_source ;
		}

		// method: get_source_spec()
		// get_source_spec() returns a human-readable specification for this source.
		std::string get_source_spec() const {
			return "reordered:" + m_source->get_source_spec() ;
		}
	private:
		
		CohortIndividualSource::UniquePtr m_source ;
		std::vector< std::size_t > m_order ;
		std::vector< std::string > m_ids ;
	} ;
}

#endif
