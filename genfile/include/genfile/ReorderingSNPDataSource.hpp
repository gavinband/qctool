
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef REORDERINGSNPDATASOURCE_HPP
#define REORDERINGSNPDATASOURCE_HPP

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include "../config.hpp"
#include <boost/optional.hpp>
#include <boost/function.hpp>
#include "genfile/SNPDataSource.hpp"

namespace genfile {

	// A SNPDataSource that represents a view of an underlying source but with a specified
	// re-ordering of samples.
	class ReorderingSNPDataSource: public SNPDataSource
	{
	public:
		ReorderingSNPDataSource(
			SNPDataSource::UniquePtr source,
			std::vector< std::size_t > const& order
		):
			m_source( source ),
			m_order( order )
		{
			assert( m_source.get() ) ;
			assert( m_order.size() == m_source->number_of_samples() ) ;
			assert( *std::min_element( m_order.begin(), m_order.end() ) == 0 ) ;
			assert( *std::max_element( m_order.begin(), m_order.end() ) == ( m_order.size() - 1 ) ) ;
			assert( std::set< std::size_t > ( m_order.begin(), m_order.end() ).size() == m_order.size() ) ;
			m_inverse_order.resize( order.size() ) ;
			for( std::size_t i = 0; i < m_order.size(); ++i ) {
				m_inverse_order[ m_order[i] ] = i ;
			}
		}

		void set_expected_ploidy( GetPloidy getPloidy ) {
			m_source->set_expected_ploidy( getPloidy ) ;
		} ;
		Metadata get_metadata() const {
			return m_source->get_metadata() ;
		} ;
		operator bool() const {
			return m_source->operator bool() ;
		} ;
		unsigned int number_of_samples() const {
			return m_source->number_of_samples() ;
		};
		OptionalSnpCount total_number_of_snps() const {
			return m_source->total_number_of_snps() ;
		} ;
		std::string get_source_spec() const {
			return "reordered:" + m_source->get_source_spec() ;
		} ;
		SNPDataSource const& get_parent_source() const {
			return *m_source ;
		}
		SNPDataSource const& get_base_source() const {
			return m_source->get_base_source() ;
		}
		
		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const {
			return "reordered:" + m_source->get_summary( prefix, column_width ) ;
		}

	protected:

		void get_snp_identifying_data_impl( 
			VariantIdentifyingData* result
		) {
			m_source->get_snp_identifying_data( result ) ;
		} ;	

		VariantDataReader::UniquePtr read_variant_data_impl() ;

		void ignore_snp_probability_data_impl() {
			m_source->ignore_snp_probability_data() ;
		}
		void reset_to_start_impl() {
			m_source->reset_to_start() ;
		}

	private:
		SNPDataSource::UniquePtr m_source ;
		// If order[i] = j then this indicates
		// that sample i in result dataset is sample j in the original dataset.
		std::vector< std::size_t > m_order ;
		// We need to compute the inverse.
		std::vector< std::size_t > m_inverse_order ;
	} ;
}

#endif
