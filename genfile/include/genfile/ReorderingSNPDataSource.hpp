
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
#include <limits>
#include "config/config.hpp"
#include <boost/optional.hpp>
#include <boost/function.hpp>
#include "genfile/SNPDataSource.hpp"

namespace genfile {

	// A SNPDataSource that represents a view of an underlying source but with a specified
	// re-ordering of samples.
	class ReorderingSNPDataSource: public SNPDataSource
	{
	public:
		static std::size_t const eNotIncluded ;
		
		static SNPDataSource::UniquePtr create(
			SNPDataSource::UniquePtr source,
			std::vector< std::size_t > const& order
		) ;

	public:
		/*
		* Present a view of a data source with samples reordered according to the 
		* second argument.
		* Argument 'order' is a mapping from target to source samples
		* i.e. if order[j] = v_j, then v_j is the index of target sample
		* j in the source data.  If target sample j is not in the source data
		* it should be given value eNotIncluded.
		*/
		ReorderingSNPDataSource(
			SNPDataSource::UniquePtr source,
			std::vector< std::size_t > const& order
		) ;

		Metadata get_metadata() const {
			return m_source->get_metadata() ;
		} ;
		operator bool() const {
			return m_source->operator bool() ;
		} ;
		unsigned int number_of_samples() const {
			return m_order.size() ;
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
		// If order[j] = v_j then this indicates
		// that sample j in result dataset is sample v_j in the original dataset.
		std::vector< std::size_t > m_order ;
		// Inverse order, stored for efficient computation.
		std::vector< std::size_t > m_inverse_order ;
		
		
		// Used in sample mapping
		std::vector< uint32_t > m_ploidies ;
		std::vector< std::pair< OrderType, ValueType > > m_types ;
		std::vector< std::vector< char > > m_entry_types ;
		std::vector< std::vector< Integer > > m_ints ;
		std::vector< std::vector< std::string > > m_strings ;
		std::vector< std::vector< double > > m_doubles ;
		
		friend struct ReorderingVariantDataReader ;
		
	} ;
}

#endif
