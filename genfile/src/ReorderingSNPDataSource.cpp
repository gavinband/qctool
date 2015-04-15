
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include "../config.hpp"
#include <boost/optional.hpp>
#include <boost/function.hpp>
#include "genfile/SNPDataSource.hpp"
#include "genfile/ReorderingSNPDataSource.hpp"
#include "genfile/VariantDataReader.hpp"

#define DEBUG_REORDERING_SNP_DATA_SOURCE 1

namespace genfile {
	namespace {
		struct ReorderingVariantDataReader: public VariantDataReader {
			struct SampleSetter: public VariantDataReader::PerSampleSetter {
				SampleSetter(
					VariantDataReader::PerSampleSetter& setter,
					std::vector< std::size_t > const& order
				):
					m_setter( setter ),
					m_order( order )
				{
#if DEBUG_REORDERING_SNP_DATA_SOURCE
					std::cerr << "ReorderingVariantDataReader::SampleSetter::reordering:\n" ;
					for( std::size_t i = 0; i < m_order.size(); ++i ) {
						std::cerr << " -- source " << i << " -> " << m_order[i] << ".\n" ;
					}
#endif					
				}

				void set_number_of_samples( std::size_t n ) {
					m_setter.set_number_of_samples( n ) ;
				} ;
				void set_number_of_alleles( std::size_t n ) {
					m_setter.set_number_of_alleles( n ) ;
				}
				bool set_sample( std::size_t i ) {
					return m_setter.set_sample( m_order[i] ) ;
				}
				void set_number_of_entries( std::size_t n ) {
					return m_setter.set_number_of_entries( n ) ;
				}
				void set_order_type( OrderType const order_type, ValueType const value_type ) {
					m_setter.set_order_type( order_type, value_type ) ;
				} ;
				void operator()( MissingValue const value ) {
					m_setter( value ) ;
				}
				void operator()( std::string& value ) {
					m_setter( value ) ;
				}
				void operator()( Integer const value ) {
					m_setter( value ) ;
				}
				void operator()( double const value ) {
					m_setter( value ) ;
				}
					
			private:
				VariantDataReader::PerSampleSetter& m_setter ;
				std::vector< std::size_t > const& m_order ;
			} ;
			
			ReorderingVariantDataReader( VariantDataReader::UniquePtr reader, std::vector< std::size_t > const& order ):
				m_reader( reader ),
				m_order( order )
			{
				// Don't do any order checks, we assume order is an ordering.
			}

			VariantDataReader& get( std::string const& spec, VariantDataReader::PerSampleSetter& setter ) {
				m_reader->get( spec, SampleSetter( setter, m_order ) ) ;
				// Do something here.
				return *this ;
			}
			bool supports( std::string const& spec ) const {
				return m_reader->supports( spec ) ;
			} ;
			void get_supported_specs( SpecSetter setter ) const {
				return m_reader->get_supported_specs( setter ) ;
			} ;
			std::size_t get_number_of_samples() const {
				return m_reader->get_number_of_samples() ;
			}

		private:
			VariantDataReader::UniquePtr m_reader ;
			// If order[i] = j then this indicates
			// that sample i in result dataset is sample j in the original dataset.
			std::vector< std::size_t > const& m_order ;
		} ;
	}

	VariantDataReader::UniquePtr ReorderingSNPDataSource::read_variant_data_impl() {
		return VariantDataReader::UniquePtr(
			new ReorderingVariantDataReader( m_source->read_variant_data(), m_inverse_order )
		) ;
	}
}
