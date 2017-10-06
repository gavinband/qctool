
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include "config/config.hpp"
#include <boost/optional.hpp>
#include <boost/function.hpp>
#include "genfile/SNPDataSource.hpp"
#include "genfile/ReorderingSNPDataSource.hpp"
#include "genfile/VariantDataReader.hpp"

// #define DEBUG_REORDERING_SNP_DATA_SOURCE 1

namespace genfile {
	namespace {
		struct ReorderingVariantDataReader: public VariantDataReader {
			struct SampleSetter: public VariantDataReader::PerSampleSetter {

				enum EntryTypes { eMissing = 0, eString = 1, eInteger = 2, eDouble = 3 } ;

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

				void initialise( std::size_t nSamples, std::size_t nAlleles ) {
					m_number_of_alleles = nAlleles ;
					m_ploidies.resize( nSamples, 0 ) ;
					m_types.resize( nSamples ) ;
					m_entry_types.resize( nSamples ) ;
					m_ints.resize( nSamples ) ;
					m_strings.resize( nSamples ) ;
					m_doubles.resize( nSamples ) ;
				} ;
				bool set_sample( std::size_t i ) {
					m_storage_i = m_order[i] ;
					return true ;
				}
				void set_number_of_entries( uint32_t ploidy, std::size_t n, OrderType const order_type, ValueType const value_type ) {
					m_ploidies[ m_storage_i ] = ploidy ;
					m_types[ m_storage_i ] = std::make_pair( order_type, value_type ) ;
					m_entry_types[ m_storage_i ].clear() ;	
					m_entry_types[ m_storage_i ].resize( n, eMissing ) ;
					m_ints[ m_storage_i ].resize( n ) ;	
					m_strings[ m_storage_i ].resize( n ) ;	
					m_doubles[ m_storage_i ].resize( n ) ;	
				}
				void set_value( std::size_t i, MissingValue const value ) {
					m_entry_types[ m_storage_i ][i] = eMissing ;
				}
				void set_value( std::size_t i, std::string& value ) {
					m_entry_types[ m_storage_i ][i] = eString ;
					m_strings[ m_storage_i ][i] = value ;
				}
				void set_value( std::size_t i, Integer const value ) {
					m_entry_types[ m_storage_i ][i] = eInteger ;
					m_ints[ m_storage_i ][i] = value ;
				}
				void set_value( std::size_t i, double const value ) {
					m_entry_types[ m_storage_i ][i] = eDouble ;
					m_doubles[ m_storage_i ][i] = value ;
				}
				void finalise()  {
					m_setter.initialise( m_entry_types.size(), m_number_of_alleles ) ;
					for( std::size_t i = 0; i < m_entry_types.size(); ++i ) {
						if( m_setter.set_sample( i ) ) {
							m_setter.set_number_of_entries( m_ploidies[i], m_entry_types[i].size(), m_types[i].first, m_types[i].second ) ;
							for( std::size_t j = 0; j < m_entry_types[i].size(); ++j ) {
								switch( m_entry_types[i][j] ) {
									case eMissing: m_setter.set_value( j, MissingValue() ) ; break ;
									case eString: m_setter.set_value( j, m_strings[i][j] ) ; break ;
									case eInteger: m_setter.set_value( j, m_ints[i][j] ) ; break ;
									case eDouble: m_setter.set_value( j, m_doubles[i][j] ) ; break ;
									default: assert(0) ;
								}
							}
						}
					}
					m_setter.finalise() ;
				}
			
			private:
				VariantDataReader::PerSampleSetter& m_setter ;
				std::size_t m_number_of_alleles ;
				std::vector< std::size_t > const& m_order ;
				std::vector< uint32_t > m_ploidies ;
				std::vector< std::pair< OrderType, ValueType > > m_types ;
				std::vector< std::vector< char > > m_entry_types ;
				std::vector< std::vector< Integer > > m_ints ;
				std::vector< std::vector< std::string > > m_strings ;
				std::vector< std::vector< double > > m_doubles ;
				std::size_t m_storage_i ;
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
