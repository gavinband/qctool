
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

	size_t const ReorderingSNPDataSource::eNotIncluded = std::numeric_limits< std::size_t >::max() ;

	SNPDataSource::UniquePtr ReorderingSNPDataSource::create(
		SNPDataSource::UniquePtr source,
		std::vector< std::size_t > const& order
	) {
		bool order_is_identity = true ;
		if( source->number_of_samples() == order.size() ) {
			for( std::size_t i = 0; i < order.size(); ++i ) {
				if( order[i] != i ) {
					order_is_identity = false ;
					break ;
				}
			}
		}
		
		if( order_is_identity ) {
			return source ;
		} else {
			return SNPDataSource::UniquePtr(
				new ReorderingSNPDataSource( source, order )
			) ;
		}
	}
	
	ReorderingSNPDataSource::ReorderingSNPDataSource(
		SNPDataSource::UniquePtr source,
		std::vector< std::size_t > const& order
	):
		m_source( source ),
		m_order( order ),
		m_inverse_order( m_source->number_of_samples(), eNotIncluded )
	{
		assert( m_source.get() ) ;
		// Check all values are eNotIncluded or in range
		// And check all included values are unique, i.e. no duplicating source samples.
		assert( *std::min_element( m_order.begin(), m_order.end() ) >= 0 ) ;
		std::set< std::size_t > uniqueSourceSamples ;
		std::size_t sampleCount = 0 ;
		for( std::size_t i = 0; i < m_order.size(); ++i ) {
			if( m_order[i] != eNotIncluded ) {
				assert( m_order[i] < m_source->number_of_samples() ) ;
				uniqueSourceSamples.insert( m_order[i] ) ;
				assert( uniqueSourceSamples.size() == ++sampleCount ) ;
			}
		}
		/* Compute the inverse ordering, i.e. mapping from source to target samples.
		Samples not in output are given the value eNotIncluded */
		for( std::size_t i = 0; i < m_order.size(); ++i ) {
			if( m_order[i] != eNotIncluded ) {
				m_inverse_order[ m_order[i] ] = i ;
			}
		}
		
		/* Create appropriate storage. */
		m_ploidies.resize( m_order.size() ) ;
		m_types.resize( m_order.size() ) ;
		m_entry_types.resize( m_order.size() ) ;
		m_ints.resize( m_order.size() ) ;
		m_strings.resize( m_order.size() ) ;
		m_doubles.resize( m_order.size() ) ;
	}
	
	struct ReorderingVariantDataReader: public VariantDataReader {
		struct SampleSetter: public VariantDataReader::PerSampleSetter {

			enum EntryTypes { eMissing = 0, eString = 1, eInteger = 2, eDouble = 3 } ;

			/*
			* source_to_target is a vector of indices (o_i)
			* where o_i is the index at which sample i will appear in the output.
			* Ordering must cover a contiguous range of indices.
			* Any value equal to max value of std::size_t will be ignored, not appear in output.
			* Any output sample not set is output as missing data.
			*/
			SampleSetter(
				VariantDataReader::PerSampleSetter& setter,
				std::vector< std::size_t > const& target_to_source,
				std::vector< std::size_t > const& source_to_target,
				std::vector< uint32_t >& ploidies,
				std::vector< std::pair< OrderType, ValueType > >& types,
				std::vector< std::vector< char > >& entry_types,
				std::vector< std::vector< Integer > >& ints,
				std::vector< std::vector< std::string > >& strings,
				std::vector< std::vector< double > >& doubles
			):
				m_setter( setter ),
				m_target_to_source( target_to_source ),
				m_source_to_target( source_to_target ),
				m_ploidies( ploidies ),
				m_types( types ),
				m_entry_types( entry_types ),
				m_ints( ints ),
				m_strings( strings ),
				m_doubles( doubles )
			{
#if DEBUG_REORDERING_SNP_DATA_SOURCE
				std::cerr << "ReorderingVariantDataReader::SampleSetter::reordering:\n" ;
				for( std::size_t i = 0; i < m_source_to_target.size(); ++i ) {
					std::cerr << " -- source " << i << " -> " << m_source_to_target[i] << ".\n" ;
				}
#endif					
			}

			void initialise( std::size_t nSamples, std::size_t nAlleles ) {
				m_number_of_alleles = nAlleles ;
				// Set ploidy to 0, types to completely missing, clear entry types for each sample
				std::fill( m_ploidies.begin(), m_ploidies.end(), 0 ) ;
				std::fill( m_types.begin(), m_types.end(), std::make_pair( eUnknownOrderType, eUnknownValueType ) ) ;
				for( std::size_t i = 0 ; i < m_entry_types.size(); ++i ) {
					m_entry_types[i].clear() ;
				}
			} ;
			bool set_sample( std::size_t i ) {
				m_storage_i = m_source_to_target[i] ;
				return (m_storage_i != ReorderingSNPDataSource::eNotIncluded) ;
			}
			void set_number_of_entries( uint32_t ploidy, std::size_t n, OrderType const order_type, ValueType const value_type ) {
				if( m_storage_i != ReorderingSNPDataSource::eNotIncluded ) {
					m_ploidies[ m_storage_i ] = ploidy ;
					m_types[ m_storage_i ] = std::make_pair( order_type, value_type ) ;
					m_entry_types[ m_storage_i ].clear() ;	
					m_entry_types[ m_storage_i ].resize( n, eMissing ) ;
					m_ints[ m_storage_i ].resize( n ) ;	
					m_strings[ m_storage_i ].resize( n ) ;	
					m_doubles[ m_storage_i ].resize( n ) ;	
				}
			}
			void set_value( std::size_t i, MissingValue const value ) {
				if( m_storage_i != ReorderingSNPDataSource::eNotIncluded ) {
					m_entry_types[ m_storage_i ][i] = eMissing ;
				}
			}
			void set_value( std::size_t i, std::string& value ) {
				if( m_storage_i != ReorderingSNPDataSource::eNotIncluded ) {
					m_entry_types[ m_storage_i ][i] = eString ;
					m_strings[ m_storage_i ][i] = value ;
				}
			}
			void set_value( std::size_t i, Integer const value ) {
				if( m_storage_i != ReorderingSNPDataSource::eNotIncluded ) {
					m_entry_types[ m_storage_i ][i] = eInteger ;
					m_ints[ m_storage_i ][i] = value ;
				}
			}
			void set_value( std::size_t i, double const value ) {
				if( m_storage_i != ReorderingSNPDataSource::eNotIncluded ) {
					m_entry_types[ m_storage_i ][i] = eDouble ;
					m_doubles[ m_storage_i ][i] = value ;
				}
			}
			void finalise() {
				m_setter.initialise( m_entry_types.size(), m_number_of_alleles ) ;
				for( std::size_t i = 0; i < m_target_to_source.size(); ++i ) {
					if( m_setter.set_sample( i ) ) {
						if( m_target_to_source[i] != ReorderingSNPDataSource::eNotIncluded ) {
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
				}
				m_setter.finalise() ;
			}

		private:
			VariantDataReader::PerSampleSetter& m_setter ;
			std::vector< std::size_t > const& m_target_to_source ;
			std::vector< std::size_t > const& m_source_to_target ;

			std::size_t m_number_of_alleles ;
			std::vector< uint32_t >& m_ploidies ;
			std::vector< std::pair< OrderType, ValueType > >& m_types ;
			std::vector< std::vector< char > >& m_entry_types ;
			std::vector< std::vector< Integer > >& m_ints ;
			std::vector< std::vector< std::string > >& m_strings ;
			std::vector< std::vector< double > >& m_doubles ;

			std::size_t m_storage_i ;
		} ;
	
		ReorderingVariantDataReader(
			ReorderingSNPDataSource& source,
			VariantDataReader::UniquePtr reader
		):
			m_source( source ),
			m_reader( reader )
		{
			// Don't do any order checks, we assume order is an ordering.
		}

		VariantDataReader& get( std::string const& spec, VariantDataReader::PerSampleSetter& setter ) {
			m_reader->get(
				spec,
				SampleSetter(
					setter,
					m_source.m_order,
					m_source.m_inverse_order,
					m_source.m_ploidies,
					m_source.m_types,
					m_source.m_entry_types,
					m_source.m_ints,
					m_source.m_strings,
					m_source.m_doubles
				)
			) ;
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
		ReorderingSNPDataSource& m_source ;
		VariantDataReader::UniquePtr m_reader ;
	} ;

	VariantDataReader::UniquePtr ReorderingSNPDataSource::read_variant_data_impl() {
		return VariantDataReader::UniquePtr(
			new ReorderingVariantDataReader(
				*this,
				m_source->read_variant_data()
			)
		) ;
	}
}
