
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "genfile/vcf/CallReader.hpp"
#include "genfile/vcf/Types.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/Error.hpp"

namespace genfile {
	namespace vcf {
		namespace impl {
			std::vector< std::string > parse_format( std::string const& format ) {
				std::vector< std::string > result = string_utils::split( format, ":" ) ;
				std::map< std::string, std::size_t > elt_counts ;
				for( std::size_t i = 0; i < result.size(); ++i ) {
					if( ++elt_counts[ result[i] ] > 1 ) {
						throw BadArgumentError( "genfile::vcf::impl::parse_format()", "format = \"" + format + "\"" ) ;
					}
					// GT must be first if specified
					if( result[i] == "GT" && i > 0 ) {
						throw BadArgumentError( "genfile::vcf::impl::parse_format()", "format = \"" + format + "\"" ) ;
					}
				}
				return result ;
			}

			// Return a vector of pointers to VCFEntryType objects in the supplied map,
			// with the ith VCFEntryType appropriate for the ith format element.
			std::vector< VCFEntryType const* > get_entries_by_position(
				std::vector< std::string > const& format_elts,
				boost::ptr_map< std::string, VCFEntryType > const& entry_types
			)
			{
				typedef boost::ptr_map< std::string, VCFEntryType >::const_iterator TypeIterator ;
				std::vector< VCFEntryType const* > result( format_elts.size() ) ;
				for( std::size_t i = 0; i < format_elts.size(); ++i ) {
					TypeIterator where = entry_types.find( format_elts[i] ) ;
					if( where == entry_types.end() ) {
						throw BadArgumentError( "genfile::vcf::impl::get_entries_in_format()", "format_elts" ) ;
					}
					result[i] = where->second ;
				}
				return result ;
			}
		}
		
		CallReader::CallReader(
			std::size_t number_of_samples,
			std::size_t number_of_alleles,
			std::string const& format,
			std::string const& data,
			boost::ptr_map< std::string, VCFEntryType > const& entry_types
		):
			m_number_of_samples( number_of_samples ),
			m_number_of_alleles( number_of_alleles ),
			m_format_elts( impl::parse_format( format ) ),
			m_data( data ),
			m_entry_types( entry_types ),
			m_entries_by_position( impl::get_entries_by_position( m_format_elts, entry_types )),
			m_strict_mode( true )
		{
			if( m_number_of_alleles == 0 ) {
				throw BadArgumentError( "genfile::vcf::CallReader::CallReader()", "number_of_alleles = " + string_utils::to_string( number_of_alleles ) ) ;
			}
			if( m_number_of_samples == 0 ) {
				throw BadArgumentError( "genfile::vcf::CallReader::CallReader()", "number_of_samples = " + string_utils::to_string( number_of_samples ) ) ;
			}
		}

		void CallReader::get_format_elts( boost::function< void ( std::string ) > setter ) const {
			for( std::size_t i = 0; i < m_format_elts.size(); ++i ) {
				setter( m_format_elts[i] ) ;
			}
		}

		void CallReader::set_strict_mode( bool value ) {
			m_strict_mode = value ;
		}

		namespace impl {
			struct CallReaderGenotypeSetter: public CallReader::Setter {
				CallReaderGenotypeSetter(
					std::vector< vcf::EntrySetter::Integer >& entries,
					std::vector< std::size_t >& ploidies
				):
					m_entries( entries ),
					m_ploidies( ploidies ),
					m_sample_i( 0 ),
					m_current_i( 0 ),
					m_order_type( eUnorderedList )
				{}

				void initialise( std::size_t nSamples, std::size_t nAlleles ) {
					assert( nAlleles == 2 ) ;
					m_ploidies.resize( nSamples ) ;
					m_entries.reserve( nSamples * nAlleles ) ;
				}

				bool set_sample( std::size_t sample_i ) {
					assert( sample_i < m_ploidies.size() ) ;
					m_sample_i = sample_i ;
					return true ;
				}

				void set_number_of_entries( uint32_t ploidy, std::size_t n, OrderType const order_type, ValueType const value_type ) {
					if( order_type != ePerOrderedHaplotype && order_type != ePerUnorderedHaplotype ) {
						throw BadArgumentError(
							"genfile::vcf::impl::CallReaderGenotypeSetter::set_number_of_entries()",
							"order_type",
							"Expected order_type = ePerOrderedHaplotype or ePerUnorderedHaplotype for genotype field."
						) ;
					}
					if( value_type != eAlleleIndex ) {
						throw BadArgumentError(
							"genfile::vcf::impl::CallReaderGenotypeSetter::set_ordeset_number_of_entriesr_type()",
							"value_type",
							"Expected value_type = eAlleleIndex for genotype field."
						) ;
					}
					if( ploidy != n ) {
						throw BadArgumentError(
							"genfile::vcf::impl::CallReaderGenotypeSetter::set_orderset_number_of_entries_type()",
							"n",
							"Expected n == ploidy"
						) ;
					}
					m_ploidies[ m_sample_i ] = n ;
					m_order_type = order_type ;
					m_entries.resize( m_current_i + n ) ;
				}
			
				void set_value( MissingValue const value ) {
					assert( m_current_i < m_entries.size() ) ;
					m_entries[ m_current_i++ ] = -1 ;
				}

				void set_value( Integer const value ) {
					assert( m_current_i < m_entries.size() ) ;
					m_entries[ m_current_i++ ] = value ;
				}
				
				OrderType get_order_type() const {
					return m_order_type ;
				}
				
				void finalise() {}

			private:
				std::vector< vcf::EntrySetter::Integer >& m_entries ;
				std::vector< std::size_t >& m_ploidies ;
				std::size_t m_sample_i ;
				std::size_t m_current_i ;
				OrderType m_order_type ;
			} ;
		}

		CallReader& CallReader::get( std::string const& spec, Setter& setter ) {
			boost::ptr_map< std::string, VCFEntryType >::const_iterator entry_type_i = m_entry_types.find( spec ) ;
			if( entry_type_i == m_entry_types.end() ) {
				throw BadArgumentError( "genfile::vcf::CallReader::operator()", "spec = \"" + spec + "\"" ) ;
			}
			
			if( m_components.empty() ) {
				split_data() ;
			}
			assert( m_component_counts.size() == m_number_of_samples ) ;

			std::vector< std::string >::const_iterator where = std::find( m_format_elts.begin(), m_format_elts.end(), spec ) ;

			
			if( where != m_format_elts.end() ) {
				// Find out if we need to parse the GT field...
				if(
					( spec == "GT" || entry_type_i->second->check_if_requires_ploidy() )
					&&
					( m_ploidy.size() != m_number_of_samples )
				) {
					load_genotypes() ;
				}

				if( spec == "GT" ) {
					assert( m_ploidy.size() > 0 ) ;
					std::size_t index = 0 ;
					setter.initialise( m_number_of_samples, m_number_of_alleles ) ;
					for( std::size_t sample_i = 0; sample_i < m_number_of_samples; ++sample_i ) {
						std::size_t const ploidy = m_ploidy[ sample_i ] ;
						assert( m_genotype_calls.size() >= ( index + ploidy ) ) ;
						setter.set_sample( sample_i ) ;
						setter.set_number_of_entries(
							ploidy, ploidy,
							m_order_types[ sample_i ],
							eAlleleIndex
						) ;
						for( std::size_t const current_index = index; index < ( current_index + ploidy ); ++index ) {
							if( m_genotype_calls[ index ] == -1 ) {
								setter.set_value( MissingValue() ) ;
							} else {
								setter.set_value( m_genotype_calls[ index ] ) ;
							}
						}
					}
				}
				else {
					setter.initialise( m_number_of_samples, m_number_of_alleles ) ;
					for(
						std::size_t sample_i = 0, component_index = 0;
						sample_i < m_number_of_samples;
						component_index += m_component_counts[ sample_i++ ]
					) {
						setter.set_sample( sample_i ) ;
						assert( m_components.size() >= ( component_index + m_component_counts[ sample_i ] ) ) ;
						set_values(
							sample_i,
							&m_components[ 0 ] + component_index,
							&m_components[ 0 ] + ( component_index + m_component_counts[ sample_i ] ),
							( where - m_format_elts.begin()),
							*(entry_type_i->second),
							setter
						) ;
					}
				}
			}

			return *this ;
		}

		void CallReader::split_data() {
			std::vector< string_utils::slice > elts ;
			elts.reserve( m_number_of_samples ) ;
			string_utils::slice( m_data ).split( "\t", &elts ) ;
			if( elts.size() != m_number_of_samples ) {
				throw MalformedInputError( "(data)", 0 ) ;
			}
			m_components.reserve( m_number_of_samples ) ;
			for( std::size_t sample_i = 0; sample_i < elts.size(); ++sample_i ) {
				std::size_t const current_size = m_components.size() ;
				if( elts[ sample_i ].size() > 0 ) {
					elts[ sample_i ].split( ":", &m_components ) ;
				}
				m_component_counts.push_back( m_components.size() - current_size ) ;
				if( m_component_counts.back() > m_format_elts.size() ) {
					throw MalformedInputError( "(data)", 0, sample_i ) ;
				}
			}
		}
		
		void CallReader::load_genotypes() {
			// Find the GT field in the format string...
			std::size_t const GT_field_pos = std::find( m_format_elts.begin(), m_format_elts.end(), "GT" ) - m_format_elts.begin() ;
			// ...it is required to be the first field.
			if( GT_field_pos != 0 ) {
				throw MalformedInputError( "(data)", 0 ) ;
			}

			m_ploidy.resize( m_number_of_samples ) ;
			m_order_types.resize( m_number_of_samples ) ;
			m_genotype_calls.reserve( m_number_of_samples * 2 ) ;
			std::size_t total_number_of_calls = 0 ;
			impl::CallReaderGenotypeSetter genotype_setter( m_genotype_calls, m_ploidy ) ;
			for(
				std::size_t sample_i = 0, component_index = 0;
				sample_i < m_number_of_samples;
				component_index += m_component_counts[ sample_i++ ]
			) {
				try {
					genotype_setter.set_sample( sample_i ) ;
					m_genotype_call_entry_type.parse( m_components[ component_index + GT_field_pos ], m_number_of_alleles, eUnknownPloidy, genotype_setter ) ;
					total_number_of_calls += m_ploidy[ sample_i ] ;
					m_order_types[ sample_i ] = genotype_setter.get_order_type() ;
				}
				catch( string_utils::StringConversionError const& ) {
					throw MalformedInputError( "(data)", 0, sample_i ) ;
				}
				catch( BadArgumentError const& e ) {
					throw MalformedInputError( "(data)", 0, sample_i ) ;
				}
			}
			assert( m_genotype_calls.size() == total_number_of_calls ) ;
		}
		
		void CallReader::set_values(
				std::size_t sample_i,
				string_utils::slice const* begin_components,
				string_utils::slice const* end_components,
				std::size_t field_i,
				VCFEntryType const& entry_type,
				Setter& setter
		) {
			try {
				unsafe_set_values(
					sample_i,
					begin_components,
					end_components,
					field_i,
					entry_type,
					setter
				) ;
			}
			catch( string_utils::StringConversionError const& ) {
				throw MalformedInputError( "(data)", 0, sample_i ) ;	
			}
			catch( BadArgumentError const& ) {
				throw MalformedInputError( "(data)", 0, sample_i ) ;
			}
		}

		void CallReader::unsafe_set_values(
			std::size_t sample_i,
			string_utils::slice const* begin_components,
			string_utils::slice const* end_components,
			std::size_t field_i,
			VCFEntryType const& entry_type,
			Setter& setter
		) {
			// GT should be handled elsewhere.
			assert( m_format_elts[ field_i ] != "GT" ) ;
			// Decide if element is trailing (so not specified).
			bool const elt_is_trailing = (( begin_components + field_i ) >= end_components ) ;
			// std::cerr << "CallReader::unsafe_set_values(): sample_i=" << sample_i << ", field_i=" << field_i << ".\n" ;
			uint32_t ploidy = eUnknownPloidy ;
			if( entry_type.check_if_requires_ploidy() ) {
				assert( m_ploidy.size() == m_number_of_samples ) ;
				ploidy = m_ploidy[ sample_i ] ;
			}
			if( elt_is_trailing ) {
				entry_type.get_missing_value( m_number_of_alleles, ploidy, setter ) ;
			}
			else {
				try {
					entry_type.parse( *(begin_components + field_i ), m_number_of_alleles, ploidy, setter ) ;
				}
				catch( BadArgumentError const& ) {
					if( m_strict_mode ) {
						throw ;
					} else {
						// parse error, set missing value.
						entry_type.get_missing_value( m_number_of_alleles, ploidy, setter ) ;
					}
				}
			}
		}
	}
}
