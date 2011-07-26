#include "genfile/vcf/CallReader.hpp"
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
			m_entries_by_position( impl::get_entries_by_position( m_format_elts, entry_types ))
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

		namespace impl {
			struct CallReaderGenotypeSetter {
				CallReaderGenotypeSetter( std::vector< std::vector< Entry > >& result ):
					m_result( result )
				{}
				
				void operator()( std::size_t sample_i, std::vector< Entry > const& values ) {
					assert( sample_i < m_result.size() ) ;
					m_result[sample_i] = values ;
				}
			private:
				std::vector< std::vector< Entry > >& m_result ;
			} ;
		}

		CallReader& CallReader::get( std::string const& spec, Setter setter ) {
			boost::ptr_map< std::string, VCFEntryType >::const_iterator entry_type_i = m_entry_types.find( spec ) ;
			if( entry_type_i == m_entry_types.end() ) {
				throw BadArgumentError( "genfile::vcf::CallReader::operator()", "spec = \"" + spec + "\"" ) ;
			}

			if( m_components.empty() ) {
				std::vector< string_utils::slice > elts = string_utils::slice( m_data ).split( "\t" ) ;
				if( elts.size() != m_number_of_samples ) {
					throw MalformedInputError( "(data)", 0 ) ;
				}
				m_components.resize( m_number_of_samples ) ;
				for( std::size_t sample_i = 0; sample_i < elts.size(); ++sample_i ) {
					if( elts[ sample_i ].size() > 0 ) {
						m_components[ sample_i ] = elts[ sample_i ].split( ":" ) ;
						if( m_components[ sample_i ].size() > m_format_elts.size() ) {
							throw MalformedInputError( "(data)", 0, sample_i ) ;
						}
					}
				}
			}
			std::vector< std::string >::const_iterator where = std::find( m_format_elts.begin(), m_format_elts.end(), spec ) ;

			if( where != m_format_elts.end() ) {
				// Find out if we need to parse the GT field...
				if(
					( spec == "GT" || entry_type_i->second->check_if_requires_ploidy() )
					&&
					( m_genotype_calls.size() != m_number_of_samples )
				) {
					m_genotype_calls.resize( m_number_of_samples ) ;
					// Find the GT field in the format string...
					std::size_t const GT_field_pos = std::find( m_format_elts.begin(), m_format_elts.end(), "GT" ) - m_format_elts.begin() ;
					// ...it is required to be the first field.
					if( GT_field_pos != 0 ) {
						throw MalformedInputError( "(data)", 0 ) ;
					}
					impl::CallReaderGenotypeSetter genotype_setter( m_genotype_calls ) ;
					for( std::size_t sample_i = 0; sample_i < m_components.size(); ++sample_i ) {
						m_genotype_calls[ sample_i ] = m_genotype_call_entry_type.parse( m_components[ sample_i ][ GT_field_pos ], m_number_of_alleles ) ;
					}
				}
				if( spec == "GT" ) {
					for( std::size_t sample_i = 0; sample_i < m_components.size(); ++sample_i ) {
						setter( sample_i, m_genotype_calls[ sample_i ] ) ;
					}
				}
				else {
					for( std::size_t sample_i = 0; sample_i < m_components.size(); ++sample_i ) {
						set_values( sample_i, m_components[ sample_i ], ( where - m_format_elts.begin()), *(entry_type_i->second), setter ) ;
					}
				}
			}

			return *this ;
		}

		void CallReader::set_values(
				std::size_t sample_i,
				std::vector< string_utils::slice > const& components,
				std::size_t field_i,
				VCFEntryType const& entry_type,
				Setter const& setter
		) {
			try {
				unsafe_set_values(
					sample_i,
					components,
					field_i,
					entry_type,
					setter
				) ;
			}
			catch( string_utils::StringConversionError const& ) {
				throw MalformedInputError( "(data)", 0, sample_i ) ;
			}
			catch( BadArgumentError const& e ) {
				throw MalformedInputError( "(data)", 0, sample_i ) ;
			}
		}

		void CallReader::unsafe_set_values(
			std::size_t sample_i,
			std::vector< string_utils::slice > const& components,
			std::size_t field_i,
			VCFEntryType const& entry_type,
			Setter const& setter
		) {
			// GT should be handled elsewhere.
			assert( m_format_elts[ field_i ] != "GT" ) ;
			// Decide if element is trailing (so not specified).
			bool const elt_is_trailing = ( field_i >= components.size() ) ;
			if( entry_type.check_if_requires_ploidy() ) {
				assert( m_genotype_calls.size() == m_number_of_samples ) ;
				std::size_t ploidy = m_genotype_calls[ sample_i ].size() ;
				if( elt_is_trailing ) {
					setter( sample_i, entry_type.get_missing_value( m_number_of_alleles, ploidy ) ) ;
				}
				else {
					setter( sample_i, entry_type.parse( components[ field_i ], m_number_of_alleles, ploidy )) ;
				}
			}
			else {
				if( elt_is_trailing ) {
					setter( sample_i, entry_type.get_missing_value( m_number_of_alleles ) ) ;
				}
				else {
					setter( sample_i, entry_type.parse( components[ field_i ], m_number_of_alleles )) ;
				}
			}
		}
	}
}
