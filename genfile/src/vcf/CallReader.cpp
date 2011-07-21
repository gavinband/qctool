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
				m_genotype_calls.resize( m_number_of_samples ) ;
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
				set_values( m_components, ( where - m_format_elts.begin()), *(entry_type_i->second), setter ) ;
			}

			return *this ;
		}

		void CallReader::set_values(
			std::vector< std::vector< string_utils::slice > > const& elts,
			std::size_t const field_i,
			VCFEntryType const& entry_type,
			Setter const& setter
		) {
			assert( elts.size() == m_number_of_samples ) ;
			assert( field_i < m_format_elts.size() ) ;
			for( std::size_t sample_i = 0; sample_i < elts.size(); ++sample_i ) {
				set_values( sample_i, elts[ sample_i ], field_i, entry_type, setter ) ;
			}
		}
		
		void CallReader::set_values(
				std::size_t sample_i,
				std::vector< string_utils::slice > const& components,
				std::size_t element_pos,
				VCFEntryType const& entry_type,
				Setter const& setter
		) {
			try {
				unsafe_set_values(
					sample_i,
					components,
					element_pos,
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
			std::size_t element_pos,
			VCFEntryType const& entry_type,
			Setter const& setter
		) {
			// parse genotype call if either 1. the GT field is requested or 2. the requested field requires
			// ploidy information to parse.
			bool const need_genotype_calls = ( m_format_elts[ element_pos ] == "GT" || entry_type.check_if_requires_ploidy() ) ;

			if( need_genotype_calls && m_genotype_calls[ sample_i ].empty() ) {
				// Find the GT field in the format string: must either be 0 or one-past-the-end if not present.
				std::size_t const GT_field_pos = std::find( m_format_elts.begin(), m_format_elts.end(), "GT" ) - m_format_elts.begin() ;
				assert( GT_field_pos == 0 || GT_field_pos == m_format_elts.size() ) ;
				// GT field must be present in the format and in the data
				if( GT_field_pos != 0 || components.size() == 0 ) {
					throw MalformedInputError( "(data)", 0, sample_i ) ;
				}
				m_genotype_calls[ sample_i ] = m_genotype_call_entry_type.parse( components[ 0 ], m_number_of_alleles ) ;
			}

			// Decide if element is trailing (so not specified).
			bool const elt_is_trailing = ( element_pos >= components.size() ) ;

			if( m_format_elts[ element_pos ] == "GT" ) {
				assert( !elt_is_trailing ) ;
				setter( sample_i, m_genotype_calls[ sample_i ] ) ;
			}
			else if( entry_type.check_if_requires_ploidy() ) {
				std::size_t ploidy = m_genotype_calls[ sample_i ].size() ;
				if( elt_is_trailing ) {
					setter( sample_i, entry_type.get_missing_value( m_number_of_alleles, ploidy ) ) ;
				}
				else {
					setter( sample_i, entry_type.parse( components[ element_pos ], m_number_of_alleles, ploidy )) ;
				}
			}
			else {
				if( elt_is_trailing ) {
					setter( sample_i, entry_type.get_missing_value( m_number_of_alleles ) ) ;
				}
				else {
					setter( sample_i, entry_type.parse( components[ element_pos ], m_number_of_alleles )) ;
				}
			}
		}
	}
}
