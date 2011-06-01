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

		CallReader& CallReader::operator()( std::string const& spec, Setter setter ) {
			boost::ptr_map< std::string, VCFEntryType >::const_iterator entry_type_i = m_entry_types.find( spec ) ;
			if( entry_type_i == m_entry_types.end() ) {
				throw BadArgumentError( "genfile::vcf::CallReader::operator()", "spec = \"" + spec + "\"" ) ;
			}
			std::vector< std::string >::const_iterator where = std::find( m_format_elts.begin(), m_format_elts.end(), spec ) ;
			m_setters.insert( std::make_pair( where - m_format_elts.begin(), std::make_pair( setter, entry_type_i->second ))) ;
			return *this ;
		}

		CallReader::~CallReader() {
			// Empty data string means one individual with no data specified.
			// It is potentially valid if there's no GT field.
			set_values( string_utils::slice( m_data ).split( "\t" ), m_setters ) ;
		}

		void CallReader::set_values( std::vector< string_utils::slice > const& elts, Setters const& setters ) const {
			if( elts.size() != m_number_of_samples ) {
				throw MalformedInputError( "(data)", 0 ) ;
			}
			for( std::size_t i = 0; i < elts.size(); ++i ) {
				set_values( i, elts[i], m_setters ) ;
			}
		}
		
		void CallReader::set_values( std::size_t individual_i, string_utils::slice const& elt, Setters const& setters ) const {
			std::vector< string_utils::slice > components ;
			// Special case: if elt is blank, treat as no data at all for the individual.
			if( elt.size() > 0 ) {
				components = elt.split( ":" ) ;
			}
			if( components.size() > m_format_elts.size() ) {
				throw MalformedInputError( "(data)", 0, individual_i ) ;
			}
			
			try {
				// We will parse the genotype calls, if present, lazily, i.e when and if we need them.
				// We need them if any type needs the ploidy or if the GT field is requested.
				GenotypeCallVCFEntryType genotype_entry_type ;
				// Find the GT field in the format string: must either be 0 or one-past-the-end if not present.
				std::size_t const GT_field_pos = std::find( m_format_elts.begin(), m_format_elts.end(), "GT" ) - m_format_elts.begin() ;
				assert( GT_field_pos == 0 || GT_field_pos == m_format_elts.size() ) ;
				std::vector< Entry > genotype_calls ;
				std::size_t ploidy = 0 ;

				Setters::const_iterator i = m_setters.begin(), end_i = m_setters.end() ;
				for( ; i != end_i; ++i ) {
					std::size_t element_pos = i->first ;
					bool const elt_is_in_format = ( element_pos < m_format_elts.size() ) ;
					if( elt_is_in_format ) {
						bool const elt_is_trailing = ( element_pos >= components.size() ) ;
						bool const need_genotype_calls =
							( elt_is_in_format && m_format_elts[ element_pos ] == "GT" )
							|| i->second.second->check_if_requires_ploidy() ;
					
						// parse genotype call if either 1. the GT field is requested or 2. the requested field requires
						// ploidy information to parse.
						if( need_genotype_calls ) {
							// GT field must be present in the format and in the data
							if( GT_field_pos != 0 || components.size() == 0 ) {
								throw MalformedInputError( "(data)", 0, individual_i ) ;
							}
							genotype_calls = genotype_entry_type.parse( components[ 0 ], m_number_of_alleles ) ;
							ploidy = genotype_calls.size() ;
						}

						if( m_format_elts[ element_pos ] == "GT" ) {
							assert( !elt_is_trailing ) ;
							i->second.first.operator()( individual_i, genotype_calls ) ;
						}
						else if( i->second.second->check_if_requires_ploidy() ) {
							if( elt_is_trailing ) {
								i->second.first.operator()( individual_i, i->second.second->get_missing_value( m_number_of_alleles, ploidy ) ) ;
							}
							else {
								i->second.first.operator()( individual_i, i->second.second->parse( components[ element_pos ], m_number_of_alleles, ploidy )) ;
							}
						}
						else {
							if( elt_is_trailing ) {
								i->second.first.operator()( individual_i, i->second.second->get_missing_value( m_number_of_alleles ) ) ;
							}
							else {
								i->second.first.operator()( individual_i, i->second.second->parse( components[ element_pos ], m_number_of_alleles )) ;
							}
						}
					}
				}
			}
			catch( string_utils::StringConversionError const& ) {
				throw MalformedInputError( "(data)", 0, individual_i ) ;
			}
			catch( BadArgumentError const& e ) {
				throw MalformedInputError( "(data)", 0, individual_i ) ;
			}
		}
	}
}
