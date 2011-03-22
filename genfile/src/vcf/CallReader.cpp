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
				}
				if( result.empty() || result[0] != "GT" ) {
					throw BadArgumentError( "genfile::vcf::impl::parse_format()", "format = \"" + format + "\"" ) ;
					
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
			std::size_t number_of_alleles,
			std::string const& format,
			std::string const& data,
			boost::ptr_map< std::string, VCFEntryType > const& entry_types
		):
			m_number_of_alleles( number_of_alleles ),
			m_format_elts( impl::parse_format( format ) ),
			m_data( data ),
			m_entry_types( entry_types ),
			m_entries_by_position( impl::get_entries_by_position( m_format_elts, entry_types ))
		{
			if( m_number_of_alleles == 0 ) {
				throw BadArgumentError( "genfile::vcf::CallReader::CallReader()", "number_of_alleles = " + string_utils::to_string( number_of_alleles ) ) ;
			}

			m_genotype_entry_type.reset( new GenotypeCallVCFEntryType ) ;
		}

		CallReader& CallReader::operator()( std::string const& spec, Setter setter ) {
			boost::ptr_map< std::string, VCFEntryType >::const_iterator entry_type_i = m_entry_types.find( spec ) ;
			if( entry_type_i == m_entry_types.end() ) {
				throw BadArgumentError( "genfile::vcf::CallReader::operator()", "spec = \"" + spec + "\"" ) ;
			}
			std::vector< std::string >::const_iterator where = std::find( m_format_elts.begin(), m_format_elts.end(), spec ) ;
			m_setters.insert( std::make_pair( where - m_format_elts.begin(), setter )) ;
			return *this ;
		}

		CallReader::~CallReader() {
			set_values( string_utils::slice( m_data ).split( "\t" ), m_setters ) ;
		}

		void CallReader::set_values( std::vector< string_utils::slice > const& elts, Setters const& setters ) const {
			for( std::size_t i = 0; i < elts.size(); ++i ) {
				set_values( i, elts[i], m_setters ) ;
			}
		}
		
		void CallReader::set_values( std::size_t individual_i, string_utils::slice const& elt, Setters const& setters ) const {
			std::vector< string_utils::slice > components = elt.split( ":" ) ;
			if( components.empty() || components.size() > m_entries_by_position.size() ) {
				throw MalformedInputError( "(data)", 0, individual_i ) ;
			}
			
			// First read the genotype calls, which must be present.
			try {
				std::vector< Entry > genotypes = m_genotype_entry_type->parse( components[0], m_number_of_alleles ) ;
				std::size_t const ploidy = genotypes.size() ;
			
				Setters::const_iterator i = m_setters.begin(), end_i = m_setters.end() ;
				for( ; i != end_i; ++i ) {
					std::size_t element_pos = i->first ;
					if( element_pos >= components.size() ) {
						i->second.operator()( individual_i, m_entries_by_position[ element_pos ]->get_missing_value( m_number_of_alleles, ploidy ) ) ;
					}
					else {
						i->second.operator()( individual_i, m_entries_by_position[ element_pos ]->parse( components[ element_pos ], m_number_of_alleles, ploidy )) ;
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
