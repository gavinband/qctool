#include "genfile/vcf/CallReader.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/Error.hpp"

namespace genfile {
	namespace vcf {
		CallReader::CallReader(
			std::string const& format,
			std::string const& data,
			boost::ptr_map< std::string, VCFEntryType > const& entry_types
		):
			m_format_elts( string_utils::split( format, ":" )),
			m_data( data ),
			m_entry_types( entry_types ),
			m_entries_by_position( get_entries_by_position( m_format_elts, entry_types ))
		{}

		// Return a vector of pointers to VCFEntryType objects in the supplied map,
		// with the ith VCFEntryType appropriate for the ith format element.
		std::vector< VCFEntryType const* > CallReader::get_entries_by_position(
			std::vector< std::string > const& format_elts,
			boost::ptr_map< std::string, VCFEntryType > const& entry_types
		) const
		{
			typedef boost::ptr_map< std::string, VCFEntryType >::const_iterator TypeIterator ;
			std::vector< VCFEntryType const* > result( format_elts.size() ) ;
			for( std::size_t i = 0; i < format_elts.size(); ++i ) {
				TypeIterator where = entry_types.find( format_elts[i] ) ;
				if( where == entry_types.end() ) {
					throw BadArgumentError( "genfile::vcf::CallReader::get_entries_in_format()", "format_elts" ) ;
				}
				result[i] = where->second ;
			}
			return result ;
		}

		CallReader& CallReader::operator()( std::string const& spec, Setter setter ) {
			std::vector< std::string >::const_iterator where = std::find( m_format_elts.begin(), m_format_elts.end(), spec ) ;
			if( where == m_format_elts.end() ) {
				throw BadArgumentError( "genfile::vcf::CallReader::operator()", "spec = \"" + spec + "\"" ) ;
			}
			else if( m_setters.find( where - m_format_elts.begin() ) != m_setters.end() ) {
				throw DuplicateKeyError( "genfile::vcf::CallReader::operator()", "spec = \"" + spec + "\"" ) ;
			}
			m_setters.insert( std::make_pair( where - m_format_elts.begin(), setter )) ;
			return *this ;
		}

		CallReader::~CallReader() {
			set_values( string_utils::split_and_strip( m_data, "\t", "\r\n" ), m_setters ) ;
		}

		void CallReader::set_values( std::vector< std::string > const& elts, Setters const& setters ) const {
			for( std::size_t i = 0; i < elts.size(); ++i ) {
				set_values( i, elts[i], m_setters ) ;
			}
		}

		void CallReader::set_values( std::size_t individual_i, std::string const& elt, Setters const& setters ) const {
			std::vector< std::string > components = string_utils::split( elt, ":" ) ;
			if( components.size() > m_entries_by_position.size() ) {
				throw MalformedInputError( "(unknown source)", 0 ) ;
			}
			Setters::const_iterator i = m_setters.begin(), end_i = m_setters.end() ;
			for( ; i != end_i; ++i ) {
				std::size_t element_pos = i->first ;
				if( element_pos >= components.size() ) {
					i->second.operator()( individual_i, m_entries_by_position[ element_pos ]->missing() ) ;
				}
				assert( element_pos < m_entries_by_position.size() ) ;
				i->second.operator()( individual_i, m_entries_by_position[ element_pos ]->parse( components[ element_pos ] )) ;
			}
		}
	}
}
