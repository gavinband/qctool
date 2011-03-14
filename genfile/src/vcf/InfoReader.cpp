#include <map>
#include <string>
#include <boost/ptr_container/ptr_map.hpp>
#include "genfile/vcf/InfoReader.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/Error.hpp"

namespace genfile {
	namespace vcf {
		namespace impl {
			std::map< std::string, std::string > parse_data( std::string const& data ) {
				std::vector< std::string > elts = string_utils::split( data, ";" ) ;
				std::map< std::string, std::string > result ;
				for( std::size_t i = 0; i < elts.size(); ++i ) {
					std::map< std::string, std::string >::const_iterator where = result.find( elts[i] ) ;
					if( where != result.end() ) {
						throw BadArgumentError( "genfile::vcf::InfoReader::parse_data()", "data = \"" + data + "\"" ) ;
					}
					std::size_t pos = elts[i].find( '=' ) ;
					if( pos == std::string::npos ) {
						result[ elts[i] ] = "(no value supplied)" ;
					}
					else {
						result[ elts[i].substr( 0, pos ) ] = elts[i].substr( pos+1, elts[i].size() ) ;
					}
				}
				return result ;
			}
		}

		InfoReader::InfoReader(
			std::size_t number_of_alleles,
			std::string const& data,
			boost::ptr_map< std::string, VCFEntryType > const& entry_types
		):
			m_number_of_alleles( number_of_alleles ),
			m_data( impl::parse_data( data ) ),
			m_entry_types( entry_types )
		{}
			
		InfoReader& InfoReader::operator()( std::string const& spec, Setter setter ) {
			if( !m_setters.insert( std::make_pair( spec, setter )).second ) {
				throw BadArgumentError( "genfile::vcf::InfoReader::operator()", "spec = \"" + spec + "\"" ) ;
			}
			return *this ;
		}

		InfoReader::~InfoReader() {
			set_values() ;
		}
		
		void InfoReader::set_values() const {
			for( Setters::const_iterator i = m_setters.begin(); i != m_setters.end(); ++i ) {
				std::map< std::string, std::string >::const_iterator where = m_data.find( i->first ) ;
				boost::ptr_map< std::string, VCFEntryType >::const_iterator entry_type_i = m_entry_types.find( i->first ) ;
				assert( entry_type_i != m_entry_types.end() ) ;
				if( where == m_data.end() ) {
					i->second( entry_type_i->second->get_missing_value( m_number_of_alleles ) ) ;
				}
				else {
					i->second( entry_type_i->second->parse( where->second, m_number_of_alleles )) ;
				}
			}
		}
	}
}