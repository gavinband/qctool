#include <string>
#include <vector>
#include <cassert>
#include <iostream>
#include <set>
#include <algorithm>
#include "genfile/Pedigree.hpp"
#include "genfile/FromPedFilePedigree.hpp"
#include "genfile/snp_data_utils.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/Error.hpp"

namespace genfile {
	FromPedFilePedigree::FromPedFilePedigree( std::string const& filename ):
		m_filename( filename )
	{
		std::auto_ptr< std::istream > stream = open_text_file_for_input(
			filename,
			get_compression_type_indicated_by_filename( filename )
		) ;
		load_data_from_stream( *stream ) ;
	}

	FromPedFilePedigree::FromPedFilePedigree( std::istream& stream ):
		m_filename( "(unnamed stream)" )
	{
		load_data_from_stream( stream ) ;
	}

	void FromPedFilePedigree::load_data_from_stream( std::istream& stream ) {
		std::string line ;
		std::vector< std::string > individuals ;	
		std::map< std::string, std::string > families ;
		std::map< std::string, std::vector< std::string > > parents ;
		std::map< std::string, std::vector< std::string > > children ;
		std::map< std::string, Sex > sexes ;

		std::size_t line_count = 0 ;
		for( ; std::getline( stream, line ); ++line_count ) {
			std::vector< std::string > elts = string_utils::split_and_strip_discarding_empty_entries( line, " \t", "\r\t " ) ;
			// First 6 columns have the relevant info.
			// See http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml
			std::string const& id = elts[1] ;
			individuals.push_back( id ) ;
			families[ id ] = elts[0] ;

			if( elts[2] == "NA" || elts[2] == "?" || elts[2] == "0" ) {
				assert( elts[3] == elts[2] ) ;
				parents[ id ].push_back( "0" ) ;
				parents[ id ].push_back( "0" ) ;
			}
			else {
				if( elts[3] == "NA" || elts[3] == "?" || elts[3] == "0" ) {
					throw MalformedInputError( m_filename, line_count, 3 ) ;
				}
				parents[ id ].push_back( elts[2] ) ;
				parents[ id ].push_back( elts[3] ) ;

				children[ elts[2] ].push_back( id ) ;
				children[ elts[3] ].push_back( id ) ;
			}

			if( elts[4] == "1" || string_utils::to_lower( elts[4] ) == "m" ) {
				sexes[id] = eMale ;
			}
			else if( elts[4] == "2" || string_utils::to_lower( elts[4] ) == "f" ) {
				sexes[id] = eFemale ;
			}
			else if( elts[4] == "other" ) {
				sexes[id] = eUnknown ;
			}
			else {
				throw MalformedInputError( m_filename, line_count, 4 ) ;
			}
		}
		
		{
			std::sort( individuals.begin(), individuals.end() ) ;
			// check for uniqueness
			if( std::unique( individuals.begin(), individuals.end() ) != individuals.end() ) {
				throw MalformedInputError( m_filename, line_count ) ;
			}
		}
		
		for( std::map< std::string, std::vector< std::string > >::iterator i = parents.begin(); i != parents.end(); ++i ) {
			std::sort( i->second.begin(), i->second.end() ) ;
			assert( i->second.size() == 2 ) ;
			if( std::unique( i->second.begin(), i->second.end() ) != i->second.end() ) {
				assert( i->second[0] == "0" && i->second[1] == "0" ) ;
			}
		}

		for( std::map< std::string, std::vector< std::string > >::iterator i = children.begin(); i != children.end(); ++i ) {
			std::sort( i->second.begin(), i->second.end() ) ;
			assert( std::unique( i->second.begin(), i->second.end() ) == i->second.end() ) ;
		}
		
		m_individuals = individuals ;
		m_families = families ;
		m_parents = parents ;
		m_children = children ;
		m_sexes = sexes ;
	}

	std::size_t FromPedFilePedigree::get_number_of_individuals() const {
		return m_individuals.size() ;
	}

	std::string FromPedFilePedigree::get_id_of( std::size_t i ) const {
		assert( i < m_individuals.size() ) ;
		return m_individuals[i] ;
	}
	
	std::string FromPedFilePedigree::get_family_of( std::string const& individual ) const {
		return get_map_value( m_families, individual ) ;
	}

	std::vector< std::string > FromPedFilePedigree::get_parents_of( std::string const& individual ) const {
		return FromPedFilePedigree::get_map_value( m_parents, individual ) ;
	}

	std::vector< std::string > FromPedFilePedigree::get_children_of( std::string const& individual ) const {
		return get_map_value( m_children, individual ) ;
	}

	Pedigree::Sex FromPedFilePedigree::get_sex_of( std::string const& individual ) const {
		return get_map_value( m_sexes, individual ) ;
	}
	
	std::vector< std::string > FromPedFilePedigree::get_map_value(
		std::map< std::string, std::vector< std::string > > const& map,
		std::string const& key
	) {
		std::map< std::string, std::vector< std::string > >::const_iterator
			where = map.find( key ) ;
		assert( where != map.end() ) ;
		return where->second ;
	}

	std::string FromPedFilePedigree::get_map_value(
		std::map< std::string, std::string > const& map,
		std::string const& key
	) {
		std::map< std::string, std::string >::const_iterator where = map.find( key ) ;
		assert( where != map.end() ) ;
		return where->second ;
	}

	Pedigree::Sex FromPedFilePedigree::get_map_value(
		std::map< std::string, Sex > const& map,
		std::string const& key
	) {
		std::map< std::string, Sex >::const_iterator where = map.find( key ) ;
		assert( where != map.end() ) ;
		return where->second ;
	}
	
}
