
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <map>
#include <vector>
#include "genfile/snp_data_utils.hpp"
#include "genfile/wildcard.hpp"
#include "InputToOutputFilenameMapper.hpp"

void InputToOutputFilenameMapper::add_filename_pair( std::string const& path_to_existing_files, std::string const& filename_template ) {
	m_existing_files = genfile::wildcard::find_files_by_chromosome(
		path_to_existing_files,
		genfile::wildcard::eALL_CHROMOSOMES
	) ;

	if( m_existing_files.size() == 0 ) {
		throw genfile::FileNotFoundError( path_to_existing_files ) ;
	}

	if( filename_template != "" ) {
		std::vector< genfile::wildcard::FilenameMatch > constructed_filenames = genfile::wildcard::construct_corresponding_filenames(
			m_existing_files,
			filename_template,
			'#'
		) ;
		assert( constructed_filenames.size() == m_existing_files.size() ) ;

		// Construct the filename->filename mapping.
		std::string last_constructed_filename = "" ;
		for( std::size_t j = 0; j < m_existing_files.size(); ++j ) {
			if( constructed_filenames[j].filename() != last_constructed_filename ) {
				m_constructed_filenames.push_back( constructed_filenames[j].filename() ) ;
				last_constructed_filename = constructed_filenames[j].filename() ;
			}
			m_filename_correspondence[ j ] = m_constructed_filenames.size() - 1;
		}
	}
}

void InputToOutputFilenameMapper::add_filename_pairs( std::vector< std::string > const& paths_to_existing_files, std::vector< std::string > const& filename_templates ) {
	assert( paths_to_existing_files.size() == filename_templates.size() ) ;
	for( std::size_t i = 0; i < paths_to_existing_files.size(); ++i ) {
		add_filename_pair( paths_to_existing_files[i], filename_templates[i] ) ;
	}
}

std::size_t InputToOutputFilenameMapper::filename_corresponding_to( std::size_t index ) const {
	std::map< std::size_t, std::size_t >::const_iterator where = m_filename_correspondence.find( index ) ;
	assert( where != m_filename_correspondence.end() ) ;
	return where->second ;
}
