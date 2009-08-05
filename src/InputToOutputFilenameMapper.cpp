#include <string>
#include <map>
#include <vector>
#include "InputToOutputFilenameMapper.hpp"
#include "wildcard.hpp"

void InputToOutputFilenameMapper::add_filename_pair( std::string const& path_to_existing_files, std::string const& filename_template ) {
	std::vector< wildcard::FilenameMatch >
		existing_files = wildcard::find_files_matching_path_with_integer_wildcard( path_to_existing_files, '#', 1, 100 ) ;

	if( filename_template == "" ) {
		for( std::size_t j = 0; j < existing_files.size(); ++j ) {
			m_names_of_existing_files.push_back( existing_files[j].filename() ) ;
		}
	}
	else {
		std::vector< wildcard::FilenameMatch >
			constructed_filenames = wildcard::construct_corresponding_filenames( existing_files, filename_template, '#' ) ;
		assert( constructed_filenames.size() == existing_files.size() ) ;

		// Construct our filename data structures
		std::string last_constructed_filename = "" ;
		for( std::size_t j = 0; j < existing_files.size(); ++j ) {
			m_names_of_existing_files.push_back( existing_files[j].filename() ) ;
			if( constructed_filenames[j].filename() != last_constructed_filename ) {
				m_constructed_filenames.push_back( constructed_filenames[j].filename() ) ;
				last_constructed_filename = constructed_filenames[j].filename() ;
			}
			m_filename_correspondence[ m_names_of_existing_files.size() - 1 ] = m_constructed_filenames.size() - 1;
		}
	}
}

std::size_t InputToOutputFilenameMapper::filename_corresponding_to( std::size_t index ) const {
	return m_filename_correspondence.find( index )->second ;
}
