#ifndef GEN_TOOLS_InputToOutputFilenameMapper_HPP
#define GEN_TOOLS_InputToOutputFilenameMapper_HPP

#include <string>
#include <map>
#include <vector>
#include <cassert>
#include "wildcard.hpp"

class InputToOutputFilenameMapper
{
public:
	void add_filename_pair( std::string const& path_to_existing_files, std::string const& filename_template ) ;
	void add_filename_pairs( std::vector< std::string > const& paths_to_existing_files, std::vector< std::string > const& filename_templates ) ;


	std::size_t number_of_input_files() const { return m_existing_files.size() ; }
	std::size_t number_of_output_filenames() const { return m_constructed_filenames.size() ; }

	std::string const& input_file( std::size_t i ) const {
		return m_existing_files[i].filename() ;
	}

	std::string const& output_filename( std::size_t i ) const {
		return m_constructed_filenames[i] ;
	}

	std::string const& matched_part( std::size_t i ) const {
		return m_existing_files[i].match() ;
	}

	std::vector< std::string > input_files() const {
		std::vector< std::string > result( m_existing_files.size() ) ;
		for( std::size_t i = 0; i < result.size(); ++i ) {
			result[i] = m_existing_files[i].filename() ;
		}
		return result ;
	}

	std::vector< std::string > const& output_filenames() const {
		return m_constructed_filenames ;
	}

	typedef std::map< std::size_t, std::size_t > FilenameCorrespondence ;
	std::size_t filename_corresponding_to( std::size_t index ) const ;

private:

	std::vector< wildcard::FilenameMatch > m_existing_files ;
	std::vector< std::string > m_constructed_filenames ;
	std::map< std::size_t, std::size_t > m_filename_correspondence ;
} ;


#endif
