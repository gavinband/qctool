#ifndef GEN_TOOLS_InputToOutputFilenameMapper_HPP
#define GEN_TOOLS_InputToOutputFilenameMapper_HPP

#include <string>
#include <map>
#include <vector>
#include <cassert>

class InputToOutputFilenameMapper
{
public:
	void add_filename_pair( std::string const& path_to_existing_files, std::string const& filename_template ) ;

	std::vector< std::string > const& input_files() const { return m_names_of_existing_files ; }
	std::vector< std::string > const& output_filenames() const { return m_constructed_filenames ; }
	typedef std::map< std::size_t, std::size_t > FilenameCorrespondence ;
	std::size_t index_of_filename_corresponding_to( std::size_t index ) const ;

private:

	std::vector< std::string > m_paths_to_existing_files ;
	std::vector< std::string > m_templates ;
	
	std::vector< std::string > m_names_of_existing_files ;
	std::vector< std::string > m_constructed_filenames ;
	std::map< std::size_t, std::size_t > m_filename_correspondence ;
} ;


#endif
