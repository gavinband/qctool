#ifndef GEN_TOOLS_FILE_BACKUP_CREATOR_HPP
#define GEN_TOOLS_FILE_BACKUP_CREATOR_HPP

#include <string>
#include <map>

#include "FileUtil.hpp"


class FileBackupCreator {
public:
	void backup_file_if_necessary( std::string const& filename ) ;
	std::map< std::string, std::string > const& backed_up_files() const { return m_backed_up_files ; }	
protected:
	void backup_file_if_necessary( std::string const& filename, std::string const& backup_filename ) ;
private:
	std::map< std::string, std::string > m_backed_up_files ;	
} ;

#endif
