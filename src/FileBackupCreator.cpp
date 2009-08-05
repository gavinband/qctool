#include <string>

#include "FileUtil.hpp"
#include "FileBackupCreator.hpp"

void FileBackupCreator::backup_file_if_necessary( std::string const& filename ) {
	std::string backup_filename ;
	backup_filename = tmpnam(0) ;
	backup_file_if_necessary( filename, backup_filename ) ;
}

void FileBackupCreator::backup_file_if_necessary( std::string const& filename, std::string const& backup_filename ) {
	if( filename != "" && exists( filename )) {
		rename( filename, backup_filename ) ;
		m_backed_up_files[ filename ] = backup_filename ;
	}
}

