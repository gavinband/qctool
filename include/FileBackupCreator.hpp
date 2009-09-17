#ifndef GEN_TOOLS_FILE_BACKUP_CREATOR_HPP
#define GEN_TOOLS_FILE_BACKUP_CREATOR_HPP

#include <string>
#include <map>

#include "FileUtil.hpp"


class FileBackupCreator {
public:
	~FileBackupCreator() {};
	void backup_file_if_necessary( std::string const& filename ) ;
	std::map< std::string, std::string > const& backed_up_files() const { return m_backed_up_files ; }	
protected:
	virtual std::string backup_file_if_necessary_impl( std::string const& filename ) = 0 ;
private:
	std::map< std::string, std::string > m_backed_up_files ;	
} ;

class ToTmpFileBackupCreator: public FileBackupCreator
{
private:
	std::string backup_file_if_necessary_impl( std::string const& filename ) ;
	bool backup_file_if_necessary( std::string const& filename, std::string const& backup_filename ) ;
} ;

class ToNumberedFileBackupCreator: public FileBackupCreator
{
public:
	ToNumberedFileBackupCreator( std::size_t max_number_of_backups = 10 ): m_max_number_of_backups( max_number_of_backups ) {}
	
private:
	std::string backup_file_if_necessary_impl( std::string const& filename ) ;
	std::size_t m_max_number_of_backups ;
} ;

#endif
