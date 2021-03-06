
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <string>
#include <cstdio>
#include "FileBackupCreator.hpp"
#include "genfile/string_utils/string_utils.hpp"
#include <boost/filesystem.hpp>

namespace BFS = boost::filesystem ;

void FileBackupCreator::backup_file_if_necessary( std::string const& filename ) {
	std::string backup_filename = backup_file_if_necessary_impl( filename ) ;
	if( backup_filename != "" ) {
		m_backed_up_files[ filename ] = backup_filename ;
	}
}

std::string ToTmpFileBackupCreator::backup_file_if_necessary_impl( std::string const& filename ) {
	std::string backup_filename ;
	backup_filename = tmpnam(0) ;
	if( backup_file_if_necessary( filename, backup_filename )) {
		return backup_filename ;
	}
	else {
		return "" ;
	}
}

bool ToTmpFileBackupCreator::backup_file_if_necessary( std::string const& filename, std::string const& backup_filename ) {
	if( filename != "" && BFS::exists( filename )) {
		BFS::rename( filename, backup_filename ) ;
		return true ;
	}
	return false ;
}

std::string ToNumberedFileBackupCreator::backup_file_if_necessary_impl( std::string const& filename ) {
	std::vector< std::string > filename_stack ;
	filename_stack.push_back( filename ) ;

	// Generate a list of numbered filenames whose numbers we'll increment.
	// If we hit the limit set by m_max_number_of_backups, we delete the last of these.
	while( BFS::exists( filename_stack.back() )) {
		std::string next_backup = filename + ".~" + genfile::string_utils::to_string( filename_stack.size() ) ;
		filename_stack.push_back( next_backup ) ;
		if( filename_stack.size() > m_max_number_of_backups ) {
			if( BFS::exists( next_backup )) {
				BFS::remove( next_backup ) ;
			}
			break ;
		}
	}
	
	// Now increment the numbers
	for( std::size_t i = filename_stack.size(); i > 1; --i ) {
		BFS::rename( filename_stack[ i - 2 ], filename_stack[ i - 1 ] ) ;
	}
	
	return (filename_stack.size() > 1) ? filename_stack[1] : "" ;
}
