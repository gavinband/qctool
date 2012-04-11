
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_FILEUTIL_HPP__
#define QCTOOL_FILEUTIL_HPP__

#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <iterator>
#include <iostream>
#include <exception>
#include <memory>

typedef std::auto_ptr< std::istream > INPUT_FILE_PTR ;
typedef std::auto_ptr< std::ostream > OUTPUT_FILE_PTR ;

struct FileError: public std::exception
{
public:
	FileError( std::string const& filename ): m_filename( filename ) {}
	virtual ~FileError() throw() {}
	
	char const* what() const throw() { return "FileError" ; }
	std::string const& filename() const { return m_filename ; }

private:
	std::string const m_filename ;
} ;

// thrown to indicate that an output file with wildcard appeared, but the corresponding input
// file had no wildcard.
struct FileNotOpenedError: public FileError
{
	FileNotOpenedError( std::string const& filename ): FileError( filename ) {}
	char const* what() const throw() { return "FileNotOpenedError" ; }
} ;     

struct InvalidFileFormatError: public FileError
{
	InvalidFileFormatError( std::string const& filename ): FileError( filename ) {}
	char const* what() const throw() { return "InvalidFileFormatError" ; }
} ;

//
enum FileCompressionType {e_None = 0x1, e_Gzip = 0x2, e_Bzip2 = 0x4, e_FileCompressionMask = 0xf } ;
enum FileModeType { e_TextMode = 0x10, e_BinaryMode = 0x20, e_FileModeMask = 0xf0 } ;

FileCompressionType determine_file_compression( std::string const& filename ) ;
FileModeType determine_file_mode( std::string const& filename ) ;

// Return a stream representing a given input file, optionally with decompression
INPUT_FILE_PTR
open_file_for_input( std::string const& filename, int mode_flags ) ;

// Return a stream representing a given input file, attempting to auto-detect any compression used.
INPUT_FILE_PTR
open_file_for_input( std::string const& filename ) ;

// Return a stream representing a given output file, optionally with compression.
OUTPUT_FILE_PTR
open_file_for_output( std::string const& filename, int mode_flags ) ;

// Return a stream representing a given output file, guessing compression option from the filename.
OUTPUT_FILE_PTR
open_file_for_output( std::string const& filename ) ;

bool exists( std::string const& filename ) ;
bool is_regular( std::string const& filename ) ;
void copy_file( std::string const& filename1, std::string const& filename2 ) ;
void rename( std::string const& filename1, std::string const& filename2 ) ;
void remove( std::string const& filename ) ;
#endif

