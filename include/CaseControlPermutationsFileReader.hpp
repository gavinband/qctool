
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef CASE_CONTROL_PERMUTATIONS_FILE_READER_HPP
#define CASE_CONTROL_PERMUTATIONS_FILE_READER_HPP

#include <exception>
#include <string>
#include <fstream>
#include <vector>
#include "endianness_utils.hpp"

struct PermutationFileFormatError: public std::exception
{
	char const* what() const throw() { return "PermutationFileFormatError" ; }
} ;

struct PermutationFileEmptyRowError: public std::exception
{
	char const* what() const throw() { return "PermutationFileEmptyRowError" ; }
} ;

struct PermutationFileEmptyError: public std::exception
{
	char const* what() const throw() { return "PermutationFileEmptyError" ; }
} ;

struct PermutationFileNumbersOfZeroesDifferError: public std::exception
{
	char const* what() const throw() { return "PermutationFileNumbersOfZeroesDifferError" ; }
} ;

struct PermutationFileMalformedError: public std::exception
{
	char const* what() const throw() { return "PermutationFileMalformedError" ; }
} ;

struct PermutationFileFirstPermutationMalformedError: public std::exception
{
	char const* what() const throw() { return "PermutationFileFirstPermutationMalformedError" ; }
} ;

struct CaseControlPermutationsFileReader 
{
	 CaseControlPermutationsFileReader( std::string const& permutations_filename )
	: m_permutations( read_permutations_file( permutations_filename ))
	{
		check_numbers_of_zeroes() ;
	}
	
	std::size_t number_of_zeroes() const { return count_number_of_zeroes_in_permutation( 0 ) ; }
	std::vector< std::vector< char > > permutations() const { return m_permutations ; }
	
private:

	std::vector< std::vector< char > > read_permutations_file( std::string const& permutations_filename ) {
		std::ifstream file( permutations_filename.c_str(), std::ios::binary ) ;
		std::string line ;
		std::vector< std::vector< char > > permutations ;

		uint32_t offset, number_of_zeroes, number_of_ones, number_of_permutations ;
		
		read_little_endian_integer( file, &offset ) ;
		file.ignore( offset ) ;
		read_little_endian_integer( file, &number_of_zeroes ) ;
		read_little_endian_integer( file, &number_of_ones ) ;
		read_little_endian_integer( file, &number_of_permutations ) ;

		std::size_t size_of_permutations = number_of_zeroes + number_of_ones ;

		for( uint32_t i = 0; i < number_of_permutations; ++i ) {
			permutations.push_back( std::vector< char >( size_of_permutations )) ;
			file.read( &(permutations.back()[0]), size_of_permutations ) ;
			if( !file || (static_cast< std::size_t > (file.gcount()) != size_of_permutations )) {
				throw PermutationFileMalformedError() ;
			}
		}
		
		if( permutations.empty() ) {
			throw PermutationFileEmptyError() ;
		}
		
		return permutations ;
	}
	
	void check_numbers_of_zeroes() const {
		std::size_t number_of_zeroes = count_number_of_zeroes_in_permutation( 0 ) ;
		for( std::size_t i = 1; i < m_permutations.size(); ++i ) {
			if( count_number_of_zeroes_in_permutation( i ) != number_of_zeroes ) {
				throw PermutationFileNumbersOfZeroesDifferError() ;
			}
		}		
	}
	
	std::size_t count_number_of_zeroes_in_permutation( std::size_t index ) const {
		std::size_t number_of_zeroes = 0 ;
		for( std::size_t i = 0; i < m_permutations[ index ].size(); ++i ) {
			if( m_permutations[ index ][ i ] == 0 ) {
				++number_of_zeroes ;
			}
		}
		return number_of_zeroes ;
	}
	
private:
	
	std::vector< std::vector< char > > m_permutations ;
} ;


#endif
