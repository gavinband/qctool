#ifndef CASE_CONTROL_PERMUTATIONS_FILE_READER_HPP
#define CASE_CONTROL_PERMUTATIONS_FILE_READER_HPP

#include <exception>
#include <string>
#include <fstream>
#include <vector>

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

struct PermutationFileRowSizesDifferError: public std::exception
{
	char const* what() const throw() { return "PermutationFileRowSizesDifferError" ; }
} ;

struct PermutationFileNumbersOfZeroesDifferError: public std::exception
{
	char const* what() const throw() { return "PermutationFileNumbersOfZeroesDifferError" ; }
} ;

struct PermutationFileMalformedEntryError: public std::exception
{
	char const* what() const throw() { return "PermutationFileMalformedEntryError" ; }
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
		std::ifstream file( permutations_filename.c_str() ) ;
		std::string line ;
		std::vector< std::vector< char > > permutations ;

		std::size_t size_of_permutations = 0u ;

		while( std::getline( file, line )) {
			if( line.substr( 0, 2 ) == "//" ) {
				// skip the comment line
			}
			else {
				// line is a whitespace-separated list of 0s and 1s
				std::vector< std::string > entries = split_and_strip( line, " " ) ;
				if( entries.size() == 0 ) {
					throw PermutationFileEmptyRowError() ;
				}
				else if( size_of_permutations == 0 ) {
					size_of_permutations = entries.size() ;
				}
				else if( size_of_permutations != entries.size() ){
					throw PermutationFileRowSizesDifferError() ;
				}

				permutations.push_back( std::vector< char >( size_of_permutations )) ;
				for( std::size_t i = 0; i < size_of_permutations ; ++i ) {
					if( entries[i].size() == 1 ) {
						permutations.back()[i] = entries[i][0] ;
					}
					else {
						throw PermutationFileMalformedEntryError() ;
					}	
				}
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
