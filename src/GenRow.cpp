#include <sstream>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>
#include <cmath>

#include "../config.hpp"
#include "GToolException.hpp"
#include "GenRow.hpp"
#include "Whitespace.hpp"
#include "GenotypeProportions.hpp"

std::size_t GenRow::number_of_columns() const {
	return (m_genotype_proportions.size() * 3) + std::size_t(5) ;
}


std::size_t GenRow::number_of_samples() const {
	return m_genotype_proportions.size() ;
}

GenotypeProportions const& GenRow::genotype_proportions_for_sample( std::size_t sampleNumber ) const {
	assert( sampleNumber < number_of_samples() ) ;
	return m_genotype_proportions[ sampleNumber ] ;
}

void GenRow::reserveSpaceForNSamples( std::size_t number_of_samples ) {
	m_genotype_proportions.reserve( number_of_samples ) ;
}

bool GenRow::operator==( GenRow const& right ) const {
	return ( m_SNPID == right.m_SNPID )
		&& ( m_RSID == right.m_RSID )
		&& ( m_SNP_position == right.m_SNP_position )
		&& ( m_1st_allele == right.m_1st_allele )
		&& ( m_2nd_allele == right.m_2nd_allele )
		&& ( m_genotype_proportions == right.m_genotype_proportions ) ;

}

double read_float( std::istream& aStream ) {
#if GTOOL_USE_FAST_FLOAT_PARSER
	std::string float_str ;
	aStream >> float_str ;
	std::string::const_iterator i = float_str.begin() ,
		end_i = float_str.end() ;

	unsigned long integer_part = 0.0, fractional_part = 0ul;
	unsigned long* this_part = &integer_part ;
	int fractional_part_length = 0 ;
	bool good = true ;

	for( ; good && i != end_i; ++i ) {
		if( *i == '.' ) {
			this_part = &fractional_part ;
			fractional_part_length = std::distance( i, end_i ) - 1 ;
		}
		else if(( '0' <= *i ) && ( '9' >= *i ) ) {
			(*this_part) *= 10 ;
			(*this_part) += (*i - '0') ;
		}
		else {
			good = false ;
		}
	}

	double result ;

	if( good ) {
		result = static_cast< double >( integer_part ) + 
			+ ( static_cast< double >( fractional_part ) / std::pow( 10.0, fractional_part_length )) ;
	}
	else {
		// Failed to parse.  Try to parse using usual method.
		std::istringstream inStream( float_str ) ;
		inStream >> result ;
	}
#else
	float result ;
	aStream >> result ;
#endif
	aStream.peek() ; // flag eof
	return result ;
}


void GenRow::read_from_text_stream( std::istream& inStream ) {
	// For speed, first read a line from the file.
	std::string line ;
	std::getline( inStream, line ) ;
	std::istringstream aStream( line ) ;
	Whitespace whitespace ;

	aStream >> m_SNPID ;
	aStream >> m_RSID ;
	aStream >> m_SNP_position ;
	aStream >> m_1st_allele ;
	aStream >> m_2nd_allele ;

	m_genotype_proportions.clear() ;	

	do {
		GenotypeProportions sample_allele_proportions( 0.0, 0.0, 0.0 ) ;

		sample_allele_proportions.proportion_of_AA() = read_float( aStream ) ;
		sample_allele_proportions.proportion_of_AB() = read_float( aStream ) ;
		sample_allele_proportions.proportion_of_BB() = read_float( aStream ) ;
		m_genotype_proportions.push_back( sample_allele_proportions ) ;
		aStream >> whitespace ;
	}
	while( !aStream.eof() && static_cast<char>( aStream.peek()) != '\n' ) ;

	// At the end of the row, after any whitespace, we allow either a newline or eof.
	if( !aStream.eof() ) {
		if( static_cast<char>( aStream.peek() ) == '\n' ) {
			aStream.get() ;
			aStream.peek() ; // ensure eof is flagged if we've reached the end of the file.
		}
		else {
			throw BadRowFormatException( "Expected newline or eof at end of row." ) ;
		}
	}

	inStream.peek() ; // flag eof
	
}

void GenRow::write_to_text_stream( std::ostream& aStream ) const {
	aStream
		<< SNPID() << " "
		<< RSID() << " "
		<< SNP_position() << " "
		<< first_allele() << " "
		<< second_allele() << " " ;

	for( std::size_t i = 0 ; i < number_of_samples() ; ++i ) {
		if( i > 0 ) {
			aStream << " " ;
		}
		
		aStream
			<< genotype_proportions_for_sample(i).AA() << " "
			<< genotype_proportions_for_sample(i).AB() << " "
			<< genotype_proportions_for_sample(i).BB() ;

	}

	aStream << "\n" ;
}

void GenRow::read_from_binary_stream( std::istream& aStream ) {
	std::size_t a_number_of_samples ;
	// read data into local variables so not to destroy this object's state.
	std::string SNPID ;
	std::string RSID ;
	int SNP_position ;
	char first_allele, second_allele ;
	std::vector< GenotypeProportions > genotype_proportions ;

	aStream >> a_number_of_samples ;
	aStream >> SNPID ;
	aStream >> RSID ;
	aStream >> SNP_position ;
	aStream >> first_allele ;
	aStream >> second_allele ;

	char space_char ;
	aStream.read( &space_char, 1 ) ;
	if( space_char != ' ' ) {
		throw BadRowFormatException( "Expected space before genotype proportions." ) ;
	}
	genotype_proportions.resize( a_number_of_samples ) ;

	for( std::size_t i = 0; i < a_number_of_samples; ++i ) {
		double AA, AB, BB ;
		aStream.read( reinterpret_cast< char* >(&AA), sizeof(double)) ;
		aStream.read( reinterpret_cast< char* >(&AB), sizeof(double)) ;
		aStream.read( reinterpret_cast< char* >(&BB), sizeof(double)) ;
		genotype_proportions[i] = GenotypeProportions( AA, AB, BB ) ;
	}
	
	m_SNPID = SNPID ;
	m_RSID = RSID ;
	m_SNP_position = SNP_position ;
	m_1st_allele = first_allele ;
	m_2nd_allele = second_allele ;
	m_genotype_proportions = genotype_proportions ;
	
	aStream.peek() ;	
}

void GenRow::write_to_binary_stream( std::ostream& aStream ) const {
	aStream
		<< number_of_samples() << " "
		<< SNPID() << " "
		<< RSID() << " "
		<< SNP_position() << " "
		<< first_allele() << " "
		<< second_allele() << " " ;

	for( std::size_t i = 0 ; i < number_of_samples() ; ++i ) {
		double
			AA = genotype_proportions_for_sample(i).AA(),
			AB = genotype_proportions_for_sample(i).AB(),
			BB = genotype_proportions_for_sample(i).BB() ;
		
		aStream.write( reinterpret_cast< char* >(&AA), sizeof( double )) ;
		aStream.write( reinterpret_cast< char* >(&AB), sizeof( double )) ;
		aStream.write( reinterpret_cast< char* >(&BB), sizeof( double )) ;
	}
}


std::istream& operator>>( std::istream& inStream, GenRow& aRow ) {
	aRow.read_from_text_stream( inStream ) ;
	return inStream ;
}

std::ostream& operator<<( std::ostream& aStream, GenRow const& aRow ) {
	aRow.write_to_text_stream( aStream ) ;
	return aStream ;
}


