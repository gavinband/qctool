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

std::istream& operator>>( std::istream& inStream, GenRow& aRow ) {
	// For speed, first read a line from the file.
	std::string line ;
	std::getline( inStream, line ) ;
	std::istringstream aStream( line ) ;
	Whitespace whitespace ;

	aStream >> aRow.m_SNPID ;
	aStream >> aRow.m_RSID ;
	aStream >> aRow.m_SNP_position ;
	aStream >> aRow.m_1st_allele ;
	aStream >> aRow.m_2nd_allele ;

	aRow.m_genotype_proportions.clear() ;	

	do {
		GenotypeProportions sample_allele_proportions( 0.0, 0.0, 0.0 ) ;

		sample_allele_proportions.proportion_of_AA() = read_float( aStream ) ;
		sample_allele_proportions.proportion_of_AB() = read_float( aStream ) ;
		sample_allele_proportions.proportion_of_BB() = read_float( aStream ) ;
		aRow.m_genotype_proportions.push_back( sample_allele_proportions ) ;
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
	return aStream ;
}


std::ostream& operator<<( std::ostream& aStream, GenRow const& aRow ) {
	aStream
		<< aRow.SNPID() << " "
		<< aRow.RSID() << " "
		<< aRow.SNP_position() << " "
		<< aRow.first_allele() << " "
		<< aRow.second_allele() << " " ;

	for( std::size_t i = 0 ; i < aRow.number_of_samples() ; ++i ) {
		if( i > 0 ) {
			aStream << " " ;
		}
		
		aStream
			<< aRow.genotype_proportions_for_sample(i).proportion_of_AA() << " "
			<< aRow.genotype_proportions_for_sample(i).proportion_of_AB() << " "
			<< aRow.genotype_proportions_for_sample(i).proportion_of_BB() ;

	}

	return aStream << std::endl ;
}


