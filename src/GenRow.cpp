#include <sstream>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>

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


std::istream& operator>>( std::istream& aStream, GenRow& aRow ) {
	// prevent whitespace skipping.
	
	aStream >> std::noskipws ;
	Whitespace whitespace( " \t" ) ;

	aStream >> aRow.m_SNPID >> whitespace ;
	if( whitespace.empty() ) throw BadRowFormatException( "Expected space after SNPID" ) ;
	aStream >> aRow.m_RSID >> whitespace ;
	if( whitespace.empty() ) throw BadRowFormatException( "Expected space after RSID" ) ;
	aStream >> aRow.m_SNP_position >> whitespace ;
	if( whitespace.empty() ) throw BadRowFormatException( "Expected space after SNP position" ) ;
	aStream >> aRow.m_1st_allele >> whitespace ;
	if( whitespace.empty() ) throw BadRowFormatException( "Expected space after first allele" ) ;
	aStream >> aRow.m_2nd_allele >> whitespace ;
	if( whitespace.empty() ) throw BadRowFormatException( "Expected space after second allele" ) ;

	aRow.m_genotype_proportions.clear() ;	

	do {
		GenotypeProportions sample_allele_proportions( 0.0, 0.0, 0.0 ) ;

		aStream >> sample_allele_proportions.proportion_of_AA() >> whitespace ;
		if( whitespace.empty() ) throw BadRowFormatException( "Expected space after first allele proportion" ) ;
		aStream >> sample_allele_proportions.proportion_of_AB() >> whitespace ;
		if( whitespace.empty() ) throw BadRowFormatException( "Expected space after second allele proportion" ) ;
		aStream >> sample_allele_proportions.proportion_of_BB() >> whitespace ;

		aRow.m_genotype_proportions.push_back( sample_allele_proportions ) ;
	}
	while( !whitespace.empty() && !aStream.eof() && static_cast<char>( aStream.peek()) != '\n' ) ;

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

	// restore whitespace skipping behaviour
	return aStream >> std::skipws ;
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


