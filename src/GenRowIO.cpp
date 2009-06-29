#include <string>
#include <sstream>
#include <iostream>
#include <cassert>
#include <cmath>
#include <boost/bind.hpp>

#include "../config.hpp"
#include "GToolException.hpp"
#include "GenRow.hpp"
#include "Whitespace.hpp"
#include "GenotypeProportions.hpp"
#include "genbin.hpp"


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

// The following two functions implement the genbin specification described here:
// http://www.well.ox.ac.uk/~gav/binary_file_format.html

void GenRow::read_from_binary_stream( std::istream& aStream ) {
	genbin::read_snp_block(
		aStream,
		boost::bind< void >( &GenRow::set_number_of_samples, this, _1 ),
		boost::bind< void >( &GenRow::set_SNPID, this, _1 ),
		boost::bind< void >( &GenRow::set_RSID, this, _1 ),
		boost::bind< void >( &GenRow::set_SNP_position, this, _1 ),
		boost::bind< void >( &GenRow::set_alleles, this, _1, _2 ),
		boost::bind< void >( &GenRow::set_genotype_probabilities, this, _1, _2, _3, _4 )
	) ;
}

void GenRow::write_to_binary_stream( std::ostream& aStream ) const {
	genbin::write_snp_block(
		aStream,
		m_genotype_proportions.size(),
		m_SNPID,
		m_RSID,
		m_SNP_position,
		m_1st_allele,
		m_2nd_allele,
		boost::bind< double >( &GenRow::get_AA_probability, this, _1 ),
		boost::bind< double >( &GenRow::get_AB_probability, this, _1 ),
		boost::bind< double >( &GenRow::get_BB_probability, this, _1 )
	) ;		
}

std::istream& operator>>( std::istream& inStream, GenRow& aRow ) {
	aRow.read_from_text_stream( inStream ) ;
	return inStream ;
}

std::ostream& operator<<( std::ostream& aStream, GenRow const& aRow ) {
	aRow.write_to_text_stream( aStream ) ;
	return aStream ;
}


