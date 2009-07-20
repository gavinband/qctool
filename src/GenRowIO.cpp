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
	return result ;
}


std::istream& GenRowIdentifyingData::read_from_text_stream( std::istream& aStream ) {
	aStream >> m_SNPID ;
	aStream >> m_RSID ;
	aStream >> m_SNP_position ;
	aStream >> m_1st_allele ;
	return aStream >> m_2nd_allele ;
}	

std::istream& GenRow::read_from_text_stream( std::istream& inStream ) {
	// For speed, first read a line from the file.
	std::string line ;
	std::getline( inStream, line ) ;
	std::istringstream aStream( line ) ;

	GenRowIdentifyingData::read_from_text_stream( aStream ) ;

	if( aStream ) {
		set_number_of_samples( 0 ) ;	
		int count = 0 ;
	
		while( aStream ) {
			GenotypeProportions sample_genotype_proportions( 0.0, 0.0, 0.0 ) ;
			sample_genotype_proportions.AA() = read_float( aStream ) ;
			if( !aStream ) break ;
			++count ;
			sample_genotype_proportions.AB() = read_float( aStream ) ;
			if( !aStream ) break ;
			++count ;
			sample_genotype_proportions.BB() = read_float( aStream ) ;
			if( !aStream ) break ;
			++count ;
			add_genotype_proportions( sample_genotype_proportions ) ;
			aStream.peek() ; // flag eof if we reached it.
		} ;

		// At the end of the row, after any whitespace, we allow either a newline or eof.
		if( count % 3 != 0 ) {
			throw BadRowFormatException( "Expected line to contain a full set of probabilities (AA, AB and BB)." ) ;
		}
	}
	else {
		inStream.setstate( std::ios::failbit ) ;
	}

	return inStream ;
}

std::ostream& GenRow::write_to_text_stream( std::ostream& aStream ) const {
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
			<< get_AA_probability(i) << " "
			<< get_AB_probability(i) << " "
			<< get_BB_probability(i) ;

	}

	aStream << "\n" ;
	
	return aStream ;
}

std::istream& operator>>( std::istream& inStream, GenRow& aRow ) {
	aRow.read_from_text_stream( inStream ) ;
	return inStream ;
}

std::ostream& operator<<( std::ostream& aStream, GenRow const& aRow ) {
	aRow.write_to_text_stream( aStream ) ;
	return aStream ;
}


