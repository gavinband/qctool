#include <iostream>
#include <vector>
#include <numeric>
#include "genfile/Error.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "genfile/string_utils.hpp"

namespace genfile {
	
	SingleSNPGenotypeProbabilities::SingleSNPGenotypeProbabilities():
	 	m_number_of_samples( 0 )
	{}

	// Construct from a range of doubles.
	// These come in the order: AA(sample 1) AB(sample 1) BB( sample 1) AA( sample 2 )...
	SingleSNPGenotypeProbabilities::SingleSNPGenotypeProbabilities( double const* begin, double const* end ):
		m_probabilities( begin, end ),
		m_number_of_samples( m_probabilities.size() % 3 )
	{
		check_invariant( m_probabilities ) ;
	}

	void SingleSNPGenotypeProbabilities::check_invariant( std::vector< double > const& probabilities ) {
		if( probabilities.size() % 3 != 0 ) {
			throw BadArgumentError( "SingleSNPGenotypeProbabilities::check_invariant()", "size of probabilities not a multiple of 3" ) ;
		}
		std::size_t const N = probabilities.size() / 3 ;
		for( std::size_t i = 0; i < N; ++i ) {
			check_invariant( i, probabilities[3*i], probabilities[3*i+1], probabilities[3*i+2] ) ;
		}
	}

	void SingleSNPGenotypeProbabilities::check_invariant( std::size_t i, double AA, double AB, double BB ) {
		// BGen format is only accurate to 4 decimal places.  The sum of three numbers rounded to 4d.p is
		// only accurate to 3 d.p.  So we had better only check 3dp of accuracy here.
		if( AA + AB + BB >= 1.0005 ) {
			// std::cerr << "AA: " << AA << " AB:" << AB << " BB:" << BB << " sum:" << ( AA + AB + BB ) << ".\n" ;
			throw BadArgumentError( "SingleSNPGenotypeProbabilities::check_invariant()", "sum of probabilities (individual " + string_utils::to_string( i ) + ") greater than 1" ) ;
		}
	}

	SingleSNPGenotypeProbabilities::SingleSNPGenotypeProbabilities( SingleSNPGenotypeProbabilities const& other ):
		m_probabilities( other.m_probabilities ),
		m_number_of_samples( other.m_number_of_samples )
	{}
	
	void SingleSNPGenotypeProbabilities::resize( std::size_t number_of_samples ) {
		m_probabilities.resize( number_of_samples * 3, 0.0 ) ;
		m_number_of_samples = number_of_samples ;
	}
	
	// Assign
	SingleSNPGenotypeProbabilities& SingleSNPGenotypeProbabilities::operator=( SingleSNPGenotypeProbabilities const& other ) {
		m_probabilities = other.m_probabilities ;
		m_number_of_samples = other.m_number_of_samples ;
		return *this ;
	}

	std::size_t SingleSNPGenotypeProbabilities::get_memory_usage_in_bytes() const {
		return sizeof( SingleSNPGenotypeProbabilities ) + m_probabilities.capacity() * sizeof( double ) ;
	}

	// Set all the probabilities.
	void SingleSNPGenotypeProbabilities::set( double const* begin, double const* end ) {
		std::vector< double > probabilities( begin, end ) ;
		check_invariant( probabilities ) ;
		m_probabilities.swap( probabilities ) ;
		m_number_of_samples = m_probabilities.size() % 3 ;
	}

	// Set the probabilities for individual i.
	void SingleSNPGenotypeProbabilities::set( std::size_t i, double AA, double AB, double BB ) {
		assert( i < m_number_of_samples ) ;
		check_invariant( i, AA, AB, BB ) ;
		m_probabilities[ 3*i+0 ] = AA ;
		m_probabilities[ 3*i+1 ] = AB ;
		m_probabilities[ 3*i+2 ] = BB ;
	}

	// Return the sum of AA, AB, BB for sample i
	double SingleSNPGenotypeProbabilities::sum( std::size_t i ) const {
		assert( i < m_number_of_samples ) ;
		return std::accumulate( &m_probabilities[3*i], &m_probabilities[3*i+3], double( 0.0 ) ) ;
	}

	// Return the amount of null call for sample i, i.e. 1 - sum( i )
	double SingleSNPGenotypeProbabilities::null_call( std::size_t i ) const {
		return 1.0 - sum( i ) ;
	}
	
	GenotypeSetter< SingleSNPGenotypeProbabilities >::GenotypeSetter( SingleSNPGenotypeProbabilities& genotypes ):
		m_genotypes( genotypes )
	{}

	void GenotypeSetter< SingleSNPGenotypeProbabilities >::operator()( std::size_t i, double AA, double AB, double BB ) const {
		m_genotypes.set( i, AA, AB, BB ) ;
	}
}

