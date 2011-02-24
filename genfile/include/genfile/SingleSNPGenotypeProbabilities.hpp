#ifndef GENFILE_GENOTYPE_PROBABILITIES_HPP
#define GENFILE_GENOTYPE_PROBABILITIES_HPP

#include <cassert>
#include <vector>
#include "genfile/get_set.hpp"

namespace genfile {
	
	struct SingleSNPGenotypeProbabilities
	// class SingleSNPGenotypeProbabilities
	// This class holds genotype probabilities for a number of individuals at a single SNP,
	// and provides methods to access that data.
	// Invariant: the sum of probabilities for each individual is at most 1, to 3d.p accuracy.
	{
	public:
		// zero-initialised set of probabilities.
		SingleSNPGenotypeProbabilities( std::size_t number_of_samples = 0 ) ;
		// Construct from a range of doubles.
		// These come in the order: AA(sample 1) AB(sample 1) BB( sample 1) AA( sample 2 )...
		SingleSNPGenotypeProbabilities( double const* begin, double const* end ) ;
		// Copy-construct
		SingleSNPGenotypeProbabilities( SingleSNPGenotypeProbabilities const& other ) ;
		// Assign
		SingleSNPGenotypeProbabilities& operator=( SingleSNPGenotypeProbabilities const& other ) ;

		// Return the number of samples we have genotypes for.
		std::size_t get_number_of_samples() const {
			return m_number_of_samples ;
		}
		// Return the number of samples we have genotypes for.
		std::size_t size() const {
			return m_number_of_samples ;
		}

		// Return the number of bytes used
		std::size_t get_memory_usage_in_bytes() const ;
		
		// Resize to give storage for the given number of samples.
		// This has the same semantics as vector resize, i.e. shortenings will lose data,
		// largenings will keep data in the existing sample.
		void resize( std::size_t number_of_samples ) ;

		// Set all the probabilities.
		void set( double const* begin, double const* end ) ;
		// Set the probabilities for individual i.
		void set( std::size_t i, double AA, double AB, double BB ) ;

		double AA( std::size_t i ) const {
			return operator()( i, 0 ) ;
		}
		double AB( std::size_t i ) const {
			return operator()( i, 1 ) ;
		}
		double BB( std::size_t i ) const {
			return operator()( i, 2 ) ;
		}
		// Return genotype probability g for individual i.
		// g is 0:AA 1:AB 2:BB
		double operator()( std::size_t i, std::size_t g ) const {
			assert( i < m_number_of_samples ) ;
			assert( g < 3 ) ;
			return m_probabilities[ 3*i + g ] ;
		}
		// Return the sum of AA, AB, BB for sample i
		double sum( std::size_t i ) const ;
		// Return the amount of null call for sample i, i.e. 1 - sum( i )
		double null_call( std::size_t i ) const ; 

	private:
		std::size_t m_number_of_samples ;
		std::vector< double > m_probabilities ;
		void check_invariant( std::vector< double > const& probabilities ) ;
		void check_invariant( std::size_t i, double AA, double AB, double BB ) ;
	} ;
	
	// setter for SingleSNPGenotypeProbabilities
	template<>
	struct GenotypeSetter< SingleSNPGenotypeProbabilities >
	{
		GenotypeSetter( SingleSNPGenotypeProbabilities& genotypes ) ;
		void operator()( std::size_t i, double, double, double ) const ;	

		private:
			SingleSNPGenotypeProbabilities& m_genotypes ;
	} ;
}

#endif
