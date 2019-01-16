
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_SHOTGUNSTOCHASTICSEARCH_HPP
#define METRO_SHOTGUNSTOCHASTICSEARCH_HPP

#include <vector>
#include <unordered_map>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>

namespace metro {
	/*
	* This class performs a FINEMAP-like 'shotgun stochastic search'.
	* Given an integer N, starting with the empty state, this class searches by generating
	* all possible changes to the current state by adding / changing / removing any of the
	* N possible indicators.  It then computes a likelihood using the user-supplied
	* likelihood function.  Finally, it updates the current state by moving to one of the new states
	* uniformly weighted by the likelihood.
	* See Benner et al, "FINEMAP: Efficient variable selection using summary data from genome-wide association studies",
	* Bioinformatics (2016) for algorithm details.
	* This class merely implements the 'search' aspect of this algorithm.  To make efficient use of it,
	* the user must use a likelihood function that e.g. prohibits too-large numbers of predictor variables.
	*/
	struct ShotgunStochasticSearch {
		typedef std::vector< std::size_t > SelectedStates ;
		typedef boost::function< double ( SelectedStates const& ) > ComputeLL ;
		struct StateHash {
		public:
			size_t operator()(const SelectedStates &p) const {
				return boost::hash_value(p);
			}
		} ;

		typedef std::unordered_map< SelectedStates, double, StateHash > Store ;

		ShotgunStochasticSearch(
			std::size_t n,
			ComputeLL compute_ll,
			std::uint32_t rng_seed
		) ;

		/* Compute the next search update with the given likelihood function */
		SelectedStates const& update() ;

		/* Return our store of states visited so far */ 
		Store const& visited_states() const ;
		void visited_states( boost::function< void( SelectedStates const&, double ) > callback ) const ;
	
	private:
		std::size_t const m_N ;
		ShotgunStochasticSearch::ComputeLL m_compute_ll ;
		SelectedStates m_current_state ;
		Store m_lls ;
		boost::random::mt19937 m_rng;
	
	private:
		/* Generate all states obtained by deleting a member
		of the current state, and append to the given vector */
		void generate_deletions( std::vector< SelectedStates >* result ) ;

		/* Generate all states obtained by adding a member or changing a member
		of the current state, and append to the given vector */
		void generate_additions_and_changes( std::vector< SelectedStates >* result ) ;
	} ;
}

#endif


