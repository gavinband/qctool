
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include "metro/ShotgunStochasticSearch.hpp"

// #define DEBUG 1

namespace metro {
	ShotgunStochasticSearch::ShotgunStochasticSearch(
		std::size_t n,
		ShotgunStochasticSearch::ComputeLL compute_ll,
		std::uint32_t rng_seed
	):
		m_N(n),
		m_compute_ll( compute_ll ),
		m_rng( rng_seed )
	{
		assert( m_compute_ll ) ;
		m_lls.insert( std::make_pair( m_current_state, m_compute_ll( m_current_state ))) ;
	}
	
	ShotgunStochasticSearch::SelectedStates const& ShotgunStochasticSearch::update() {
		// sample states...
		std::size_t count = 0 ;
	
		std::vector< SelectedStates > new_states ;
#if DEBUG
		std::cerr << "Current state: " << m_current_state.size() << "elts: " << m_current_state << ".\n" ;
#endif
		generate_deletions( &new_states ) ;
		generate_additions_and_changes( &new_states ) ;
#if DEBUG
		std::cerr << "Generated " << new_states.size() << "new states:\n" ;
		for( std::size_t i = 0; i < new_states.size(); ++i ) {
			std::cerr << " -- " << new_states[i] << "\n" ;
		}
#endif

		std::vector< double > lls( new_states.size(), 0 ) ;
		for( std::size_t i = 0; i < new_states.size(); ++i ) {
			Store::const_iterator where = m_lls.find( new_states[i] ) ;
			if( where == m_lls.end() ) {
				lls[i] = m_compute_ll( new_states[i] ) ;
				m_lls.insert( std::make_pair( new_states[i], lls[i] ) ) ;
			} else {
				lls[i] = where->second ;
			}
		}
	
		double const total_ll = std::accumulate( lls.begin(), lls.end(), 0.0 ) ;
		// Pick a new state at random
		boost::random::discrete_distribution< std::size_t, double > dist( lls ) ;
		m_current_state = new_states[ dist(m_rng) ] ;
		
		return m_current_state ;
	}

	ShotgunStochasticSearch::Store const& ShotgunStochasticSearch::visited_states() const {
		return m_lls ;
	}
	
	void ShotgunStochasticSearch::generate_deletions(
		std::vector< ShotgunStochasticSearch::SelectedStates >* result
	) {
		SelectedStates new_state ;
		// deletions...
		for( std::size_t i = 0; i < m_current_state.size(); ++i ) {
			new_state.resize( m_current_state.size() - 1 ) ;
			std::copy( m_current_state.begin(), m_current_state.begin() + i, new_state.begin() ) ;
			std::copy( m_current_state.begin() + i + 1, m_current_state.end(), new_state.begin() + i ) ;
			result->push_back( new_state ) ;
		}
	}

	void ShotgunStochasticSearch::generate_additions_and_changes(
		std::vector< ShotgunStochasticSearch::SelectedStates >* result
	) {
		SelectedStates new_addition_state( m_current_state.size() + 1 ) ;
		SelectedStates new_change_state( m_current_state.size() ) ;
		std::size_t begin_addition = 0 ;
		std::size_t end_addition = 0 ;
	
		SelectedStates::const_iterator const begin_i = m_current_state.begin() ;
		SelectedStates::const_iterator const end_i = m_current_state.end() ;
		for( std::size_t i = 0; i <= m_current_state.size(); ++i, begin_addition = end_addition + 1 ) {
			end_addition = (i == m_current_state.size()) ? m_N : m_current_state[i] ;
			for( std::size_t a = begin_addition; a < end_addition; ++a ) {
				// In this loop, a will visit everything that can be added
				// and i will point to the element of current state that follows a.

				// Build the addition states
				std::copy( begin_i, begin_i + i, new_addition_state.begin() ) ;
				new_addition_state[i] = a ;
				std::copy( begin_i + i, end_i, new_addition_state.begin() + i + 1 ) ;
				result->push_back( new_addition_state ) ;
			
				// Build the change states
				{
					std::size_t k = 0 ;
					for( ; k < i; ++k ) {
						std::copy( begin_i, begin_i + k, new_change_state.begin() ) ;
						std::copy( begin_i + k + 1, begin_i + i, new_change_state.begin() + k ) ;
						new_change_state[i-1] = a ;
						std::copy( begin_i + i, end_i, new_change_state.begin() + i ) ;
						result->push_back( new_change_state ) ;
					}
					for( ; k < m_current_state.size(); ++k ) {
						std::copy( begin_i, begin_i + i, new_change_state.begin() ) ;
						new_change_state[i] = a ;
						std::copy( begin_i + i, begin_i + k, new_change_state.begin() + i + 1 ) ;
						std::copy( begin_i + k + 1, end_i, new_change_state.begin() + k + 1 ) ;
						result->push_back( new_change_state ) ;
					}
				}
			}
		}
	}
}



