#ifndef GENFILE_GENOME_POSITION_RANGE_HPP
#define GENFILE_GENOME_POSITION_RANGE_HPP

#include "genfile/GenomePosition.hpp"

namespace genfile {
	// Class GenomePositionRange
	// Represents a closed, nonempty range of physical positions in the genome.
	struct GenomePositionRange {
		GenomePositionRange( GenomePosition start, GenomePosition end ) ;
		GenomePositionRange( GenomePositionRange const& ) ;
		GenomePositionRange& operator=( GenomePositionRange const& other ) ;
		
		GenomePosition const& get_start() const { return m_start ; }
		GenomePosition const& get_end() const { return m_end ; }

		bool check_if_contains( GenomePosition const& position ) const ;

		bool operator==( GenomePositionRange const& other ) const ;
		// Ranges are compared by their start and then their length.
		bool operator<( GenomePositionRange const& other ) const ;
		bool operator<=( GenomePositionRange const& other ) const ;

		private:
			GenomePosition m_start ;
			GenomePosition m_end ;

			// forbid default construction.
			GenomePositionRange() ;
	} ;
	
	std::ostream& operator<<( std::ostream&, GenomePositionRange const& ) ;
}

#endif
