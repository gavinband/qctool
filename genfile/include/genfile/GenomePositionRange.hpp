
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_GENOME_POSITION_RANGE_HPP
#define GENFILE_GENOME_POSITION_RANGE_HPP

#include "genfile/GenomePosition.hpp"

namespace genfile {
	// Class GenomePositionRange
	// Represents a closed, nonempty range of physical positions in the genome.
	struct GenomePositionRange {
		// parse a range in the form chr:start-end
		// or just start-end, whence the chromosome is treated as NA
		static GenomePositionRange parse( std::string const& spec ) ;

		GenomePositionRange( Position start, Position end ) ;
		GenomePositionRange( GenomePosition start, GenomePosition end ) ;
		GenomePositionRange( Chromosome chromosome, Position start, Position end ) ;
		GenomePositionRange( GenomePositionRange const& ) ;
		GenomePositionRange& operator=( GenomePositionRange const& other ) ;
		
		Chromosome const chromosome() const { return m_start.chromosome() ; }
		GenomePosition const& start() const { return m_start ; }
		GenomePosition const& end() const { return m_end ; }

		bool contains( GenomePosition const& position ) const ;

		bool operator==( GenomePositionRange const& other ) const ;
		// Ranges are compared by their start and then their length.
		bool operator<( GenomePositionRange const& other ) const ;
		bool operator<=( GenomePositionRange const& other ) const ;

		private:
			GenomePosition m_start ;
			GenomePosition m_end ;
			bool m_have_chromosome ;

			// forbid default construction.
			GenomePositionRange() ;
	} ;
	
	std::ostream& operator<<( std::ostream&, GenomePositionRange const& ) ;
}

#endif
