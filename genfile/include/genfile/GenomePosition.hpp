#ifndef GENFILE_GENOME_POSITION_HPP
#define GENFILE_GENOME_POSITION_HPP

#include <iostream>
#include <stdint.h>
#include <utility>
#include "genfile/Chromosome.hpp"
#include "genfile/GenomePosition.hpp"

namespace genfile {
	typedef uint32_t Position ;

	struct GenomePosition
	{
		typedef genfile::Chromosome Chromosome ;
	
		GenomePosition() {}

		GenomePosition( Chromosome chromosome, Position position )
			: m_data( std::make_pair( chromosome, position ))
		{}

		GenomePosition( GenomePosition const& other )
			: m_data( other.m_data )
		{}
	
		GenomePosition& operator=( GenomePosition const& other ) {
			m_data = other.m_data ;
			return *this ;
		}

		// Load a genome position in one of the following formats:
		// CC:xxxxxx (chromosome / position pair)
		// xxxxx (unknown chromosome / position pair).
		GenomePosition( std::string const& position_spec ) ;
	
		Chromosome const& chromosome() const { return m_data.first ; }
		Chromosome& chromosome() { return m_data.first ; }
		Position const& position() const { return m_data.second ; }
		Position& position() { return m_data.second ; }
	
		friend bool operator<( GenomePosition const& left, GenomePosition const& right ) ;
		friend bool operator<=( GenomePosition const& left, GenomePosition const& right ) ;
		friend bool operator==( GenomePosition const& left, GenomePosition const& right ) ;
		friend bool operator!=( GenomePosition const& left, GenomePosition const& right ) ;

		static Position get_max_position( Chromosome chromosome ) ;
		static Position get_min_position( Chromosome chromosome ) ;
	private:
	
		std::pair< Chromosome, Position > m_data ;
	} ;

	std::ostream& operator<<( std::ostream&, GenomePosition const& ) ;
	std::istream& operator>>( std::istream& oStream, GenomePosition& pos ) ;
	
}

#endif
