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
	
		Chromosome const& chromosome() const { return m_data.first ; }
		Chromosome& chromosome() { return m_data.first ; }
		Position const& position() const { return m_data.second ; }
		Position& position() { return m_data.second ; }
	
		friend bool operator<( GenomePosition const& left, GenomePosition const& right ) ;
		friend bool operator<=( GenomePosition const& left, GenomePosition const& right ) ;
		friend bool operator==( GenomePosition const& left, GenomePosition const& right ) ;
		friend bool operator!=( GenomePosition const& left, GenomePosition const& right ) ;
	private:
	
		std::pair< Chromosome, Position > m_data ;
	} ;

	std::ostream& operator<<( std::ostream&, GenomePosition const& ) ;
}

#endif
