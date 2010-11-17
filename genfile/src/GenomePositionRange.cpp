#include "genfile/Error.hpp"
#include "genfile/GenomePosition.hpp"
#include "genfile/GenomePositionRange.hpp"

namespace genfile {
	// Represents a closed, but possibly empty, range of genome positions.
	GenomePositionRange::GenomePositionRange( GenomePosition start, GenomePosition end ):
		m_start( start ),
		m_end( end )
	{
		if( m_end < m_start ) {
			throw BadArgumentError( "genfile::GenomePositionRange::GenomePositionRange( start, end )", "start, end" ) ;
		}
		
		if(( m_start.chromosome() == Chromosome() || m_start.chromosome() == Chromosome() ) && m_start.chromosome() != m_end.chromosome() ) {
			throw BadArgumentError( "genfile::GenomePositionRange::GenomePositionRange( start, end )", "start, end" ) ;
		}
	}
	
	GenomePositionRange::GenomePositionRange( GenomePositionRange const& other ):
		m_start( other.m_start ),
		m_end( other.m_end )
	{}

	GenomePositionRange& GenomePositionRange::operator=( GenomePositionRange const& other ) {
		m_start = other.m_start ;
		m_end = other.m_end ;
		return *this ;
	}
	
	bool GenomePositionRange::check_if_contains( GenomePosition const& position ) const {
		if( m_start.chromosome() == Chromosome() ) {
			return m_start.position() <= position.position() && position.position() <= m_end.position() ;
		}
		else {
			return ( m_start <= position ) && ( position <= m_end ) ;
		}
	}
	
	bool GenomePositionRange::operator==( GenomePositionRange const& other ) const {
		return m_start == other.m_start && m_end == other.m_end ;
	}
	// Ranges are compared by their start and then their length.
	bool GenomePositionRange::operator<( GenomePositionRange const& other ) const {
		return m_start < other.m_start || ( m_start == other.m_start && m_end < other.m_end ) ;
	}
	bool GenomePositionRange::operator<=( GenomePositionRange const& other ) const {
		return m_start < other.m_start || ( m_start == other.m_start && m_end <= other.m_end ) ;
	}
	
	std::ostream& operator<<( std::ostream& out, GenomePositionRange const& range ) {
		if( range.get_start().chromosome() != Chromosome() ) {
			out << range.get_start().chromosome() << ":" ;
		}
		if( range.get_start().position() > 0 ) {
			out << range.get_start().position() ;
		}
		out << "-" ;
		if( range.get_end().chromosome() != range.get_start().chromosome() && range.get_end().chromosome() != Chromosome() ) {
			out << range.get_end().chromosome() << ":" ;
		}
		if( range.get_end().position() < GenomePosition::get_max_position( range.get_end().chromosome() )) {
			out << range.get_end().position() ;
		}
		return out ;
	}
	
}
