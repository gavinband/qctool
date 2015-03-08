
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <vector>
#include "genfile/string_utils.hpp"
#include "genfile/Error.hpp"
#include "genfile/GenomePosition.hpp"
#include "genfile/GenomePositionRange.hpp"

namespace genfile {
	GenomePositionRange GenomePositionRange::parse( std::string const& spec ) {
		std::vector< std::string > pieces = string_utils::split_and_strip( spec, ":", " " ) ;
		bool all_chromosomes = ( pieces.size() == 1 ) ;
		if ( all_chromosomes ) {
			pieces.insert( pieces.begin(), "NA" ) ;
		}
		if ( pieces.size() != 2 ) {
			throw genfile::BadArgumentError( "genfile::GenomePositionRange::parse", "spec=\"" + spec + "\"" ) ;
		}

		std::size_t separator_pos = pieces[1].find( '-' ) ;
		if ( separator_pos == std::string::npos ) {
			throw genfile::BadArgumentError( "genfile::GenomePositionRange::parse", "spec=\"" + spec + "\"" ) ;
		}

		genfile::Chromosome chromosome( pieces[0] ) ;

		genfile::Position start_position = GenomePosition::get_min_position( chromosome ) ; ;
		if ( separator_pos != 0 ) {
			start_position = genfile::string_utils::to_repr< genfile::Position >( pieces[1].substr( 0, separator_pos )) ;
		}

		genfile::Position end_position = GenomePosition::get_max_position( chromosome ) ;
		if ( separator_pos != pieces[1].size() - 1 ) {
			end_position = genfile::string_utils::to_repr< genfile::Position >( pieces[1].substr( separator_pos + 1, pieces[1].size() )) ;
		}

		if( all_chromosomes ) {
			return genfile::GenomePositionRange(
				start_position,
				end_position
			) ;
		}
		else {
			return genfile::GenomePositionRange(
				GenomePosition( chromosome, start_position ),
				GenomePosition( chromosome, end_position )
			) ;
		}
	}

	GenomePositionRange::GenomePositionRange( Position start, Position end ):
		m_start( GenomePosition( Chromosome(), start ) ),
		m_end( GenomePosition( Chromosome(), end ) ),
		m_have_chromosome( false )
	{
		if( m_end < m_start ) {
			throw BadArgumentError( "genfile::GenomePositionRange::GenomePositionRange( start, end )", "start, end" ) ;
		}
	}

	GenomePositionRange::GenomePositionRange( Chromosome chromosome, Position start, Position end ):
		m_start( GenomePosition( chromosome, start ) ),
		m_end( GenomePosition( chromosome, end ) ),
		m_have_chromosome( true )
	{
		if( m_end < m_start ) {
			throw BadArgumentError( "genfile::GenomePositionRange::GenomePositionRange( start, end )", "start, end" ) ;
		}
	}

	GenomePositionRange::GenomePositionRange( GenomePosition start, GenomePosition end ):
		m_start( start ),
		m_end( end ),
		m_have_chromosome( true )
	{
		if( m_start.chromosome() != m_end.chromosome() ) {
			throw BadArgumentError( "genfile::GenomePositionRange::GenomePositionRange( start, end )", "start, end" ) ;
		}

		if( m_end < m_start ) {
			throw BadArgumentError( "genfile::GenomePositionRange::GenomePositionRange( start, end )", "start, end" ) ;
		}
	}
	
	GenomePositionRange::GenomePositionRange( GenomePositionRange const& other ):
		m_start( other.m_start ),
		m_end( other.m_end ),
		m_have_chromosome( other.m_have_chromosome )
	{}

	GenomePositionRange& GenomePositionRange::operator=( GenomePositionRange const& other ) {
		m_start = other.m_start ;
		m_end = other.m_end ;
		return *this ;
	}
	
	bool GenomePositionRange::contains( GenomePosition const& position ) const {
		if( !m_have_chromosome ) {
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
		if( range.start().chromosome() != Chromosome() ) {
			out << range.start().chromosome() << ":" ;
		}
		if( range.start().position() > 0 ) {
			out << range.start().position() ;
		}
		out << "-" ;
		if( range.end().chromosome() != range.start().chromosome() && range.end().chromosome() != Chromosome() ) {
			out << range.end().chromosome() << ":" ;
		}
		if( range.end().position() < GenomePosition::get_max_position( range.end().chromosome() )) {
			out << range.end().position() ;
		}
		return out ;
	}
	
}
