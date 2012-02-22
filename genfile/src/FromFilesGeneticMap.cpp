#include <memory>
#include <iostream>
#include <utility>
#include <cassert>
#include <limits>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/multi_index/ordered_index.hpp>

#include "genfile/FromFilesGeneticMap.hpp"
#include "genfile/wildcard.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/string_utils/slice.hpp"
#include "genfile/Error.hpp"

namespace genfile {
	
	FromFilesGeneticMap::MapPoint::MapPoint():
		position_in_bp(0),
		position_in_cM(0),
		rate(0)
	{}
	
	FromFilesGeneticMap::MapPoint::MapPoint( Position position_in_bp_, double position_in_cM_, double rate_ ):
		position_in_bp( position_in_bp_ ),
		position_in_cM( position_in_cM_ ),
		rate( rate_ )
	{}

	FromFilesGeneticMap::MapPoint::MapPoint( MapPoint const& other ):
		position_in_bp( other.position_in_bp ),
		position_in_cM( other.position_in_cM ),
		rate( other.rate )
	{}
	
	FromFilesGeneticMap::FromFilesGeneticMap( std::vector< genfile::wildcard::FilenameMatch > const& filenames, ProgressCallback progress_callback )
	: m_progress_callback( progress_callback )
	{
		assert( filenames.size() > 0 ) ;
		setup( filenames ) ;
	}

	FromFilesGeneticMap::FromFilesGeneticMap( std::string const& filename, ProgressCallback progress_callback )
		: m_progress_callback( progress_callback )
	{
		std::vector< genfile::wildcard::FilenameMatch > matches( genfile::wildcard::find_files_by_chromosome( filename )) ;
		setup( matches ) ;
	}

	void FromFilesGeneticMap::setup( std::vector< genfile::wildcard::FilenameMatch > const& matches ) {
		boost::ptr_vector< std::istream > sources ;
		std::size_t total_rows = 0 ;
		for( std::size_t i = 0; i < matches.size(); ++i ) {
			std::string const& filename = matches[i].filename() ;
			sources.push_back( open_text_file_for_input( filename ).release() ) ;
		}
		assert( sources.size() == matches.size() ) ;

		std::size_t rows_complete = 0 ;
		for( std::size_t i = 0; i < sources.size(); ++i ) {
			std::string line ;
			std::getline( sources[i], line ) ;
			std::vector< string_utils::slice > elts = string_utils::slice( line ).split( " " ) ;

			// Do we have a chromosome column?
			if( elts.size() > 0 && ( elts[0] == "chromosome" || elts[0] == "chr" ) ) {
				add_entries( sources[i] ) ;
			}
			else {
				// Uh-oh, no chromosome.  Attempt to infer it.
				genfile::Chromosome chromosome = matches[i].match() ;
				if( chromosome == Chromosome() ) {
					std::size_t pos = matches[i].filename().find( "_chr") ;
					if( pos != std::string::npos && matches[i].filename().size() > pos + 4 ) {
						std::string chr_string ;
						chr_string.push_back( matches[i].filename()[pos+4] ) ;
						if( chr_string[0] == '0' && matches[i].filename().size() > pos + 5 ) {
							chr_string.push_back( matches[i].filename()[pos+5] ) ;
						}
						chromosome = chr_string ;
					}
				}
				add_entries( sources[i], chromosome ) ;
			}
			if( m_progress_callback ) {
				m_progress_callback( i, matches.size() ) ;
			}
		}
		if( m_progress_callback ) {
			m_progress_callback( matches.size(), matches.size() ) ;
		}
	}

	void FromFilesGeneticMap::add_entries(
		std::istream& source
	) {
		std::string chromosome ;
		Position position ;
		double rate_in_cM_per_Mb, cM ;
		while( source >> chromosome >> position >> rate_in_cM_per_Mb >> cM ) {
			PerChromosomeMap& map = m_map[  genfile::Chromosome( chromosome ) ] ;
			// PerChromosomeMap::nth_index<0>::type& index = map.get<0>() ;
			map.insert( map.end(), MapPoint( position, cM, rate_in_cM_per_Mb )) ;
			// discard rest of line.
			std::string line ;
			std::getline( source, line ) ;
		}
	}

	void FromFilesGeneticMap::add_entries(
		std::istream& source,
		Chromosome const& chromosome
	) {
		Position position ;
		double combined_rate, cM ;
		PerChromosomeMap& map = m_map[ chromosome ] ;
		while( source >> position >> combined_rate >> cM ) {
			map.insert( map.end(), MapPoint( position, cM, combined_rate )) ;
			// discard rest of line.
			std::string line ;
			std::getline( source, line ) ;
		}
	}

	double FromFilesGeneticMap::calculate_cM_between_positions( GenomePosition position1, GenomePosition position2 ) const {
		assert( position1 <= position2) ;
		if( position1.chromosome() != position2.chromosome() ) {
			return std::numeric_limits< double >::infinity() ;
		}
		else {
			return find_cM_from_beginning_of_chromosome_at_position( position2 ) - find_cM_from_beginning_of_chromosome_at_position( position1 ) ; 
		}
	}

	double FromFilesGeneticMap::find_cM_from_beginning_of_chromosome_at_position( GenomePosition position ) const {
		Map::const_iterator where = m_map.find( position.chromosome() ) ;
		if( where == m_map.end() ) {
			throw ChromosomeNotInMapError( position.chromosome() ) ;
		}
		return find_cM_from_beginning_of_chromosome_at_position( where->second, position.position() ) ;
	}

	double FromFilesGeneticMap::find_cM_from_beginning_of_chromosome_at_position( PerChromosomeMap const& map, Position position ) const {
		PerChromosomeMap::const_iterator upper_bound = map.get<0>().lower_bound( MapPoint( position, 0, 0 ) ) ;

		if( upper_bound == map.end() ) {
			// Clamp to the highest map value (note: we assume map is not empty)
			--upper_bound ;
			return upper_bound->position_in_cM ;
		}
		else if( upper_bound == map.begin() ) {
			// Clamp to the lowest map value
			return upper_bound->position_in_cM ;
		}
		else {
			PerChromosomeMap::const_iterator lower_bound = upper_bound ;
			--lower_bound ;
			return interpolate(
				lower_bound->position_in_bp,
				lower_bound->position_in_cM,
				upper_bound->position_in_bp,
				upper_bound->position_in_cM,
				position
			) ;
		}
	}

	double FromFilesGeneticMap::interpolate( double x1, double y1, double x2, double y2, double x ) const {
		assert( x >= x1 && x <= x2 ) ;
		// interpolate as (1-t) * lower value + t * upper value, t in [0,1].
		double t ;
		if( y2 == y1 ) {
			return y2 ;
		}
		else if( x2 == x1 ) {
			assert( y2 == y1 ) ;
			t = 0.0 ;
		}
		else {
			t = ( static_cast< double >( x - x1 ) / static_cast< double >( x2 - x1 )) ;
		}
		double result = ( (1.0-t) * y1 ) + ( t * y2 ) ;
		if( result < 0.0 ) {
			std::cerr << result << "!\n" ;
			assert( 0 ) ;
		}
		return result ;
	}
	
	GenomePosition FromFilesGeneticMap::find_least_physical_position( Chromosome const& chromosome, double const cM ) const {
		GenomePosition result( chromosome, 0 ) ;
		Map::const_iterator chromosome_map = m_map.find( chromosome ) ;
		if( chromosome_map == m_map.end() ) {
			throw ChromosomeNotInMapError( chromosome ) ;
		}
		typedef PerChromosomeMap::nth_index<1>::type::const_iterator Iterator ;
		PerChromosomeMap::nth_index<1>::type const& index = chromosome_map->second.get<1>() ;
		Iterator where = index.lower_bound( MapPoint( 0, cM, 0 ) ) ;
		if( where == index.end() || ( where == index.begin() && where->position_in_cM > cM ) ) {
			throw genfile::BadArgumentError( "genfile::FromFilesGeneticMap::find_least_physical_position()", "cM=" + genfile::string_utils::to_string( cM )) ;
		}
		else if( where->position_in_cM == cM ) {
			result.position() = where->position_in_bp ;
		}
		else {
			Iterator previous = where ;
			--previous ;
			result.position() = interpolate(
				previous->position_in_cM,
				previous->position_in_bp,
				where->position_in_cM,
				where->position_in_bp,
				cM
			) ;
		}
		return result ;
	}

	GenomePosition FromFilesGeneticMap::find_greatest_physical_position( Chromosome const& chromosome, double const cM ) const {
		GenomePosition result( chromosome, 0 ) ;
		Map::const_iterator chromosome_map = m_map.find( chromosome ) ;
		if( chromosome_map == m_map.end() ) {
			throw ChromosomeNotInMapError( chromosome ) ;
		}
		typedef PerChromosomeMap::nth_index<1>::type::const_iterator Iterator ;
		PerChromosomeMap::nth_index<1>::type const& index = chromosome_map->second.get<1>() ;
		Iterator where = index.upper_bound( MapPoint( 0, cM, 0 ) ) ;
		if( where == index.end() || where == index.begin() ) {
			throw genfile::BadArgumentError( "genfile::FromFilesGeneticMap::find_least_physical_position()", "cM=" + genfile::string_utils::to_string( cM )) ;
		}
		else if( where->position_in_cM == cM ) {
			result.position() = where->position_in_bp ;
		}
		else {
			Iterator previous = where ;
			--previous ;
			result.position() = interpolate(
				previous->position_in_cM,
				previous->position_in_bp,
				where->position_in_cM,
				where->position_in_bp,
				cM
			) ;
		}
		return result ;
	}
	
	double FromFilesGeneticMap::find_rate_at_position( GenomePosition position ) const {
		Map::const_iterator where = m_map.find( position.chromosome() ) ;
		if( where == m_map.end() ) {
			throw ChromosomeNotInMapError( position.chromosome() ) ;
		}
		PerChromosomeMap::const_iterator position_i = where->second.lower_bound( MapPoint( position.position(), 0, 0 ) ) ;
		if( position_i != where->second.end() && !( position.position() < position_i->position_in_bp )) {
			// found it exactly
			return position_i->position_in_cM ;
		}
		else if( position_i == where->second.begin() ) {
				return 0.0 ;
		}
		else {
			--position_i ;
			return position_i->position_in_cM ;
		}
	}
	
	std::vector< genfile::GenomePosition > FromFilesGeneticMap::get_measured_positions() const {
		std::vector< genfile::GenomePosition > result ;
		for( Map::const_iterator c_i = m_map.begin() ; c_i != m_map.end(); ++c_i ) {
			for( PerChromosomeMap::const_iterator p_i = c_i->second.begin(); p_i != c_i->second.end(); ++p_i ) {
				result.push_back( genfile::GenomePosition( c_i->first, p_i->position_in_cM )) ;
			}
		}
		return result ;
	}
	
	std::set< Chromosome > FromFilesGeneticMap::get_chromosomes() const {
		std::set< Chromosome > result ;
		for( Map::const_iterator c_i = m_map.begin() ; c_i != m_map.end(); ++c_i ) {
			result.insert( c_i->first ) ;
		}
		return result ;
	}
	
	double FromFilesGeneticMap::get_start_of_map_in_cM( Chromosome const& chromosome ) const {
		Map::const_iterator where = m_map.find( chromosome ) ;
		if( where == m_map.end() ) {
			throw ChromosomeNotInMapError( chromosome ) ;
		}
		return where->second.begin()->position_in_cM ;
	}
	
	double FromFilesGeneticMap::get_end_of_map_in_cM( Chromosome const& chromosome ) const {
		Map::const_iterator where = m_map.find( chromosome ) ;
		if( where == m_map.end() ) {
			throw ChromosomeNotInMapError( chromosome ) ;
		}
		assert( where->second.size() > 0 ) ;
		PerChromosomeMap::const_iterator last = where->second.end() ;
		--last ;
		return last->position_in_cM ;
	}
	
	double FromFilesGeneticMap::get_length_of_genome_in_cM() const {
		std::set< Chromosome > const chromosomes = get_chromosomes() ;
		double result = 0.0 ;
		for( std::set< Chromosome >::const_iterator i = chromosomes.begin(); i != chromosomes.end(); ++i ) {
			result += get_end_of_map_in_cM( *i ) ;
		}
		return result ;
	}

	std::string FromFilesGeneticMap::get_summary() const {
		std::string result = "genetic map of total length " + genfile::string_utils::to_string( get_length_of_genome_in_cM() ) ;
		return result ;	
	}
}


