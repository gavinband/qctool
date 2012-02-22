#ifndef GENFILE_FROM_FILES_GENETIC_MAP_HPP
#define GENFILE_FROM_FILES_GENETIC_MAP_HPP

#include <memory>
#include <iostream>
#include <utility>
#include <set>
#include <map>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/function.hpp>
#include "genfile/wildcard.hpp"
#include "genfile/GenomePosition.hpp"
#include "genfile/GeneticMap.hpp"

namespace genfile {
	struct FromFilesGeneticMap: public GeneticMap
	{
	public:
		typedef boost::function< void ( std::size_t, std::size_t ) > ProgressCallback ;
	public:
		FromFilesGeneticMap( std::string const& filename, ProgressCallback progress_callback = ProgressCallback() ) ;
		FromFilesGeneticMap( std::vector< genfile::wildcard::FilenameMatch > const& filenames, ProgressCallback progress_callback = ProgressCallback() ) ;
		double calculate_cM_between_positions( GenomePosition position1, GenomePosition position2 ) const ;
		double find_cM_from_beginning_of_chromosome_at_position( GenomePosition position ) const ;
		GenomePosition find_least_physical_position( Chromosome const& chromosome, double const cM ) const ;
		GenomePosition find_greatest_physical_position( Chromosome const& chromosome, double const cM ) const ;
		double find_rate_at_position( GenomePosition position ) const ;
		std::vector< genfile::GenomePosition > get_measured_positions() const ;
		std::set< Chromosome > get_chromosomes() const ;
		double get_end_of_map_in_cM( Chromosome const& chromosome ) const ;
		double get_start_of_map_in_cM( Chromosome const& chromosome ) const ;
		double get_length_of_genome_in_cM() const ;

	private:
		void setup( std::vector< genfile::wildcard::FilenameMatch > const& ) ;
		//typedef std::map< Position, std::pair< double, double > > PerChromosomeMap ;
		struct MapPoint {
			MapPoint() ;
			MapPoint( Position position_in_bp, double position_in_cM, double rate ) ;
			MapPoint( MapPoint const& other ) ;

			Position position_in_bp ;
			double position_in_cM ;
			double rate ;

			struct CompareByPhysicalPosition {
				typedef Position first_argument_type ;
				typedef Position second_argument_type ;
				typedef bool result_type ;
				bool operator()( MapPoint const& left, MapPoint const& right ) const {
					return ( left.position_in_bp < right.position_in_bp )
						|| ( left.position_in_bp == right.position_in_bp && left.position_in_cM < right.position_in_cM ) ;
				}
			} ;

			struct CompareByRecombinationPosition {
				typedef MapPoint first_argument_type ;
				typedef MapPoint second_argument_type ;
				typedef bool result_type ;
				bool operator()( MapPoint const& left, MapPoint const& right ) const {
					return ( left.position_in_cM < right.position_in_cM )
						|| ( left.position_in_cM == right.position_in_cM && left.position_in_bp < right.position_in_bp ) ;
				}
			} ;
		} ;
		
		typedef boost::multi_index_container<
			MapPoint,
			boost::multi_index::indexed_by<
				boost::multi_index::ordered_unique< boost::multi_index::identity< MapPoint >, MapPoint::CompareByPhysicalPosition >,
				boost::multi_index::ordered_unique< boost::multi_index::identity< MapPoint >, MapPoint::CompareByRecombinationPosition >
			>
		> PerChromosomeMap ;
		typedef PerChromosomeMap::nth_index<0>::type const& PhysicalToRecombinationDistanceMap ;
		typedef PerChromosomeMap::nth_index<1>::type const& RecombinationToPhysicalDistanceMap ;
		typedef std::map< Chromosome, PerChromosomeMap > Map ;

		double find_cM_from_beginning_of_chromosome_at_position( PerChromosomeMap const&, Position position ) const ;
		void add_entries( std::istream& source ) ;
		void add_entries( std::istream& source, Chromosome const& chromosome ) ;
		double interpolate( double x1, double y1, double x2, double y2, double x ) const ;

	 	Map m_map ;
		ProgressCallback const m_progress_callback ;
	} ;	
}

#endif
