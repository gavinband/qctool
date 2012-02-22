#ifndef GENFILE_GENETIC_MAP_HPP
#define GENFILE_GENETIC_MAP_HPP

#include <memory>
#include <iostream>
#include <utility>
#include <set>
#include "genfile/wildcard.hpp"
#include "genfile/GenomePosition.hpp"

namespace genfile {
	struct GeneticMapError: public std::exception { char const* what() const throw() { return "genpos::GeneticMapError" ; } } ;

	struct GeneticMapIsEmptyError: public GeneticMapError { char const* what() const throw() { return "genpos::GeneticMapIsEmptyError" ; } } ;
	struct ChromosomeNotInMapError: public GeneticMapError
	{
		ChromosomeNotInMapError( genfile::Chromosome  const& chromosome ): m_chromosome( chromosome ) {}
		~ChromosomeNotInMapError() throw() {} ;
		char const* what() const throw() { return "genpos::ChromosomeNotInMapError" ; }
		genfile::Chromosome const& chromosome() const { return m_chromosome ;}
	private:
		genfile::Chromosome const m_chromosome ;
	} ;

	struct GeneticMap
	{
		typedef std::auto_ptr< GeneticMap > UniquePtr ;
		typedef genfile::Chromosome Chromosome ;
		
		virtual ~GeneticMap() {}
		
		virtual double calculate_cM_between_positions( GenomePosition position1, GenomePosition position2 ) const = 0 ;
		virtual GenomePosition find_least_physical_position( Chromosome const& chromosome, double const cM ) const = 0 ;
		virtual GenomePosition find_greatest_physical_position( Chromosome const& chromosome, double const cM ) const = 0 ;
		virtual double find_cM_from_beginning_of_chromosome_at_position( GenomePosition position ) const = 0 ;
		virtual double find_rate_at_position( GenomePosition position ) const = 0 ;
		virtual std::set< Chromosome > get_chromosomes() const = 0 ;
		virtual double get_start_of_map_in_cM( Chromosome const& chromosome ) const = 0 ;
		virtual double get_end_of_map_in_cM( Chromosome const& chromosome ) const = 0 ;
		virtual double get_length_of_genome_in_cM() const = 0 ;
	} ;
}

#endif
