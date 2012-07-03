
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <vector>
#include <utility>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <ctime>
#include <iostream>
#include <deque>
#include <algorithm>

#include <boost/bind.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/filesystem.hpp>

#include "appcontext/ProgramFlow.hpp"
#include "appcontext/CmdLineOptionProcessor.hpp"
#include "appcontext/ApplicationContext.hpp"
#include "appcontext/get_current_time_as_string.hpp"

#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPFilteringSNPDataSource.hpp"
#include "genfile/get_list_of_snps_in_source.hpp"
#include "genfile/utility.hpp"
#include "genfile/CommonSNPFilter.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/string_utils.hpp"

#include "statfile/BuiltInTypeStatSourceChain.hpp"

#include "db/SQLite3Connection.hpp"
#include "genfile/FromFilesGeneticMap.hpp"
#include "FileUtil.hpp"
#include "components/SNPSummaryComponent/DBOutputter.hpp"

namespace globals {
	std::string const program_name = "inthinnerator" ;
}

struct InthinneratorOptionProcessor: public appcontext::CmdLineOptionProcessor
{
	// Methods needed for CmdLineOptionProcessor::process()
	std::string get_program_name() const { return globals::program_name ; }
	
	void declare_options( OptionProcessor& options ) {
		options.set_help_option( "-help" ) ;

		options.declare_group( "Genetic map options" ) ;
		options[ "-map" ]
			.set_description( "Set the path of the genetic map panel to use." )
			.set_takes_single_value()
			.set_default_value( std::string( getenv( "HOME" )) + "/.association-tools/genetic_map_b36.bnv" ) ;

		options.declare_group( "Genotype file options" ) ;
		options[ "-g" ]
			.set_description( "Specify a file containing the SNPs to operate on." )
			.set_is_required()
			.set_takes_values( 1 )
			.set_maximum_multiplicity( 100 ) ;

		options.declare_group( "SNP selection options" ) ;
		options[ "-excl-rsids" ]
			.set_description( "Specify a file containing a whitespace-separated list of SNP rsids."
				" SNPs with ids in this file will be excluded from the analysis." )
			.set_takes_values( 1 )
			.set_maximum_multiplicity( 100 ) ;
		options[ "-incl-rsids" ]
			.set_description( "Specify a file containing a whitespace-separated list of SNP rsids."
				" SNPs with ids not in this file will be excluded from the analysis." )
			.set_takes_single_value() ;
		options[ "-excl-snpids" ]
			.set_description( "Specify a file containing a whitespace-separated list of SNP SNPIDs."
			" SNPs with ids in this file will be excluded from the analysis." )
			.set_takes_values( 1 )
			.set_maximum_multiplicity( 100 ) ;
		options[ "-incl-snpids" ]
			.set_description( "Specify a file containing a whitespace-separated list of SNP SNPIDs."
			" SNPs with ids not in this file will be excluded from the analysis." )
			.set_takes_single_value() ;
		options[ "-incl-range" ]
			.set_description( "Specify a range of SNPs (or comma-separated list of ranges of SNPs) to operate on. "
				"Each range should be in the format CC:xxxx-yyyy where CC is the chromosome and xxxx and yyyy are the "
				"start and end coordinates, or just xxxx-yyyy which matches that range from all chromosomes. "
				"You can also omit either of xxxx or yyyy to get all SNPs from the start or to the end of a chromosome." )
			.set_takes_single_value() ;
		options[ "-excl-range" ]
			.set_description( "Specify a range of SNPs (or comma-separated list of ranges of SNPs) to exclude from operation. "
				"Each range should be in the format CC:xxxx-yyyy where CC is the chromosome and xxxx and yyyy are the "
				"start and end coordinates, or just xxxx-yyyy which matches that range from all chromosomes. "
				"You can also omit either of xxxx or yyyy to get all SNPs from the start or to the end of a chromosome." )
			.set_takes_single_value() ;

		options.declare_group( "SNP thinning options" ) ;
		options[ "-min-distance" ]
			.set_description( "Specify a minimum acceptable distance between SNPs."
				" This must be in one of the forms \"<X>cM\", \"<X>M\", \"<X>bp\", \"<X>kb\", or \"<X>Mb\"."
				" The thinned list of SNPs will contain no SNPs within this distance of each other."
				" For recombination-distance offsets, a physical margin can also be specified as in \"<X>cM+<Y>kb\"." )
			.set_takes_single_value()
			.set_default_value( "0.01cM" ) ;
		options[ "-strategy" ]
			.set_description( "Specify the SNP thinning strategy if not using -rank."
				" This can be \"random\" or \"first\"." )
			.set_takes_single_value()
			.set_default_value( "random" ) ;
		options[ "-rank" ]
			.set_description( "Specify name of a file containing numerical ranks of SNPs. "
				"SNPs will be picked in nonincreasing order of rank. "
				"This file must have first five columns SNPID, rsid, chromosome, position, allele1, allele2."
			)
			.set_takes_single_value() ;
		options[ "-rank-col" ]
			.set_description( "Specify the name of the column in the file for -rank containing the ranks." )
			.set_takes_single_value() ;
		options[ "-missing-code" ]
			.set_description( "Specify a comma-separated list of strings to be treated as missing ranks. "
				"Missing ranks are always picked last." )
			.set_takes_single_value()
			.set_default_value( "NA" ) ;

		options.option_excludes_option( "-rank", "-strategy" ) ;
		options.option_excludes_option( "-strategy", "-rank" ) ;
		options.option_implies_option( "-rank", "-rank-col" ) ;
		options.option_implies_option( "-missing-code", "-rank" ) ;

		options[ "-N" ]
			.set_description( "Specify a number of thinnings to perform.  This must be at least 1."
				" If this is larger than one, output files will be numbered in the form <filename>.####" )
			.set_takes_single_value()
			.set_default_value( 1 ) ;

		options[ "-max-picks" ]
			.set_description( "Specify a number of SNPs to pick in each thinning."
				" By default we choose as many SNPs as it's possible to choose subject to the minimum distance constraint." )
			.set_takes_single_value()
			.set_default_value( std::numeric_limits< std::size_t >::max() ) ;

		options.declare_group( "Output file options" ) ;
		options[ "-o" ]
			.set_description( "Specify the output filename stub." )
			.set_takes_single_value() ;
		options[ "-odb" ]
			.set_description( "Specify the name of a database file to output." )
			.set_takes_single_value() ;
		options[ "-output-cols" ]
			.set_description( "Specify a comma-separated list of columns that should appear in the output files."
				" Possible columns are \"SNPID\", \"rsid\", \"chromosome\", \"position\", \"allele1\", \"allele2\", and \"cM_from_start_of_chromosome\"."
				" The special value \"all\" indicates that all available columns will be output." )
			.set_takes_single_value()
			.set_default_value( "all" ) ;
		options[ "-headers" ]
			.set_description( "Specify this to force output of column headersomit column headers in the output files." ) ;
		options[ "-no-headers" ]
			.set_description( "Specify this to suppress output of column headers in the output files." ) ;

		options.option_excludes_option( "-o", "-odb" ) ;
		options.option_excludes_option( "-odb", "-o" ) ;
		options.option_excludes_option( "-headers", "-no-headers" ) ;
		options.option_implies_option( "-headers", "-o" ) ;
		options.option_implies_option( "-no-headers", "-o" ) ;
		
		options[ "-write-excl-list" ]
			.set_description( "Specify that inthinnerator should write out inclusion lists.  On by default." ) ;
		options[ "-write-incl-list" ]
			.set_description( "Specify that inthinnerator should produce out exclusion lists.  On by default." ) ;
		
		options.declare_group( "Miscellaneous options" ) ;
		options[ "-log" ]
			.set_description( "Set the path of the log file to output." )
			.set_takes_single_value()
			.set_default_value( globals::program_name + ".log" ) ;
		options[ "-analysis-name" ]
			.set_description( "Specify a name to label results from this analysis with.  (This applies to modules which store their results in a qcdb file.)" )
			.set_takes_single_value()
			.set_default_value( globals::program_name + " analysis, started " + appcontext::get_current_time_as_string() ) ;
	}
} ;


// Base class for snp pickers
class SNPPicker
{
public:
	typedef std::auto_ptr< SNPPicker > UniquePtr ;
	virtual ~SNPPicker() {}
	virtual std::size_t pick(
		std::deque< std::size_t > const& among_these
	) const = 0 ;
	virtual void set_sort_order( std::deque< std::size_t > const& ) {}
	virtual std::string display() const = 0 ;
	virtual std::set< std::string > get_attribute_names() const = 0 ;
	virtual std::map< std::string, genfile::VariantEntry > get_attributes( std::size_t chosen_snp ) const = 0 ;
} ;

// Pick first available SNP
class FirstAvailableSNPPicker: public SNPPicker
{
public:
	FirstAvailableSNPPicker() {} ;

	virtual std::size_t pick(
		std::deque< std::size_t > const& among_these
	) const {
		assert( among_these.size() > 0 ) ;
		return *(among_these.begin()) ;
	}
	
	std::string display() const {
		return "FirstAvailableSNPPicker" ;
	}
	
	std::set< std::string > get_attribute_names() const {
		return std::set< std::string >() ;
	}
	std::map< std::string, genfile::VariantEntry > get_attributes( std::size_t chosen_snp ) const {
		return std::map< std::string, genfile::VariantEntry >() ;
	} ;
	
	
} ;

// Pick a random SNP
class RandomSNPPicker: public SNPPicker
{
private:
	typedef boost::mt19937 RNG ;
	typedef boost::uniform_int< std::size_t > Distribution ;
public:
	RandomSNPPicker():
		m_rng( new RNG( get_random_seed() ) )
	{
	}

	RandomSNPPicker( std::size_t seed ):
		m_rng( new RNG( seed ) )
	{}

	virtual std::size_t pick(
		std::deque< std::size_t > const& among_these
	) const {
		assert( among_these.size() > 0 ) ;
		Distribution distribution( 0, among_these.size() - 1 ) ;
		std::size_t choice = distribution( *m_rng ) ;
		assert( choice < among_these.size() ) ;
		std::deque< std::size_t >::const_iterator i = among_these.begin() ;
		std::advance( i, choice ) ;
		return *i ;
	}

	std::string display() const {
		return "RandomSNPPicker" ;
	}
	
	std::set< std::string > get_attribute_names() const {
		return std::set< std::string >() ;
	}
	std::map< std::string, genfile::VariantEntry > get_attributes( std::size_t chosen_snp ) const {
		return std::map< std::string, genfile::VariantEntry >() ;
	} ;

private:
	std::auto_ptr< RNG > m_rng ;
	
	std::size_t get_random_seed() const {
		std::size_t seed ;
		if( boost::filesystem::exists( "/dev/random" )) {
			std::ifstream ifs( "/dev/random" ) ;
			if( !ifs.is_open() ) {
				throw genfile::ResourceNotOpenedError( "/dev/random" ) ;
			}
			char buf[ sizeof( std::size_t ) ] ;
			ifs.read( buf, sizeof( std::size_t )) ;
			ifs.close() ;
			seed = *(reinterpret_cast< std::size_t* >( buf )) ;
		}
		else {
			seed = std::time(0) ;
		}
		return seed ;
	}
} ;

class HighestValueSNPPicker: public SNPPicker
{
private:
	struct DoubleComparator {
		bool operator()( double const a, double const b ) const {
			// Compare doubles putting NaN's last of all.
			// From http://stackoverflow.com/questions/4816156/are-ieee-floats-valid-key-types-for-stdmap-and-stdset
		    if ((a == a) && (b == b)) {
		        return a < b ;
		    }
		    if ((a != a) && (b != b)) return false ;
		    // We have one NaN and one non-NaN.
		    // Let's say NaN is less than everything
			return (a != a) ;
		}
	} ;

public:
	typedef std::auto_ptr< SNPPicker > UniquePtr ;
	HighestValueSNPPicker(
		std::vector< double > const& values
	):
		m_values( values ),
		m_value_map( DoubleComparator() )
	{
		for( std::size_t i = 0; i < m_values.size(); ++i ) {
			m_value_map.insert( std::make_pair( m_values[ i ], i ) ) ;
		}
	}

	// the purpose of set_sort_order is to tell HighestValueSNPPicker
	// the sorted order of SNPs so that it can treat the list of SNPs as sorted.
	void set_sort_order( std::deque< std::size_t > const& list ) {
		m_positions_in_sorted_list.resize( list.size(), std::numeric_limits< std::size_t >::max() ) ;
		for( std::size_t i = 0; i < list.size(); ++i ) {
			assert( list[i] < m_positions_in_sorted_list.size() ) ;
			assert( m_positions_in_sorted_list[ list[i] ] == std::numeric_limits< std::size_t >::max() ) ;
			m_positions_in_sorted_list[ list[i] ] = i ;
		}
	}

	// Compare SNP indices using their sorted order.
	bool compare_indices( std::size_t i, std::size_t j ) const {
		assert( i < m_positions_in_sorted_list.size() ) ;
		assert( j < m_positions_in_sorted_list.size() ) ;
		return m_positions_in_sorted_list[i] < m_positions_in_sorted_list[j] ;
	}
	
	std::size_t pick(
		std::deque< std::size_t > const& among_these
	) const {
		assert( among_these.size() > 0 ) ;
		ValueMap::reverse_iterator pick = m_value_map.rbegin() ;
		//std::cerr << "pick is " << pick->second << ", among_these has " << among_these.size() << " elts, \n" ;
		//for( std::size_t i = 0; i < among_these.size(); ++i ) {
		//	std::cerr << among_these[i] << "\n" ;//": " << m_positions_in_sorted_list[ among_these[i] ] << ".\n" ;
		//}
		while(
			!std::binary_search(
				among_these.begin(),
				among_these.end(),
				pick->second,
				boost::bind(
					&HighestValueSNPPicker::compare_indices,
					this,
					_1,
					_2
				)
			)
		) {
			++pick ;
			assert( pick != m_value_map.rend() ) ;
		}

        std::size_t chosen_snp = pick->second ;
		
		// The following line is the main optimisation which ensures we don't keep repeating work
		// we've already done.
		// Warning! This line assumes that the set of SNPs passed to this function
		// will always decrease (in the sense of set theory.)  If not, the assert above should
		// be triggered.
		m_value_map.erase( pick.base(), m_value_map.end() ) ;

        return chosen_snp ;
	}

	std::string display() const {
		return "HighestValueSNPPicker" ;
	} ;
	
	std::set< std::string > get_attribute_names() const {
		std::set< std::string > result ;
		result.insert( "rank" ) ;
		return result ;
	}
	
	std::map< std::string, genfile::VariantEntry > get_attributes( std::size_t chosen_snp ) const {
		std::map< std::string, genfile::VariantEntry > result ;
		result[ "rank" ] = m_values[ chosen_snp ] ;
		return result ;
	}
	
private:
	
	
	std::vector< double > const m_values ;
	typedef std::multimap< double, std::size_t, DoubleComparator > ValueMap ;
	mutable ValueMap m_value_map ;
	std::vector< std::size_t > m_positions_in_sorted_list ;
} ;

// Base class for too-close measurements
class ProximityTest
{
public:
	typedef std::auto_ptr< ProximityTest > UniquePtr ;
	
public:
	virtual ~ProximityTest() {}
	virtual void prepare( std::deque< std::size_t >* among_these ) const = 0 ;
	virtual void remove_snps_too_close_to( std::size_t chosen_snp_i, std::deque< std::size_t >* among_these ) const = 0 ;	
	virtual std::string display() const = 0 ;
	virtual std::set< std::string > get_attribute_names() const = 0 ;
	virtual std::map< std::string, genfile::VariantEntry > get_attributes( std::size_t ) const = 0 ;
} ;

class RecombinationDistanceProximityTest: public ProximityTest
{
public:
	RecombinationDistanceProximityTest(
		std::vector< genfile::SNPIdentifyingData > const& snps,
		genfile::GeneticMap const& genetic_map,
		double minimum_distance_in_cM,
		std::size_t margin_in_bp
	):
		m_snps( snps ),
		m_genetic_map( genetic_map ),
		m_minimum_distance_in_cM( minimum_distance_in_cM ),
		m_margin_in_bp( margin_in_bp )
	{
	}
	
	std::string display() const {
		std::string result = "RecombinationDistanceProximityTest( " + genfile::string_utils::to_string( m_minimum_distance_in_cM ) + "cM" ;
		if( m_margin_in_bp > 0 ) {
			result += "+" + genfile::string_utils::to_string( m_margin_in_bp ) + "bp" ;
		}
		result += " )" ;
		return result ;
	}

	void prepare( std::deque< std::size_t >* among_these ) const {
		// Sort by chromosome / recombination distance
		//std::sort( among_these->begin(), among_these->end(), boost::bind( &RecombinationDistanceProximityTest::compare_recombination_positions, this, _1, _2 )) ;
		std::sort( among_these->begin(), among_these->end(), boost::bind( &RecombinationDistanceProximityTest::compare_physical_positions, this, _1, _2 )) ;
	}

	void remove_snps_too_close_to( std::size_t chosen_snp, std::deque< std::size_t >* among_these ) const {
		std::pair< genfile::Chromosome, double >
			lower_bound = get_recombination_position( m_snps[ chosen_snp ].get_position() ),
			upper_bound = lower_bound ;
		lower_bound.second -= m_minimum_distance_in_cM ;
		lower_bound.second = std::max( lower_bound.second, 0.0 ) ;
		upper_bound.second += m_minimum_distance_in_cM ;

		// Find an iterator to the lowest SNP in among_these
		// such that recombination position of the (SNP + margin)
		// is greater than or equal to lower_bound
		std::deque< std::size_t >::iterator
			lower_bound_i = std::lower_bound(
				among_these->begin(),
				among_these->end(),
				lower_bound,
				boost::bind(
			 		&RecombinationDistanceProximityTest::compare_to_recombination_position,
					this,	
					_1,
					_2,
					m_margin_in_bp
				)
		) ;

		// Find an iterator to the lowest SNP in among_these
		// such that recombination position of the (SNP - margin)
		// is greater than upper_bound
		std::deque< std::size_t >::iterator
			upper_bound_i = std::upper_bound(
				among_these->begin(),
				among_these->end(),
				upper_bound,
				boost::bind(
			 		&RecombinationDistanceProximityTest::compare_recombination_position_to,
					this,
					_1,
					_2,
					m_margin_in_bp
				)
		) ;
		// We should always have the chosen SNP itself.
		assert( std::distance( lower_bound_i, upper_bound_i ) >= 1 ) ;

		among_these->erase( lower_bound_i, upper_bound_i ) ;
	}

	std::set< std::string > get_attribute_names() const {
		std::set< std::string > result ;
		result.insert( "region_lower_bp" ) ;
		result.insert( "region_upper_bp" ) ;
		result.insert( "region_lower_cm" ) ;
		result.insert( "region_upper_cm" ) ;
		return result ;
	}

	std::map< std::string, genfile::VariantEntry > get_attributes( std::size_t chosen_snp ) const {
		std::pair< genfile::Chromosome, double >
			lower_bound = get_recombination_position( m_snps[ chosen_snp ].get_position() ),
			upper_bound = lower_bound ;

		genfile::GenomePosition lower_bound_bp, upper_bound_bp ;

		lower_bound.second -= m_minimum_distance_in_cM ;
		if( lower_bound.second >= m_genetic_map.get_start_of_map_in_cM( lower_bound.first )) {
		 	lower_bound_bp = m_genetic_map.find_least_physical_position( lower_bound.first, lower_bound.second ) ;
		}
		else {
			lower_bound_bp = m_snps[ chosen_snp ].get_position() ;
		}
		lower_bound_bp.position() -= std::min( m_margin_in_bp, lower_bound_bp.position() ) ;

		upper_bound.second += m_minimum_distance_in_cM ;
		if( upper_bound.second <=  m_genetic_map.get_end_of_map_in_cM( lower_bound.first )) {
			upper_bound_bp = m_genetic_map.find_greatest_physical_position( upper_bound.first, upper_bound.second ) ;
			upper_bound_bp.position() += m_margin_in_bp ;
		}
		else {
			upper_bound_bp.chromosome() = upper_bound.first ;
			upper_bound_bp.position() = std::numeric_limits< genfile::Position >::max() ;
		}

		double lower_bound_cM = get_recombination_position( lower_bound_bp ).second ;
		double upper_bound_cM = get_recombination_position( upper_bound_bp ).second ;

		std::map< std::string, genfile::VariantEntry > result ;
		result[ "region_lower_bp" ] = int( lower_bound_bp.position()  ) ;
		result[ "region_upper_bp" ] = int( upper_bound_bp.position() ) ;
		result[ "region_lower_cm" ] = lower_bound_cM ;
		result[ "region_upper_cm" ] = upper_bound_cM ;
		return result ;
	}

private:
	std::vector< genfile::SNPIdentifyingData > const m_snps ;
	genfile::GeneticMap const& m_genetic_map ;
	double const m_minimum_distance_in_cM ;
	genfile::Position const m_margin_in_bp ;

private:
	std::pair< genfile::Chromosome, double > get_recombination_position(
		genfile::GenomePosition const& position
	) const {
		return std::make_pair(
			position.chromosome(),
			m_genetic_map.find_cM_from_beginning_of_chromosome_at_position( position )
		) ;
	}

	bool compare_recombination_positions( std::size_t a, std::size_t b ) const {
		return get_recombination_position( m_snps[a].get_position() ) < get_recombination_position( m_snps[b].get_position() ) ;
	}

	bool compare_physical_positions( std::size_t a, std::size_t b ) const {
		return m_snps[a].get_position() < m_snps[b].get_position() ;
	}
	
	// Return true if the recombination of the position of the snp with index a
	// plus the margin, specified in kb, has recombination position less than the given one.
	bool compare_to_recombination_position(
		std::size_t a,
		std::pair< genfile::Chromosome, double > const& chromosome_and_offset,
		genfile::Position margin_in_bp
	) const {
		genfile::GenomePosition position = m_snps[a].get_position() ;
		position.position() += margin_in_bp ;
		return get_recombination_position( position ) < chromosome_and_offset ;
	}

	// Return true if the given position is less than that of the position of
	// the snp with index b minus the margin.
	bool compare_recombination_position_to(
		std::pair< genfile::Chromosome, double > const& chromosome_and_offset,
		std::size_t b,
		genfile::Position margin_in_bp
	) const {
		genfile::GenomePosition position = m_snps[b].get_position() ;
		position.position() -= std::min( position.position(), margin_in_bp ) ;
		return chromosome_and_offset < get_recombination_position( position ) ;
	}

} ;

class PhysicalDistanceProximityTest: public ProximityTest
{
public:
	PhysicalDistanceProximityTest(
		std::vector< genfile::SNPIdentifyingData > const& snps,
		double minimum_distance_in_base_pairs
	):
		m_snps( snps ),
		m_minimum_distance_in_base_pairs( minimum_distance_in_base_pairs )
	{}
	
	void prepare( std::deque< std::size_t >* among_these ) const {
		// Sort by chromosome / recombination distance
		std::sort( among_these->begin(), among_these->end(), boost::bind( &PhysicalDistanceProximityTest::compare, this, _1, _2 )) ;
	}

	void remove_snps_too_close_to( std::size_t chosen_snp, std::deque< std::size_t >* among_these ) const {
		genfile::GenomePosition
			lower_bound = m_snps[ chosen_snp ].get_position(),
			upper_bound = lower_bound ;

		if( lower_bound.position() > m_minimum_distance_in_base_pairs ) {
			lower_bound.position() -= m_minimum_distance_in_base_pairs ;
		}
		else {
			lower_bound.position() = 0 ;
		}
		
		upper_bound.position() += m_minimum_distance_in_base_pairs ;

		std::deque< std::size_t >::iterator
			lower_bound_i = std::upper_bound(
				among_these->begin(),
				among_these->end(),
				lower_bound,
				boost::bind(
			 		&PhysicalDistanceProximityTest::compare_position_to_b,
					this,
					_1,
					_2
				)
		) ;

		std::deque< std::size_t >::iterator
			upper_bound_i = std::lower_bound(
				among_these->begin(),
				among_these->end(),
				upper_bound,
				boost::bind(
			 		&PhysicalDistanceProximityTest::compare_a_to_position,
					this,
					_1,
					_2
				)
		) ;

		among_these->erase( lower_bound_i, upper_bound_i ) ;
	}
	
	std::set< std::string > get_attribute_names() const {
		std::set< std::string > result ;
		result.insert( "region_lower_bp" ) ;
		result.insert( "region_upper_bp" ) ;
		return result ;
	}

	std::map< std::string, genfile::VariantEntry > get_attributes( std::size_t chosen_snp ) const {
		genfile::GenomePosition
			lower = m_snps[ chosen_snp ].get_position(),
			upper = m_snps[ chosen_snp ].get_position() ;
			
		lower.position() -= std::min( lower.position(), m_minimum_distance_in_base_pairs ) ;
		upper.position() += m_minimum_distance_in_base_pairs ;

		std::map< std::string, genfile::VariantEntry > result ;
		result[ "region_lower_bp" ] = int( lower.position() ) ;
		result[ "region_upper_bp" ] = int( upper.position() ) ;
		return result ;
	}

	
	std::string display() const {
		return "PhysicalDistanceProximityTest( " + genfile::string_utils::to_string( m_minimum_distance_in_base_pairs ) + "bp )" ;
	}
	
private:
	std::vector< genfile::SNPIdentifyingData > const m_snps ;
	genfile::Position const m_minimum_distance_in_base_pairs ;
	
	bool compare( std::size_t a, std::size_t b ) const {
		return m_snps[ a ].get_position() < m_snps[ b ].get_position() ;
	}
	
	bool compare_a_to_position( std::size_t a, genfile::GenomePosition position ) const {
		return m_snps[ a ].get_position() < position ;
	}

	bool compare_position_to_b( genfile::GenomePosition position, std::size_t b ) const {
		return position < m_snps[ b ].get_position() ;
	}
	
} ;

class InthinneratorApplication: public appcontext::ApplicationContext
{
public:
	InthinneratorApplication( int argc, char** argv ):
		ApplicationContext( globals::program_name, std::auto_ptr< appcontext::OptionProcessor >( new InthinneratorOptionProcessor() ), argc, argv, "-log" )
	{
		process() ;
	}

private:
	enum DistanceType { ePhysical = 0, eRecombination = 1 } ;
	DistanceType m_minimum_distance_type ;
	double m_minimum_distance ;
	ProximityTest::UniquePtr m_proximity_test ;
	SNPPicker::UniquePtr m_snp_picker ;
	bool m_write_incl_list, m_write_excl_list ;

private:
	void process() {
		try {
			unsafe_process() ;
		}
		catch( genfile::FileNotFoundError const& e ) {
			get_ui_context().logger() << "!! Error (" << e.what() << "): the file \"" << e.filespec() << "\" could not be found.\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
		catch( genfile::ChromosomeNotInMapError const& e ) {
			get_ui_context().logger() << "!! Error (" << e.what() << "): chromosome " << e.chromosome() << " was not in the genetic map.\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
		catch( genfile::MalformedInputError const& e ) {
			get_ui_context().logger() << "!! Error (" << e.what() << "): " << e.format_message() << ".\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
		catch( genfile::BadArgumentError const& e ) {
			get_ui_context().logger() << "!! Error (" << e.what() << "): Bad argument ("  << e.arguments() << ") to function " << e.function() << ".\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
		catch( statfile::FileNotOpenedError const& e ) {
			get_ui_context().logger() << "!! Error (" << e.what() << "): The file \"" << e.filename() << "\" could not be opened.\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
		catch( FileNotOpenedError const& e ) {
			get_ui_context().logger() << "!! Error (" << e.what() << "): The file \"" << e.filename() << "\" could not be opened.\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
	}
	
	void unsafe_process() {
		if( !options().check( "-o" ) && ! options().check( "-odb" ) ) {
			get_ui_context().logger() << "You must supply either -o or -odb.\n" ;
			throw appcontext::HaltProgramWithReturnCode( -1 ) ;
		}
		genfile::GeneticMap::UniquePtr map = load_genetic_map() ;
		get_ui_context().logger() << "Loaded: " << map->get_summary() << "\n";

		std::vector< genfile::SNPIdentifyingData > snps = get_list_of_snps( options().get_value< std::string >( "-g" ), *map ) ;
		std::vector< double > recombination_offsets = get_recombination_offsets( snps, *map ) ;
		m_proximity_test = get_proximity_test( snps, *map ) ;
		m_snp_picker = get_snp_picker( snps ) ;

		// Write an inclusion list either if user specified -write-incl-list, or didn't specify any output.
		m_write_incl_list = options().check_if_option_was_supplied( "-write-incl-list" ) || !options().check_if_option_was_supplied( "-write-excl-list" ) ;
		// Write an exclusion list either if user specified -write-excl-list, or didn't specify any output.
		m_write_excl_list = options().check_if_option_was_supplied( "-write-excl-list" ) || !options().check_if_option_was_supplied( "-write-incl-list" ) ;
		
		write_summary( *map, snps, *m_proximity_test, *m_snp_picker ) ;
		
		perform_thinnings( snps, recombination_offsets ) ;
	}
	
	genfile::GeneticMap::UniquePtr load_genetic_map() {
		UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Loading genetic map" ) ;
		genfile::GeneticMap::UniquePtr map(
			new genfile::FromFilesGeneticMap(
				options().get_value< std::string >( "-map" ),
				boost::ref( progress_context )
			)
		) ;
		return map ;
	}

	std::vector< genfile::SNPIdentifyingData > get_list_of_snps( std::string const& filename, genfile::GeneticMap const& map ) const {
		std::vector< genfile::SNPIdentifyingData > filtered_snps ;
		try {
			genfile::SNPDataSource::UniquePtr source ;
			{
				UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Opening genotype files" ) ;
				source = genfile::SNPDataSource::create_chain(
					genfile::wildcard::find_files_by_chromosome( filename ),
					boost::ref( progress_context )
				) ;
			}
			UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Loading SNP list" ) ;
			filtered_snps = genfile::get_list_of_snps_in_source( *source, boost::ref( progress_context ) ) ;
		}
		catch( std::exception const& e ) {
			get_ui_context().logger() << "File \"" << filename << "\" is not a GEN-style file.  Trying text format with columns SNPID rsid chromosome position...\n" ;
			// Not a GEN format file.  Try a different format.
			statfile::BuiltInTypeStatSourceChain::UniquePtr chain(
				statfile::BuiltInTypeStatSourceChain::open(
					genfile::wildcard::find_files_by_chromosome( filename )
				)
			) ;

			UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Loading SNP list" ) ;

			std::string SNPID ;
			std::string rsid ;
			std::string chromosome ;
			genfile::Position position ;
			while(
				(*chain)
					>> SNPID
					>> rsid
					>> chromosome
					>> position
			) {
				filtered_snps.push_back(
					genfile::SNPIdentifyingData(
						SNPID,
						rsid,
						genfile::GenomePosition( chromosome, position ),
						'?',
						'?'
					)
				) ;
				(*chain) >> statfile::ignore_all() ;
				progress_context.notify_progress( chain->number_of_rows_read(), chain->number_of_rows() ) ;
			}
		}
		
		genfile::CommonSNPFilter filter ;
		if( options().has_value( "-excl-rsids" )) {
			std::vector< std::string > filenames = options().get_values< std::string >( "-excl-rsids" ) ;
			for( std::size_t i = 0; i < filenames.size(); ++i ) {
				filter.exclude_snps_in_file( filenames[i], genfile::CommonSNPFilter::RSIDs ) ;
			}
		}
		if( options().has_value( "-incl-rsids" )) {
			filter.exclude_snps_not_in_file( options().get< std::string >( "-incl-rsids" ), genfile::CommonSNPFilter::RSIDs ) ;
		}
		if( options().has_value( "-excl-snpids" )) {
			std::vector< std::string > filenames = options().get_values< std::string >( "-excl-rsids" ) ;
			for( std::size_t i = 0; i < filenames.size(); ++i ) {
				filter.exclude_snps_in_file( filenames[i], genfile::CommonSNPFilter::SNPIDs ) ;
			}
		}
		if( options().has_value( "-incl-snpids" )) {
			filter.exclude_snps_not_in_file( options().get< std::string >( "-incl-snpids" ), genfile::CommonSNPFilter::SNPIDs ) ;
		}

		if( options().check_if_option_was_supplied( "-incl-range" )) {
			std::vector< std::string > specs = genfile::string_utils::split_and_strip_discarding_empty_entries( options().get< std::string >( "-incl-range" ), ",", " \t" ) ;
			for ( std::size_t i = 0; i < specs.size(); ++i ) {
				filter.include_snps_in_range(
					genfile::GenomePositionRange::parse( specs[i] )
				) ;
			}
		}
		
		if( options().check_if_option_was_supplied( "-excl-range" )) {
			std::vector< std::string > specs = genfile::string_utils::split_and_strip_discarding_empty_entries( options().get< std::string >( "-excl-range" ), ",", " \t" ) ;
			for ( std::size_t i = 0; i < specs.size(); ++i ) {
				filter.exclude_snps_in_range(
					genfile::GenomePositionRange::parse( specs[i] )
				) ;
			}
		}
		

		filter.exclude_chromosomes_not_in_set( map.get_chromosomes() ) ;
		
		std::vector< std::size_t > indices_of_included_snps = filter.get_indices_of_filtered_in_snps( filtered_snps ) ;

		filtered_snps = genfile::utility::select_entries( filtered_snps, indices_of_included_snps ) ;

		return filtered_snps ;
	}
	
	std::vector< double > get_recombination_offsets( std::vector< genfile::SNPIdentifyingData > const& snps, genfile::GeneticMap const& map ) const {
		UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Computing recombination positions" ) ;
		std::vector< double > offsets( snps.size() ) ;
		for( std::size_t i = 0; i < snps.size(); ++i ) {
			offsets[i] = map.find_cM_from_beginning_of_chromosome_at_position( snps[i].get_position() ) ;
			progress_context( i+1, snps.size() ) ;
		}
		return offsets ;
	}

	ProximityTest::UniquePtr get_proximity_test(
		std::vector< genfile::SNPIdentifyingData > const& snps,
		genfile::GeneticMap const& genetic_map
	) const {
		ProximityTest::UniquePtr test ;
		std::string distance_spec = options().get< std::string >( "-min-distance" ) ;
		std::vector< std::string > elts = genfile::string_utils::split_and_strip( distance_spec, "+" ) ;

		if( elts.size() < 1 || elts.size() > 2 ) {
			genfile::BadArgumentError( "InthinneratorApplication::get_proximity_test", "-min-distance = \"" + distance_spec + "\"" ) ;
		}
		
		std::vector< std::pair< double, std::string > > parsed_elts( elts.size() ) ;
		for( std::size_t i = 0; i < elts.size(); ++i ) {
			parsed_elts[i] = parse_distance_piece( elts[i] ) ;
		}

		if( parsed_elts.size() == 2 ) {
			// Assume <x>cM+20bp format.
			if( parsed_elts[0].second != "cM" || parsed_elts[1].second != "bp" ) {
				genfile::BadArgumentError( "InthinneratorApplication::get_proximity_test", "-min-distance = \"" + distance_spec + "\"" ) ;
			}
			test.reset(
				new RecombinationDistanceProximityTest(
					snps,
					genetic_map,
					parsed_elts[0].first,
					parsed_elts[1].first
				)
			) ;
		}
		else if( parsed_elts[0].second == "cM" ) {
			test.reset(
				new RecombinationDistanceProximityTest(
					snps,
					genetic_map,
					parsed_elts[0].first,
					0
				)
			) ;
		}
		else if( parsed_elts[0].second == "bp" ) {
			test.reset(
				new PhysicalDistanceProximityTest(
					snps,
					parsed_elts[0].first
				)
			) ;
		}
		else {
			throw genfile::BadArgumentError( "InthinneratorApplication::get_proximity_test", "-min-distance = \"" + distance_spec + "\"" ) ;
		}
		return test ;
	}

	std::pair< double, std::string > parse_distance_piece( std::string const piece ) const {
		std::pair< double, std::string > result ;
		if( piece.compare( piece.size() - 2, 2, "cM" ) == 0 ) {
			result.first = genfile::string_utils::to_repr< double >( piece.substr( 0, piece.size() - 2 )) ;
			result.second = "cM" ;
		}
		else if( piece.compare( piece.size() - 1, 1, "M" ) == 0 ) {
			result.first = 100.0 * genfile::string_utils::to_repr< double >( piece.substr( 0, piece.size() - 1 )) ;
			result.second = "cM" ;
		}
		else if( piece.compare( piece.size() - 2, 2, "bp" ) == 0 ) {
			result.first = genfile::string_utils::to_repr< double >( piece.substr( 0, piece.size() - 2 )) ;
			result.second = "bp" ;
		}
		else if( piece.compare( piece.size() - 2, 2, "kb" ) == 0 ) {
			result.first = 1000.0 * genfile::string_utils::to_repr< double >( piece.substr( 0, piece.size() - 2 )) ;
			result.second = "bp" ;
		}
		else if( piece.compare( piece.size() - 2, 2, "Mb" ) == 0 ) {
			result.first = 1000000.0 * genfile::string_utils::to_repr< double >( piece.substr( 0, piece.size() - 2 )) ;
			result.second = "bp" ;
		}
		return result ;
	}

	SNPPicker::UniquePtr get_snp_picker( std::vector< genfile::SNPIdentifyingData > const& snps ) const {
		SNPPicker::UniquePtr picker ;
		if( options().check_if_option_was_supplied( "-rank" )) {
			std::vector< double > values( snps.size() ) ;
			std::map< genfile::SNPIdentifyingData, double > value_map = load_ranks(
				options().get_value< std::string >( "-rank" ),
				options().get_value< std::string >( "-rank-col" ),
				genfile::string_utils::split_and_strip_discarding_empty_entries( options().get< std::string >( "-missing-code" ), ",", " \t" )
			) ;
			for( std::size_t i = 0; i < values.size(); ++i ) {
				std::map< genfile::SNPIdentifyingData, double >::const_iterator where = value_map.find( snps[i] ) ;
				if( where == value_map.end() ) {
					throw genfile::BadArgumentError( "InthinneratorApplication::get_snp_picker()", "rank file" ) ;
				}
				values[i] = value_map[ snps[i] ] ;
			}
			picker.reset( new HighestValueSNPPicker( values )) ;
		}
		else {
			std::string strategy = options().get< std::string >( "-strategy" ) ;
			if( strategy == "random" ) {
				picker.reset( new RandomSNPPicker ) ;
			}
			else if( strategy == "first" ) {
				picker.reset( new FirstAvailableSNPPicker ) ;
			}
			else {
				throw genfile::BadArgumentError( "InthinneratorApplication::get_snp_picker", "-strategy = \"" + strategy + "\"" ) ;
			}
		}
		return picker ;
	}

	std::map< genfile::SNPIdentifyingData, double > load_ranks( std::string const& filename, std::string const& column_name, std::vector< std::string > const& missing_values ) const {
		statfile::BuiltInTypeStatSourceChain::UniquePtr chain(
			statfile::BuiltInTypeStatSourceChain::open(
				genfile::wildcard::find_files_by_chromosome(
					filename
				)
			)
		) ;

		std::size_t rank_column_index ;
		double sign ;
		
		if( chain->has_column( column_name )) {
			rank_column_index = chain->index_of_column( column_name ) ;
			sign = 1.0 ;
		}
		else if( column_name.size() > 0 && column_name[0] == '-' ) {
			rank_column_index = chain->index_of_column( column_name.substr( 1, column_name.size() ) ) ;
			sign = -1.0 ;
		}
		else {
			throw genfile::BadArgumentError( "load_ranks()", "column_name=\"" + column_name + "\"" ) ;
		}
		
		return load_ranks( *chain, rank_column_index, sign, std::set< std::string >( missing_values.begin(), missing_values.end() ) ) ;
	}

	std::map< genfile::SNPIdentifyingData, double > load_ranks( statfile::BuiltInTypeStatSourceChain& chain, std::size_t rank_column_index, double sign, std::set< std::string > const& missing_values ) const {
		if( rank_column_index < 4 ) {
			throw genfile::BadArgumentError( "InthinneratorApplication::load_ranks()", "rank_column_index" ) ;
		}

		UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Loading ranks" ) ;
		
		std::map< genfile::SNPIdentifyingData, double > result ;
		std::string SNPID ;
		std::string rsid ;
		std::string chromosome ;
		genfile::Position position ;
		std::string rank_string ;
		double rank ;
		while(
			chain
				>> SNPID
				>> rsid
				>> chromosome
				>> position
				>> statfile::ignore( rank_column_index - 4 )
				>> rank_string
		) {
			if( missing_values.find( rank_string ) != missing_values.end() ) {
				rank = std::numeric_limits< double >::quiet_NaN() ;
			}
			else {
				rank = genfile::string_utils::to_repr< double >( rank_string ) ;
			}
			result[ genfile::SNPIdentifyingData( SNPID, rsid, genfile::GenomePosition( chromosome, position ), '?', '?' ) ] = sign * rank ;
			chain >> statfile::ignore_all() ;
			progress_context.notify_progress( chain.number_of_rows_read(), chain.number_of_rows() ) ;
		}
		return result ;
	}

	void write_summary(
		genfile::GeneticMap const& map,
		std::vector< genfile::SNPIdentifyingData > const& snps_in_source,
		ProximityTest const& proximity_test,
		SNPPicker const& snp_picker
	) const {
		get_ui_context().logger() << "Genetic map file(s): " << options().get_value< std::string >( "-map" ) << ".\n" ;
		get_ui_context().logger() << "Genotype file(s): " << options().get_value< std::string >( "-g" ) << ".\n" ;
		get_ui_context().logger() << "Genetic map covers these chromosome(s): " ;
		{
			std::set< genfile::Chromosome > chromosomes = map.get_chromosomes() ;
			for( std::set< genfile::Chromosome >::const_iterator i = chromosomes.begin(); i != chromosomes.end(); ++i ) {
				get_ui_context().logger() << *i << " " ;
			}
			get_ui_context().logger() << "\n" ;
		}
		
		get_ui_context().logger() << "- There are " << snps_in_source.size() << " genotyped SNPs in these chromosomes." ;
		if( snps_in_source.empty() ) {
			get_ui_context().logger() << ".\n" ;
		}
		else {
			get_ui_context().logger() << "  The first few are:\n" ;
			for( std::size_t i = 0; i < std::min( snps_in_source.size(), std::size_t(3) ); ++i ) {
				get_ui_context().logger() << "  " << snps_in_source[i] << "...\n" ;
			}
			get_ui_context().logger() << "  .\n  .\n  .\n" ;
		}
		
		get_ui_context().logger()
			<< "Throwing out SNPs based on proximity test: "
			<< proximity_test.display() << ".\n" ;
		
		get_ui_context().logger()
			<< "Picking SNPs using SNP picker: "
			<< snp_picker.display() << ".\n" ;
		
	}

	void perform_thinnings(
		std::vector< genfile::SNPIdentifyingData > const& snps,
		std::vector< double > const& recombination_offsets
	) const {
		std::size_t const N = options().get_value< std::size_t >( "-N" ) ;
		std::size_t const max_num_picks = options().get_value< std::size_t >( "-max-picks" ) ;
		assert( N > 0 ) ;
		std::size_t const number_of_digits = std::max( std::size_t( std::log10( N ) ), std::size_t( 3u )) ;

		std::string filename_stub ;
		bool db = false ;
		if( options().check( "-o" ) ) {
			filename_stub = options().get< std::string >( "-o" ) ;
		} else if( options().check( "-odb" )) {
			filename_stub = options().get< std::string >( "-odb" ) ;
			db = true ;
		}

		for( std::size_t i = 0; i < N; ++i ) {
			get_ui_context().logger() << "Picking " << (i+1) << " of " << N << "..." ;
			std::set< std::size_t > picked_snps = pick_snps( snps, max_num_picks ) ;
			get_ui_context().logger() << picked_snps.size() << " SNPs picked.\n" ;
			if( db ) {
				write_db( i, snps, recombination_offsets, picked_snps, filename_stub ) ;
			}
			else {
				std::ostringstream filenamestr ;
				filenamestr << filename_stub << "." << std::setw( number_of_digits ) << std::setfill( '0' ) << i ;
				write_output_files( snps, recombination_offsets, picked_snps, filenamestr.str() ) ;
			}
		}
	}

	std::set< std::size_t > pick_snps(
		std::vector< genfile::SNPIdentifyingData > const& snps,
		std::size_t const max_picks
	) const {
		std::deque< std::size_t > remaining_snps(
			boost::counting_iterator< std::size_t >( 0 ),
			boost::counting_iterator< std::size_t >( snps.size() )
		) ;
		
		m_proximity_test->prepare( &remaining_snps ) ;
		m_snp_picker->set_sort_order( remaining_snps ) ;

		std::set< std::size_t > picked_snps ;
		{
			UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Picking SNPs" ) ;
			while( picked_snps.size() < max_picks && remaining_snps.size() > 0 ) {
				std::size_t picked_snp = m_snp_picker->pick( remaining_snps ) ;
				m_proximity_test->remove_snps_too_close_to( picked_snp, &remaining_snps ) ;
				picked_snps.insert( picked_snp ) ;
				progress_context.notify_progress( snps.size() - remaining_snps.size(), snps.size() ) ;
			}
		}
		return picked_snps ;
	}
	
	void write_db(
		std::size_t N,
		std::vector< genfile::SNPIdentifyingData > const& snps,
		std::vector< double > const& recombination_offsets,
		std::set< std::size_t > const& indices_of_snps_to_output,
		std::string const& filename
	) const {
		using genfile::string_utils::to_string ;
		std::vector< std::string > const output_columns = get_output_columns() ;
		snp_summary_component::DBOutputter::UniquePtr
			outputter = snp_summary_component::DBOutputter::create( filename, options().get< std::string >( "-analysis-name" ), options().get_values_as_map() ) ;

		UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Writing \"" + filename + "\"" ) ;
		std::size_t snp_index = 0 ;
		for(
			std::set< std::size_t >::const_iterator i = indices_of_snps_to_output.begin() ;
			i != indices_of_snps_to_output.end() ; 
			++i, ++snp_index
		) {
			for( std::size_t j = 0; j < output_columns.size(); ++j ) {
				std::string const column_name = genfile::string_utils::to_lower( output_columns[j] ) ;
				genfile::VariantEntry value ;
				if( column_name == "cm_from_start_of_chromosome" ) {
					value = recombination_offsets[*i] ;
				}
				else {
					std::map< std::string, genfile::VariantEntry > attributes = m_proximity_test->get_attributes( *i ) ;
					std::map< std::string, genfile::VariantEntry >::const_iterator where = attributes.find( column_name ) ;
					if( where == attributes.end() ) {
						attributes = m_snp_picker->get_attributes( *i ) ;
						where = attributes.find( column_name ) ;
					}
					if( where == attributes.end() ) {
						value = genfile::MissingValue() ;
					}
					else {
						value = where->second ;
					}
				}
				
				outputter->operator()(
					snp_index,
					snps[ snp_index ],
					"inthinnerator",
					"inthinneration " + to_string(N) + "/column_name",
					value
				) ;
			}
			progress_context.notify_progress( snp_index+1, snps.size() ) ;
		}
	}

	void write_output_files(
		std::vector< genfile::SNPIdentifyingData > const& snps,
		std::vector< double > const& recombination_offsets,
		std::set< std::size_t > const& picked_snps,
		std::string const& filename
	) const {
		if( m_write_incl_list ) {
			write_output_file(
				snps,
				recombination_offsets,
				picked_snps,
				filename + ".incl"
			) ;
		}

		if( m_write_excl_list ) {
			std::set< std::size_t > unpicked_snps ;
			std::set_difference(
				boost::counting_iterator< std::size_t >( 0 ),
				boost::counting_iterator< std::size_t >( snps.size() ),
				picked_snps.begin(),
				picked_snps.end(),
				std::inserter( unpicked_snps, unpicked_snps.end() )
			) ;
			write_output_file(
				snps,
				recombination_offsets,
				unpicked_snps,
				filename + ".excl"
			) ;
		}
	}

	void write_output_file(
		std::vector< genfile::SNPIdentifyingData > const& snps,
		std::vector< double > const& recombination_offsets,
		std::set< std::size_t > const& indices_of_snps_to_output,
		std::string const& filename		
	) const {
		std::vector< std::string > const output_columns = get_output_columns() ;

		std::ofstream sink( filename.c_str() ) ;
		if( !sink.is_open() ) {
			throw genfile::ResourceNotOpenedError( filename ) ;
		}

		bool output_headers =
			(( output_columns.size() > 1 )
			|| options().check_if_option_was_supplied( "-headers" ))
			&& !options().check_if_option_was_supplied( "-no-headers" ) ;
		
		if( output_headers ) {
			sink << "# This file written by inthinnerator, " << get_current_time_as_string() << ".\n" ;
			for( std::size_t j = 0; j < output_columns.size(); ++j ) {
				if( j > 0 ) {
					sink << " " ;
				}
				sink << output_columns[j] ;
			}
			sink << "\n" ;
		}

		UIContext::ProgressContext progress_context = get_ui_context().get_progress_context( "Writing \"" + filename + "\"" ) ;
		std::size_t count = 0 ;
		for(
			std::set< std::size_t >::const_iterator i = indices_of_snps_to_output.begin() ;
			i != indices_of_snps_to_output.end() ; 
			++i
		) {
			for( std::size_t j = 0; j < output_columns.size(); ++j ) {
				std::string const column_name = genfile::string_utils::to_lower( output_columns[j] ) ;
				if( j > 0 ) {
					sink << " " ;
				}
				if( column_name == "snpid" ) {
					sink << snps[*i].get_SNPID() ;
				}
				else if( column_name == "rsid" ) {
					sink << snps[*i].get_rsid() ;
				}
				else if( column_name == "chromosome" ) {
					sink << snps[*i].get_position().chromosome() ;
				}
				else if( column_name == "position" ) {
					sink << snps[*i].get_position().position() ;
				}
				else if( column_name == "allele1" ) {
					sink << snps[*i].get_first_allele() ;
				}
				else if( column_name == "allele2" ) {
					sink << snps[*i].get_second_allele() ;
				}
				else if( column_name == "cm_from_start_of_chromosome" ) {
					sink << recombination_offsets[*i] ;
				}
				else {
					std::map< std::string, genfile::VariantEntry > attributes = m_proximity_test->get_attributes( *i ) ;
					std::map< std::string, genfile::VariantEntry >::const_iterator where = attributes.find( column_name ) ;
					if( where == attributes.end() ) {
						attributes = m_snp_picker->get_attributes( *i ) ;
						where = attributes.find( column_name ) ;
					}
					if( where == attributes.end() ) {
						sink << "NA" ;
					}
					else {
						sink << where->second ;
					}
				}
			}
			sink << "\n" ;
			progress_context.notify_progress( ++count, indices_of_snps_to_output.size() ) ;
		}
	}
	
	std::string get_current_time_as_string() const {
		time_t rawtime ;
		struct tm * timeinfo ;
		char buffer[30] ;
		std::time( &rawtime ) ;
		timeinfo = std::localtime( &rawtime ) ;
		std::strftime( buffer, 80, "%Y-%m-%d %H:%M:%S", timeinfo ) ;
		return std::string( buffer ) ;
	}
		
	std::vector< std::string > get_output_columns() const {
		std::set< std::string > admissible_cols ;
		admissible_cols.insert( "snpid" ) ;
		admissible_cols.insert( "rsid" ) ;
		admissible_cols.insert( "chromosome" ) ;
		admissible_cols.insert( "position" ) ;
		admissible_cols.insert( "allele1" ) ;
		admissible_cols.insert( "allele2" ) ;
		admissible_cols.insert( "cm_from_start_of_chromosome" ) ;

		{
			std::set< std::string > const& dynamic_cols = m_proximity_test->get_attribute_names() ;
			for(
				std::set< std::string >::const_iterator i = dynamic_cols.begin();
				i != dynamic_cols.end();
				++i
			) {
				admissible_cols.insert( genfile::string_utils::to_lower( *i ) ) ;
			}
		}

		{
			std::set< std::string > const& dynamic_cols = m_snp_picker->get_attribute_names() ;
			for(
				std::set< std::string >::const_iterator i = dynamic_cols.begin();
				i != dynamic_cols.end();
				++i
			) {
				admissible_cols.insert( genfile::string_utils::to_lower( *i ) ) ;
			}
		}

		std::vector< std::string > output_cols = genfile::string_utils::split_and_strip( options().get< std::string >( "-output-cols" ), ",", " \t") ;
		if( output_cols.size() == 1 && output_cols[0] == "all" ) {
			// construct output cols.  Do it as follows to get them in a particular order.
			output_cols.clear() ;
			output_cols = genfile::string_utils::split_and_strip(
				"SNPID,rsid,chromosome,position,allele1,allele2,cM_from_start_of_chromosome,rank,region_lower_bp,region_upper_bp,region_lower_cM,region_upper_cM",
				","
			) ;
			for( std::vector< std::string >::iterator i = output_cols.begin(); i != output_cols.end(); ) {
				if( admissible_cols.find( genfile::string_utils::to_lower( *i )) == admissible_cols.end() ) {
					output_cols.erase( i ) ;
				} else {
					++i ;
				}
			}
		}

		for( std::size_t j = 0; j < output_cols.size(); ++j ) {
			if( admissible_cols.find( genfile::string_utils::to_lower( output_cols[j] )) == admissible_cols.end() ) {
				get_ui_context().logger() << "!! You have specified the output column \"" << output_cols[j] << "\", which I don't recognise.\n" ;
				get_ui_context().logger() << "   Available output columns are:\n" ;
				for(
					std::set< std::string >::const_iterator i = admissible_cols.begin();
					i != admissible_cols.end();
					++i
				) {
					get_ui_context().logger() << "   \"" << *i << "\".\n" ;
				}
				throw genfile::BadArgumentError(
					"InthinneratorApplication::get_output_columns()",
					"output column=" + genfile::string_utils::to_string( output_cols[j] ) + "\""
				) ;
			}
		}
		return output_cols ;
	}

} ;

int main( int argc, char** argv ) {
	try {
		InthinneratorApplication( argc, argv ) ;
	}
	catch( appcontext::HaltProgramWithReturnCode const& e ) {
		return e.return_code() ;
	}
	return 0 ;
}


	
