
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_SNP_SUMMARY_COMPONENT_GENOMESEQUENCE_HPP
#define QCTOOL_SNP_SUMMARY_COMPONENT_GENOMESEQUENCE_HPP

#include <utility>
#include <map>
#include <deque>
#include <string>
#include <boost/function.hpp>
#include <boost/optional.hpp>
#include "components/SNPSummaryComponent/SNPSummaryComputation.hpp"
#include "genfile/Chromosome.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/wildcard.hpp"
#include "genfile/GenomePositionRange.hpp"

struct GenomeSequence {
public:
	typedef std::auto_ptr< GenomeSequence > UniquePtr ;
	typedef std::deque< char > ChromosomeSequence ;
	typedef ChromosomeSequence::const_iterator ConstSequenceIterator ;
	typedef std::pair< ConstSequenceIterator, ConstSequenceIterator > ConstSequenceRange ;
	typedef std::pair< genfile::GenomePositionRange, ConstSequenceRange > PhysicalSequenceRange ;
	typedef genfile::Chromosome Chromosome ;
	typedef std::pair< std::pair< genfile::Position, genfile::Position >, ChromosomeSequence > ChromosomeRangeAndSequence ;
	typedef std::map< Chromosome, ChromosomeRangeAndSequence > SequenceData ;
	typedef genfile::VariantEntry OptionalString ;
	typedef boost::function< void ( std::size_t, boost::optional< std::size_t > ) > ProgressCallback ;
public:
	static UniquePtr create( std::vector< std::string > const& fasta_filenames, ProgressCallback ) ;

public:

	GenomeSequence( std::string const& fasta_filename, ProgressCallback ) ;
	GenomeSequence( std::vector< std::string > const& fasta_filenames, ProgressCallback ) ;

	SequenceData const& sequence() const { return m_data ; }
	std::string const get_spec() const ;
	std::string get_summary( std::string const& prefix, std::size_t column_width = 80 ) const ;
	std::vector< genfile::GenomePositionRange > get_ranges() const ;

	bool has_chromosome( genfile::Chromosome const& chromosome ) const ;
	char get_base( genfile::GenomePosition const& position ) const ;
	void get_sequence( genfile::Chromosome const& chromosome, genfile::Position start, genfile::Position end, std::deque< char >* result ) const ;
	PhysicalSequenceRange get_sequence( genfile::Chromosome const& chromosome, genfile::Position start, genfile::Position end ) const ;

private:
	std::vector< std::string > const m_fasta_filenames ;
	SequenceData m_data ;
	std::map< genfile::Chromosome, std::string > m_identifiers ;
	boost::optional< std::string > m_build ;
	boost::optional< std::string > m_organism ;
	
	void load_sequence( std::vector< std::string > const& filenames, ProgressCallback callback ) ;
	void load_sequence( std::string const& name, std::string const& filename ) ;
} ;

#endif
