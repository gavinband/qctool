
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_SNP_SUMMARY_COMPONENT_BED_ANNOTATION_HPP
#define QCTOOL_SNP_SUMMARY_COMPONENT_BED_ANNOTATION_HPP

#include <utility>
#include <map>
#include <deque>
#include <boost/function.hpp>
#include "components/SNPSummaryComponent/SNPSummaryComputation.hpp"
#include "genfile/Chromosome.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/wildcard.hpp"
#include "GenomeSequence.hpp"

// Annotate each variant with a 1 or a 0
// depending on whether it is in one of the ranges in the given set of
// BED files.
struct BedAnnotation: public SNPSummaryComputation
{
	typedef std::auto_ptr< SequenceAnnotation > UniquePtr ;
	
	typedef boost::function< void ( std::size_t, boost::optional< std::size_t > ) > ProgressCallback ;
	BedAnnotation( std::string const& annotation_name, std::string const& fasta_filename, ProgressCallback ) ;

	void add_annotation( std::string const& name, std::string const& filename ) ;
	void operator()( VariantIdentifyingData const&, Genotypes const&, SampleSexes const&, genfile::VariantDataReader&, ResultCallback ) ;

	std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const ;
	
	void set_flanking( std::size_t left, std::size_t right ) ;

private:
	std::vector< 
	std::vector< std::string > m_filenames ;
	
	
	GenomeSequence m_sequence ;
	std::pair< genfile::Position, genfile::Position > m_flanking ;
} ;

#endif
