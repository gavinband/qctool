
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_SNP_SUMMARY_COMPONENT_ANCESTRAL_ALLELE_ANNOTATION_HPP
#define QCTOOL_SNP_SUMMARY_COMPONENT_ANCESTRAL_ALLELE_ANNOTATION_HPP

#include <map>
#include <deque>
#include <boost/function.hpp>
#include "components/SNPSummaryComponent/SNPSummaryComputation.hpp"
#include "genfile/Chromosome.hpp"
#include "genfile/wildcard.hpp"

struct AncestralAlleleAnnotation: public SNPSummaryComputation
{
	typedef boost::function< void ( std::size_t, boost::optional< std::size_t > ) > ProgressCallback ;
	AncestralAlleleAnnotation( std::string const& fasta_filename, ProgressCallback ) ;
	void operator()( SNPIdentifyingData const&, Genotypes const&, genfile::VariantDataReader&, ResultCallback ) ;

	std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const ;
	
private:
	std::string const m_fasta_filename ;
	std::vector< genfile::wildcard::FilenameMatch > const m_filenames ;
	typedef std::deque< char > ChromosomeSequence ;
	typedef genfile::Chromosome Chromosome ;
	typedef std::pair< std::pair< std::size_t, std::size_t >, ChromosomeSequence > ChromosomeRangeAndSequence ;
	typedef std::map< Chromosome, ChromosomeRangeAndSequence > Sequence ;
	Sequence m_sequence ;

	std::string m_organism ;
	std::string m_build ;
	
	void load_sequence( std::vector< genfile::wildcard::FilenameMatch > const& files, Sequence* sequence, ProgressCallback callback ) ;
	void load_sequence( genfile::wildcard::FilenameMatch const& file, Chromosome* chromosome, ChromosomeSequence* sequence, std::pair< std::size_t, std::size_t >* range ) ;
} ;

#endif
