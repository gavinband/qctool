
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_SNP_SUMMARY_COMPONENT_GENETIC_MAP_ANNOTATION_HPP
#define QCTOOL_SNP_SUMMARY_COMPONENT_GENETIC_MAP_ANNOTATION_HPP

#include <memory>
#include "genfile/GeneticMap.hpp"
#include "genfile/wildcard.hpp"
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"
#include "components/SNPSummaryComponent/SNPSummaryComputation.hpp"

struct GeneticMapAnnotation: public SNPSummaryComputation
{
public:
	typedef std::auto_ptr< GeneticMapAnnotation > UniquePtr ;
	typedef boost::function< void ( std::size_t, std::size_t ) > ProgressCallback ;
	static UniquePtr create(
		std::vector< genfile::wildcard::FilenameMatch > const& filenames,
		ProgressCallback = ProgressCallback()
	) ;

public:
	GeneticMapAnnotation(
		std::vector< genfile::wildcard::FilenameMatch > const& filenames,
		ProgressCallback = ProgressCallback()
	) ;
	void operator()( VariantIdentifyingData const&, Genotypes const&, SampleSexes const&, genfile::VariantDataReader&, ResultCallback ) ;
	std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const ;
private:
	genfile::GeneticMap::UniquePtr m_map ;
	std::set< genfile::Chromosome > m_map_chromosomes ;
} ;

#endif
