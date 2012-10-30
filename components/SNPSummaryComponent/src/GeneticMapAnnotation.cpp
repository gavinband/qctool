
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <memory>
#include "genfile/wildcard.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/FromFilesGeneticMap.hpp"
#include "components/SNPSummaryComponent/GeneticMapAnnotation.hpp"

GeneticMapAnnotation::UniquePtr GeneticMapAnnotation::create(
	std::vector< genfile::wildcard::FilenameMatch > const& filenames,
	ProgressCallback callback
) {
	return GeneticMapAnnotation::UniquePtr(
		new GeneticMapAnnotation( filenames, callback )
	) ;
}

GeneticMapAnnotation::GeneticMapAnnotation(
	std::vector< genfile::wildcard::FilenameMatch > const& filenames,
	ProgressCallback callback
):
	m_map(
		new genfile::FromFilesGeneticMap( filenames, callback )
	),
	m_map_chromosomes( m_map->get_chromosomes() )
{}

void GeneticMapAnnotation::operator()(
	SNPIdentifyingData const& snp,
	Genotypes const&,
	SampleSexes const&,
	genfile::VariantDataReader&,
	ResultCallback callback
) {
	if( m_map_chromosomes.find( snp.get_position().chromosome() ) != m_map_chromosomes.end() ) {
		callback(
			"cM_from_start_of_chromosome",
			m_map->find_cM_from_beginning_of_chromosome_at_position( snp.get_position() )
		) ;
	}
}

std::string GeneticMapAnnotation::get_summary( std::string const& prefix, std::size_t column_width ) const {
	std::string result = prefix + "GeneticMapAnnotation: using " + m_map->get_summary() ;
	return result ;
}
