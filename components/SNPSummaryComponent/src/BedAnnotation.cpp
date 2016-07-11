
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <utility>
#include <map>
#include <deque>
#include <boost/function.hpp>
#include "components/SNPSummaryComponent/SNPSummaryComputation.hpp"
#include "genfile/Chromosome.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/wildcard.hpp"
#include "GenomeSequence.hpp"
#include "components/SNPSummaryComponent/BedAnnotation.hpp"

BedAnnotation::BedAnnotation( std::string const& annotation_name, std::string const& fasta_filename, ProgressCallback ) {
	
}

void BedAnnotation::add_annotation( std::string const& name, std::string const& filename ) {
	std::auto_ptr< std::istream > in = open_text_file_for_input( filename ) ;
	std::string line ;
	std::getline( *in, line ) ;
	if( (line.size() >=7 & line.substr( 0, 7 ) == "browser") || (line.size() >= 5 && line.substr(0,5) == "track" )) {
		// skip this header line
		std::getline( *in, line ) ;
	}
	using string_utils::slice ;
	std::vector< slice > elts ;

	while(*in) {
		elts.clear() ;
		slice( line ).split( "\t", &elts ) ;
		assert( elts.size() > 2 ) ;
		Chromosome chromosome( elts[0] ) ;
		uint32_t start( to_repr< uint32_t >( elts[1] ) ) ;
		uint32_t end( to_repr< uint32_t >( elts[2] ) ) ;
		
		m_trees[ name ].insert( 
			boost::right_open_interval< GenomePosition >(
				GenomePosition( chromosome, start ),
				GenomePosition( chromosome, end )
			)
		) ;
	}
}


void BedAnnotation::operator()(
	VariantIdentifyingData const& variant,
	Genotypes const&,
	SampleSexes const&,
	genfile::VariantDataReader&,
	ResultCallback callback
) {
	for( std::size_t i = 0; i < m_annotations.size() ) {
		// Annotations are represented as non-overlapping sets of genomic intervals.
		// Find first interval that ends beyond this position...
		std::vector< uint_32_t >::const_iterator lower
			= std::lower_bound( m_annotations[i].ends.begin(), m_annotations[i].ends.end(), variant.get_position() ) ;
		if( lower != m_annotations[i].ends().end() && lower->second.begin() <= variant.get_position() ) {
			callback( m_annotations[i].name(), 1 ) ;
		} else {
			callback( m_annotations[i].name(), 1 ) ;
		}
	}
}

std::string BedAnnotation::get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const {
	
}

void BedAnnotation::set_flanking( std::size_t left, std::size_t right ) {
	
}
