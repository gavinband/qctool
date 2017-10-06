
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <utility>
#include <map>
#include <deque>
#include <boost/icl/interval_set.hpp>
#include <boost/function.hpp>
#include <boost/format.hpp>
#include "components/SNPSummaryComponent/SNPSummaryComputation.hpp"
#include "genfile/Chromosome.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/wildcard.hpp"
#include "genfile/FileUtils.hpp"
#include "components/SNPSummaryComponent/Bed3Annotation.hpp"

Bed3Annotation::UniquePtr Bed3Annotation::create() {
	return UniquePtr( new Bed3Annotation() ) ;
}

Bed3Annotation::Bed3Annotation() {
	
}

void Bed3Annotation::add_annotation( std::string const& name, std::string const& filename, int left_margin_bp, int right_margin_bp ) {
	std::auto_ptr< std::istream > in = genfile::open_text_file_for_input( filename ) ;
	std::string line ;
	std::getline( *in, line ) ;
	if( (line.size() >=7 && line.substr( 0, 7 ) == "browser") || (line.size() >= 5 && line.substr(0,5) == "track" )) {
		// skip this header line
		std::getline( *in, line ) ;
	}
	using genfile::string_utils::slice ;
	using genfile::string_utils::to_repr ;
	std::vector< slice > elts ;

	while(*in) {
		elts.clear() ;
		slice( line ).split( "\t", &elts ) ;
		assert( elts.size() > 2 ) ;
		genfile::Chromosome chromosome( elts[0] ) ;
		uint32_t start( to_repr< uint32_t >( elts[1] ) ) ;
		uint32_t end( to_repr< uint32_t >( elts[2] ) ) ;
		
		// Bed file is 0-based, right-open.
		// We store annotations as 1-based, closed intervals.
		// Therefore add 1 to start coord.
		Annotation::interval_type interval(
			genfile::GenomePosition( chromosome, start + 1 - left_margin_bp ),
			genfile::GenomePosition( chromosome, end + 1 + right_margin_bp )
		) ;
		AnnotationMap::iterator where = m_annotations.find( name ) ;
		if( where == m_annotations.end() ) {
			m_annotation_names.push_back( name ) ;
		}
		m_annotations[ name ].insert( interval ) ;
		std::getline( *in, line ) ;
	}
}

void Bed3Annotation::operator()(
	VariantIdentifyingData const& variant,
	Genotypes const&,
	Ploidy const&,
	genfile::VariantDataReader&,
	ResultCallback callback
) {
	for( std::size_t i = 0; i < m_annotation_names.size(); ++i ) {
		AnnotationMap::const_iterator ai = m_annotations.find( m_annotation_names[i] ) ;
		assert( ai != m_annotations.end() ) ;
		int64_t const result = boost::icl::contains( ai->second, variant.get_position() ) ? 1 : 0 ;
		callback( m_annotation_names[i], result ) ;
	}
}

std::string Bed3Annotation::get_summary( std::string const& prefix, std::size_t column_width ) const {
	std::ostringstream o ;
	o << prefix << boost::format( "Bed3Annotation with %d annotations:\n" ) % m_annotation_names.size() ;
	for( std::size_t i = 0; i < m_annotation_names.size(); ++i ) {
		AnnotationMap::const_iterator ai = m_annotations.find( m_annotation_names[i] ) ;
		assert( ai != m_annotations.end() ) ;
		o << prefix << " - " << ( boost::format( "%s (%d disjoint intervals)\n" ) % m_annotation_names[i] % ai->second.iterative_size() ).str() ;
	}
	return o.str() ;
}
