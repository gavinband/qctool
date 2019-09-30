
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_SNP_SUMMARY_COMPONENT_BED3_ANNOTATION_HPP
#define QCTOOL_SNP_SUMMARY_COMPONENT_BED3_ANNOTATION_HPP

#include <utility>
#include <map>
#include <deque>
#include <boost/function.hpp>
#include <boost/icl/interval_set.hpp>
#include "components/SNPSummaryComponent/SNPSummaryComputation.hpp"
#include "genfile/Chromosome.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/wildcard.hpp"

namespace stats {
	// Annotate each variant with a set of 1's or 0's according to whether
	// it lies in ranges in one or more BED files (plus or minus a margin).
	// This just uses the first three columns (chrom, start, end) of the BED file.
	// Filenames of BED files are passed in using the add_annotation() function.
	// This class translates between BED style 0-based, half-open coordinates
	// and the 1-based coordinates used in qctool implicitly.
	struct Bed3Annotation: public SNPSummaryComputation
	{
	public:
		typedef std::auto_ptr< Bed3Annotation > UniquePtr ;
		static UniquePtr create() ;

	public:
		Bed3Annotation() ;
		void add_annotation( std::string const& name, std::string const& filename, int left_margin_bp = 0, int right_margin_bp = 0 ) ;
		void operator()( VariantIdentifyingData const&, Genotypes const&, Ploidy const&, genfile::VariantDataReader&, ResultCallback ) ;

		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const ;
	private:
		typedef boost::icl::interval_set< genfile::GenomePosition > Annotation ;
		typedef std::map< std::string, Annotation > AnnotationMap ;
		AnnotationMap m_annotations ;
		std::vector< std::string > m_annotation_names ;
	} ;
}

#endif
