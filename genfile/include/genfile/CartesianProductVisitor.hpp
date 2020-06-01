
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_CARTESIANPRODUCTVISITOR_HPP
#define GENFILE_CARTESIANPRODUCTVISITOR_HPP 1

#include <vector>
#include "boost/optional.hpp"
#include "boost/noncopyable.hpp"
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/MultiSourceVisitor.hpp"

namespace genfile {
	struct CartesianProductVisitor: public MultiSourceVisitor {
	public:
		CartesianProductVisitor( bool lower_triangle = false ) ;

		void add_source(
			std::string const& name,
			genfile::SNPDataSource* source
		) ;

		// Specify that the visitor should skip any configuration
		// containing two variants with the same chromosome and
		// within the given distance
		void set_min_distance( int64_t distance ) ;

		// Return the next configuration via the callback
		// Callback receives:
		// - A 0-1 vector indicating which variants have changed
		// since the last call to step() (or all ones on the first call.)
		// - The list of VariantIdentifyingData structures
		// - The VariantDataReaders.
		bool step(  Callback callback ) ;

		// Return the number of combinations that will be returned,
		// if it can be computed.
		boost::optional< std::size_t > count() const ;
	
	private:
		std::vector< std::string > m_names ;
		bool const m_lower_triangle ;
		int64_t m_min_distance ;
		std::vector< genfile::SNPDataSource* > m_sources ;
		std::vector< genfile::VariantIdentifyingData > m_variants ;
		std::vector< genfile::VariantDataReader::SharedPtr > m_readers ;
		std::vector< int > m_changed ;
	
	private:
	
		bool step_impl() ;
		void reset_source( std::size_t i ) ;

		// utility function to return minimum physical distance between any pair
		// of variants on same chromosome, or maximum +pve value of int64 if all
		// on different chromosomes
		int64_t compute_min_distance(
			std::vector< genfile::VariantIdentifyingData > const& variants
		) const ;
	} ;
}

#endif
