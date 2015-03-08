
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_COHORT_INDIVIDUAL_SOURCE_CHAIN_HPP
#define GENFILE_COHORT_INDIVIDUAL_SOURCE_CHAIN_HPP

#include <memory>
#include <string>
#include <boost/ptr_container/ptr_vector.hpp>
#include "genfile/CohortIndividualSource.hpp"

namespace genfile {
	class CohortIndividualSourceChain: public CohortIndividualSource
	{
	public:
		typedef std::auto_ptr< CohortIndividualSourceChain > UniquePtr ;
	public:
		CohortIndividualSourceChain() ;
		~CohortIndividualSourceChain() ;
		
		void add_source( CohortIndividualSource::UniquePtr source, std::string const& source_name ) ;
		
		std::size_t get_number_of_individuals() const ;
		ColumnSpec get_column_spec() const ;

		Entry get_entry( std::size_t sample_i, std::string const& column_name ) const ;

		std::string get_source_spec() const ;
	private:
		boost::ptr_vector< CohortIndividualSource > m_sources ;
        std::vector< std::string > m_source_names ;
		ColumnSpec m_column_spec ;
	} ;
}

#endif

