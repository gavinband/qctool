
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_DATALINKINGCOHORTINDIVIDUALSOURCE_HPP
#define GENFILE_DATALINKINGCOHORTINDIVIDUALSOURCE_HPP

#include <vector>
#include <string>
#include <map>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/optional.hpp>
#include "genfile/CohortIndividualSource.hpp"

namespace genfile {
	// class DataLinkingCohortIndividualSource
	// This class implements a source which links data using the ID_1 columns
	// between several sources passed to it.
	// The first source, which must be passed to the constructor, is regarded as the 'master' source,
	// and contains the sample IDs which the source representes.
	// Subsequent sources may be passed to the add_source() method and can contain any sample IDs at all
	// (but sample IDs not in the first source will be ignored.)
	// It is legal to have two columns of the same name, but only the first will ever be matched by name.
	class DataLinkingCohortIndividualSource: public CohortIndividualSource
	{
	public:
		typedef std::auto_ptr< DataLinkingCohortIndividualSource > UniquePtr ;

	public:
		DataLinkingCohortIndividualSource( CohortIndividualSource::UniquePtr master ) ;
		~DataLinkingCohortIndividualSource() ;

		void add_source( CohortIndividualSource::UniquePtr source ) ;

		std::size_t get_number_of_individuals() const ;
		ColumnSpec get_column_spec() const ;
		Entry get_entry( std::size_t sample_i, std::string const& column_name ) const ;

		CohortIndividualSource const& get_base_source() const ;
		CohortIndividualSource const& get_parent_source() const ;

		std::string get_source_spec() const ;

	private:
		boost::ptr_vector< CohortIndividualSource > m_sources ;
		ColumnSpec m_column_spec ;
		typedef std::map< std::size_t, std::vector< boost::optional< std::size_t > > > IndexMap ;
		IndexMap m_sample_indices ;
	} ;
}

#endif
