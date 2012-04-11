
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_SAMPLEFILTERINGCOHORTINDIVIDUALSOURCE_HPP
#define GENFILE_SAMPLEFILTERINGCOHORTINDIVIDUALSOURCE_HPP

#include <string>
#include <set>
#include <vector>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/CompositeCohortIndividualSource.hpp"

namespace genfile {
	class SampleFilteringCohortIndividualSource: public CompositeCohortIndividualSource
	{
	public:
		typedef std::auto_ptr< SampleFilteringCohortIndividualSource > UniquePtr ;
		typedef std::auto_ptr< SampleFilteringCohortIndividualSource const > ConstUniquePtr ;
		static UniquePtr create(
			CohortIndividualSource::UniquePtr source,
			std::set< std::size_t > const& indices_of_samples_to_exclude
		) ;
		static UniquePtr create(
			CohortIndividualSource::ConstUniquePtr source,
			std::set< std::size_t > const& indices_of_samples_to_exclude
		) ;

	public:
		SampleFilteringCohortIndividualSource(
			CohortIndividualSource::ConstUniquePtr source,
			std::set< std::size_t > const& indices_of_samples_to_exclude
		) ;
	
		std::size_t get_number_of_individuals() const ;
		ColumnSpec get_column_spec() const ;
		Entry get_entry( std::size_t sample_i, std::string const& column_name ) const ;

		bool check_for_column( std::string const& column_name ) const ;

		CohortIndividualSource const& get_parent_source() const ;
		CohortIndividualSource const& get_base_source() const ;

		std::string get_source_spec() const ;
	private:
		CohortIndividualSource::ConstUniquePtr m_source ;
		std::vector< std::size_t > const m_indices_of_samples_to_include ;
		static std::vector< std::size_t > get_indices_of_samples_to_include(
			CohortIndividualSource const& source,
			std::set< std::size_t > const& indices_of_samples_to_exclude
		) ;
	} ;
}

#endif
