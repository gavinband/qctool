#include "genfile/CohortIndividualSource.hpp"
#include "genfile/SampleFilteringCohortIndividualSource.hpp"

namespace genfile {
	
	SampleFilteringCohortIndividualSource::UniquePtr SampleFilteringCohortIndividualSource::create(
		CohortIndividualSource::UniquePtr source,
		std::set< std::size_t > const& indices_of_samples_to_exclude
	) {
		CohortIndividualSource::ConstUniquePtr const_source( source.release() ) ;
		return SampleFilteringCohortIndividualSource::create( const_source, indices_of_samples_to_exclude ) ;
	}

	SampleFilteringCohortIndividualSource::UniquePtr SampleFilteringCohortIndividualSource::create(
		CohortIndividualSource::ConstUniquePtr source,
		std::set< std::size_t > const& indices_of_samples_to_exclude
	) {
		return SampleFilteringCohortIndividualSource::UniquePtr(
			new SampleFilteringCohortIndividualSource(
				source,
				indices_of_samples_to_exclude
			)
		) ;
	}

	SampleFilteringCohortIndividualSource::SampleFilteringCohortIndividualSource(
		CohortIndividualSource::ConstUniquePtr source,
		std::set< std::size_t > const& indices_of_samples_to_exclude
	):
		m_source( source ),
		m_indices_of_samples_to_include( get_indices_of_samples_to_include( *m_source, indices_of_samples_to_exclude ))
	{
	}

	std::size_t SampleFilteringCohortIndividualSource::get_number_of_individuals() const {
		return m_indices_of_samples_to_include.size() ;
	}

	CohortIndividualSource::ColumnSpec SampleFilteringCohortIndividualSource::get_column_spec() const {
		return m_source->get_column_spec() ;
	}
	
	CohortIndividualSource::Entry SampleFilteringCohortIndividualSource::get_entry( std::size_t sample_i, std::string const& column_name ) const {
		assert( sample_i < m_indices_of_samples_to_include.size() ) ;
		return m_source->get_entry( m_indices_of_samples_to_include[ sample_i ], column_name ) ;
	}

	bool SampleFilteringCohortIndividualSource::check_for_column( std::string const& column_name ) const {
		return m_source->check_for_column( column_name ) ;
	}

	CohortIndividualSource const& SampleFilteringCohortIndividualSource::get_parent_source() const {
		return *m_source ;
	}

	CohortIndividualSource const& SampleFilteringCohortIndividualSource::get_base_source() const {
		return m_source->get_base_source() ;
	}

	std::string SampleFilteringCohortIndividualSource::get_source_spec() const {
		return "sample-filtered:" + m_source->get_source_spec() ;
	}

	std::vector< std::size_t > SampleFilteringCohortIndividualSource::get_indices_of_samples_to_include(
		CohortIndividualSource const& source,
		std::set< std::size_t > const& indices_of_samples_to_exclude
	) {
		std::vector< std::size_t > indices_to_include ;
		indices_to_include.reserve( source.get_number_of_individuals() ) ;
		for( std::size_t i = 0; i < source.get_number_of_individuals(); ++i ) {
			if( indices_of_samples_to_exclude.find( i ) == indices_of_samples_to_exclude.end() ) {
				indices_to_include.push_back( i ) ;
			}
		}
		return indices_to_include ;
	}
}
